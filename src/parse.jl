# ============================================================================
# MEME File Parsing
# ============================================================================

"""
    detect_alphabet_from_file(filepath) → (sequence_type::Symbol, num_motifs::Int)

Peek into a MEME file to detect the alphabet type and count motifs.
"""
function detect_alphabet_from_file(filepath::String)::Tuple{Symbol, Int}
    sequence_type = :unknown
    num_motifs = 0

    try
        open(filepath, "r") do f
            line_count = 0
            for line in eachline(f)
                line_count += 1

                if line_count <= 100 && startswith(line, "ALPHABET=")
                    alphabet = strip(split(line, "=")[2])
                    if alphabet == "ACGT"
                        sequence_type = :dna
                    elseif alphabet == "ACGU"
                        sequence_type = :rna
                    elseif length(alphabet) >= 20
                        sequence_type = :protein
                    end
                end

                if startswith(line, "MOTIF ")
                    num_motifs += 1
                end
            end
        end
    catch e
        @warn "Could not read $filepath" exception = e
    end

    return (sequence_type, num_motifs)
end

"""
    extract_species_from_filename(filename) → Vector{String}

Extract species names from filename using Genus_species patterns.
"""
function extract_species_from_filename(filename::String)::Vector{String}
    species = String[]
    base = replace(filename, ".dna_encoded.meme" => "", ".meme" => "")

    species_pattern = r"([A-Z][a-z]+_[a-z]+)"
    for m in eachmatch(species_pattern, base)
        match_str = m.match
        if !any(x -> occursin(x, lowercase(match_str)),
                ["jaspar", "hocomoco", "cisbp", "uniprobe"])
            push!(species, replace(match_str, "_" => " "))
        end
    end

    return species
end

"""
    extract_species_from_path(folder, filename) → Vector{String}

Extract species information from folder and filename context.
"""
function extract_species_from_path(folder::String, filename::String)::Vector{String}
    species = String[]
    folder_lower = lowercase(folder)
    filename_lower = lowercase(filename)

    folder_species_map = Dict(
        "human"   => ["Homo sapiens"],
        "mouse"   => ["Mus musculus"],
        "fly"     => ["Drosophila melanogaster"],
        "worm"    => ["Caenorhabditis elegans"],
        "yeast"   => ["Saccharomyces cerevisiae"],
        "ecoli"   => ["Escherichia coli"],
        "arabd"   => ["Arabidopsis thaliana"],
        "malaria" => ["Plasmodium falciparum"],
    )

    for (key, sp) in folder_species_map
        if occursin(key, folder_lower)
            append!(species, sp)
            break
        end
    end

    append!(species, extract_species_from_filename(filename))

    cue_map = [
        ("vertebrate", "vertebrates"),
        ("fungi",      "fungi"),
        ("insect",     "insects"),
        ("plant",      "plants"),
        ("nematode",   "nematodes"),
        ("urochordate","urochordates"),
        ("diatom",     "diatoms"),
    ]

    for (cue, label) in cue_map
        if occursin(cue, filename_lower)
            push!(species, label)
        end
    end

    if occursin("prokaryote", folder_lower)
        push!(species, "prokaryotes")
    end

    # Deduplicate while preserving order
    seen = Set{String}()
    unique_species = String[]
    for s in species
        sl = lowercase(s)
        if sl ∉ seen
            push!(seen, sl)
            push!(unique_species, s)
        end
    end

    return unique_species
end

# ── Full MEME motif parser ──────────────────────────────────────────────────

"""
    parse_meme_file(filepath; kwargs...) → Vector{ParsedMotif}

Parse all motifs from a single MEME-format file.
"""
function parse_meme_file(filepath::String;
                          source_file::String = "",
                          database_name::String = "",
                          species::Vector{String} = String[],
                          sequence_type_override::String = "",
                          is_dna_encoded::Bool = false)::Vector{ParsedMotif}
    motifs = ParsedMotif[]

    alphabet = ""
    strand   = ""

    # Current motif state
    motif_id   = ""
    motif_name = ""
    width      = 0
    nsites     = 0
    evalue     = 0.0
    url        = ""
    matrix     = Vector{Float64}[]
    in_matrix  = false

    matrix_accum = Float64[]
    alength      = 0

    function flush_matrix!()
        if alength > 0 && !isempty(matrix_accum)
            for i in 1:alength:length(matrix_accum)
                stop = min(i + alength - 1, length(matrix_accum))
                push!(matrix, matrix_accum[i:stop])
            end
        end
        empty!(matrix_accum)
    end

    function flush_motif!()
        flush_matrix!()
        if !isempty(motif_id) && !isempty(matrix)
            seq_type = if !isempty(sequence_type_override)
                sequence_type_override
            elseif alphabet == "ACGU"
                "rna"
            elseif alphabet == "ACGT"
                "dna"
            elseif length(alphabet) >= 20
                "protein"
            else
                "unknown"
            end

            push!(motifs, ParsedMotif(
                motif_id, motif_name, alphabet, seq_type, strand,
                length(matrix), nsites, evalue, url,
                matrix, source_file, database_name, species, is_dna_encoded
            ))
        end
        motif_id   = ""
        motif_name = ""
        width      = 0
        nsites     = 0
        evalue     = 0.0
        url        = ""
        matrix     = Vector{Float64}[]
        empty!(matrix_accum)
        in_matrix  = false
    end

    try
        open(filepath, "r") do f
            for raw_line in eachline(f)
                line = strip(raw_line)

                # Header fields
                if startswith(line, "ALPHABET=")
                    alphabet = strip(split(line, "="; limit=2)[2])
                    continue
                end

                if startswith(line, "strands:")
                    strand = strip(line[length("strands:")+1:end])
                    continue
                end

                # New motif block
                if startswith(line, "MOTIF ")
                    flush_motif!()
                    parts = split(line; limit=3)
                    motif_id   = length(parts) >= 2 ? parts[2] : ""
                    motif_name = length(parts) >= 3 ? parts[3] : ""
                    in_matrix  = false
                    continue
                end

                # Letter-probability matrix header
                if startswith(line, "letter-probability matrix:")
                    flush_matrix!()
                    in_matrix = true
                    tokens = split(line)
                    for (j, tok) in enumerate(tokens)
                        if tok == "alength=" && j < length(tokens)
                            alength = something(tryparse(Int, tokens[j+1]), alength)
                        elseif tok == "w=" && j < length(tokens)
                            width = something(tryparse(Int, tokens[j+1]), width)
                        elseif tok == "nsites=" && j < length(tokens)
                            nsites = something(tryparse(Int, tokens[j+1]), nsites)
                        elseif tok == "E=" && j < length(tokens)
                            evalue = something(tryparse(Float64, tokens[j+1]), evalue)
                        elseif startswith(tok, "alength=") && length(tok) > 8
                            alength = something(tryparse(Int, tok[9:end]), alength)
                        elseif startswith(tok, "w=") && length(tok) > 2
                            width = something(tryparse(Int, tok[3:end]), width)
                        elseif startswith(tok, "nsites=") && length(tok) > 7
                            nsites = something(tryparse(Int, tok[8:end]), nsites)
                        elseif startswith(tok, "E=") && length(tok) > 2
                            evalue = something(tryparse(Float64, tok[3:end]), evalue)
                        end
                    end
                    continue
                end

                # URL line
                if startswith(line, "URL ")
                    url = strip(line[5:end])
                    flush_matrix!()
                    in_matrix = false
                    continue
                end

                # Matrix rows
                if in_matrix
                    if isempty(line)
                        flush_matrix!()
                        in_matrix = false
                        continue
                    end

                    any_num = false
                    for tok in split(line)
                        v = tryparse(Float64, tok)
                        if v !== nothing
                            push!(matrix_accum, v)
                            any_num = true
                        end
                    end
                    if !any_num
                        flush_matrix!()
                        in_matrix = false
                    end
                end
            end

            # Flush last motif
            flush_motif!()
        end
    catch e
        @warn "Error parsing $filepath" exception = e
    end

    return motifs
end
