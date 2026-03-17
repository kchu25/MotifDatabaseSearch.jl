#!/usr/bin/env julia
"""
parse2json.jl — Parse every MEME motif file discovered by the index into a
single JSON file containing all position-weight matrices (PWMs / PFMs).

Output structure (motifs.json):
{
  "metadata": { "created": "...", "total_motifs": N, "total_files": M },
  "motifs": [
    {
      "motif_id":       "MA0004.1",
      "motif_name":     "Arnt",
      "alphabet":       "ACGT",
      "sequence_type":  "dna",
      "strand":         "+ -",
      "width":          6,
      "nsites":         20,
      "evalue":         0.0,
      "url":            "https://...",
      "matrix":         [[0.2, 0.8, 0.0, 0.0], ...],   # rows = positions
      "source_file":    "JASPAR/JASPAR2024_CORE_vertebrates_non-redundant_v2.meme",
      "database_name":  "JASPAR",
      "species":        ["vertebrates"],
      "is_dna_encoded": false
    },
    ...
  ]
}

Usage:
    julia --project=. parse2json.jl                     # defaults
    julia --project=. parse2json.jl motifs.json          # custom output path
"""

include("motif_database_query.jl")
using .MotifDatabaseQuery
using JSON3
using Dates

# ─────────────────────────────────────────────────────────────────────────────
# MEME Parser
# ─────────────────────────────────────────────────────────────────────────────

struct ParsedMotif
    motif_id::String
    motif_name::String
    alphabet::String
    sequence_type::String
    strand::String
    width::Int
    nsites::Int
    evalue::Float64
    url::String
    matrix::Vector{Vector{Float64}}    # rows = positions, cols = alphabet letters
    source_file::String                # relative path within database root
    database_name::String
    species::Vector{String}
    is_dna_encoded::Bool
end

"""
Parse all motifs from a single MEME file.

Returns a Vector{ParsedMotif}.
"""
function parse_meme_file(filepath::String;
                          source_file::String = "",
                          database_name::String = "",
                          species::Vector{String} = String[],
                          sequence_type_override::String = "",
                          is_dna_encoded::Bool = false)::Vector{ParsedMotif}
    motifs = ParsedMotif[]
    
    alphabet = ""
    strand = ""
    
    # Current motif state
    motif_id = ""
    motif_name = ""
    width = 0
    nsites = 0
    evalue = 0.0
    url = ""
    matrix = Vector{Float64}[]
    in_matrix = false
    
    # Accumulator for matrix values (handles line-wrapped numbers)
    matrix_accum = Float64[]
    alength = 0            # alphabet length from header — used to chunk rows

    function flush_matrix!()
        # Turn the flat accumulator into rows of `alength` values each.
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
        # reset per-motif state
        motif_id = ""
        motif_name = ""
        width = 0
        nsites = 0
        evalue = 0.0
        url = ""
        matrix = Vector{Float64}[]
        empty!(matrix_accum)
        in_matrix = false
    end

    try
        open(filepath, "r") do f
            for raw_line in eachline(f)
                line = strip(raw_line)

                # ── header fields ──────────────────────────────
                if startswith(line, "ALPHABET=")
                    alphabet = strip(split(line, "="; limit=2)[2])
                    continue
                end

                if startswith(line, "strands:")
                    strand = strip(line[length("strands:")+1:end])
                    continue
                end

                # ── new motif block ────────────────────────────
                if startswith(line, "MOTIF ")
                    flush_motif!()
                    parts = split(line; limit=3)
                    motif_id   = length(parts) >= 2 ? parts[2] : ""
                    motif_name = length(parts) >= 3 ? parts[3] : ""
                    in_matrix = false
                    continue
                end

                # ── letter-probability matrix header ───────────
                if startswith(line, "letter-probability matrix:")
                    flush_matrix!()           # safety: close any open matrix
                    in_matrix = true
                    # Parse metadata: fields are "key= value" (space between = and value)
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
                        # Also handle joined form: "alength=4", "w=9", etc.
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

                # ── URL line ───────────────────────────────────
                if startswith(line, "URL ")
                    url = strip(line[5:end])
                    flush_matrix!()
                    in_matrix = false
                    continue
                end

                # ── matrix rows (with line-wrapping support) ───
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
                        in_matrix = false   # non-numeric line → end of matrix
                    end
                end
            end

            # flush last motif in file
            flush_motif!()
        end
    catch e
        @warn "Error parsing $filepath" exception = e
    end
    
    return motifs
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function main()
    output_path = length(ARGS) >= 1 ? ARGS[1] : "motifs.json"
    database_root = "./motif_databases.12.27/motif_databases"
    index_file = "motif_index.json"
    
    # ── load or build index ────────────────────────────────────
    println("─"^70)
    println("PARSE MEME → JSON")
    println("─"^70)
    
    index = MotifDatabaseIndex(database_root)
    if isfile(index_file)
        println("Loading file index from $index_file ...")
        load_index!(index, index_file)
    else
        println("Building file index ...")
        build_index!(index, peek_files=true, verbose=true)
        save_index(index, index_file)
    end
    
    # ── parse every file ───────────────────────────────────────
    all_motifs = Dict{String, Any}[]
    n_files = length(index.files)
    
    println("\nParsing $n_files MEME files ...")
    
    t0 = time()
    for (i, mf) in enumerate(index.files)
        if i % 200 == 0
            println("  $i / $n_files files ...")
        end
        
        # relative path for compact storage
        rel = relpath(mf.path, database_root)
        
        parsed = parse_meme_file(mf.path;
            source_file      = rel,
            database_name    = mf.database_name,
            species          = mf.species,
            sequence_type_override = string(mf.sequence_type),
            is_dna_encoded   = mf.is_dna_encoded
        )
        
        for pm in parsed
            push!(all_motifs, Dict{String, Any}(
                "motif_id"       => pm.motif_id,
                "motif_name"     => pm.motif_name,
                "alphabet"       => pm.alphabet,
                "sequence_type"  => pm.sequence_type,
                "strand"         => pm.strand,
                "width"          => pm.width,
                "nsites"         => pm.nsites,
                "evalue"         => pm.evalue,
                "url"            => pm.url,
                "matrix"         => pm.matrix,
                "source_file"    => pm.source_file,
                "database_name"  => pm.database_name,
                "species"        => pm.species,
                "is_dna_encoded" => pm.is_dna_encoded
            ))
        end
    end
    elapsed = round(time() - t0; digits=1)
    
    println("\nParsed $(length(all_motifs)) motifs from $n_files files in $(elapsed)s")
    
    # ── write JSON ─────────────────────────────────────────────
    println("Writing $output_path ...")
    
    output = Dict{String, Any}(
        "metadata" => Dict{String, Any}(
            "created"      => string(now()),
            "total_motifs" => length(all_motifs),
            "total_files"  => n_files,
            "database_root" => database_root,
        ),
        "motifs" => all_motifs
    )
    
    open(output_path, "w") do io
        JSON3.write(io, output)
    end
    
    fsize = round(filesize(output_path) / 1024 / 1024; digits=1)
    println("Done — $(output_path) ($(fsize) MB, $(length(all_motifs)) motifs)")
end

main()
