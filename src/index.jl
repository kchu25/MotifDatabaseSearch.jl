# ============================================================================
# Index Building — scan extracted MEME database folder
# ============================================================================

"""
    build_index!(index; peek_files=true, verbose=true)

Scan the database folder and populate the `MotifDatabaseIndex` with all
discovered MEME files and their metadata.
"""
function build_index!(index::MotifDatabaseIndex;
                      peek_files::Bool = true,
                      verbose::Bool = true)
    if verbose
        @info "Scanning motif databases in: $(index.database_root)"
    end

    # Find all .meme files recursively
    meme_files = String[]
    for (root, _, files) in walkdir(index.database_root)
        for file in files
            if endswith(file, ".meme")
                push!(meme_files, joinpath(root, file))
            end
        end
    end

    if verbose
        @info "Found $(length(meme_files)) MEME files"
    end

    for (i, filepath) in enumerate(meme_files)
        if verbose && i % 100 == 0
            @info "Processing file $i / $(length(meme_files)) …"
        end

        rel_path = relpath(filepath, index.database_root)
        parts    = splitpath(rel_path)
        folder   = length(parts) > 1 ? joinpath(parts[1:end-1]...) : ""
        filename = basename(filepath)

        is_dna_encoded = occursin(".dna_encoded.", filename)
        database_name  = length(parts) > 1 ? parts[1] : ""

        if peek_files
            sequence_type, num_motifs = detect_alphabet_from_file(filepath)
        else
            sequence_type = :unknown
            num_motifs    = 0
            if occursin("protein", lowercase(folder))
                sequence_type = :protein
            else
                sequence_type = :dna
            end
        end

        # Override for RNA-oriented databases
        db_lower = lowercase(database_name)
        if db_lower in ("cisbp-rna", "mirbase", "rna")
            sequence_type = :rna
        end

        species = extract_species_from_path(folder, filename)

        motif_file = MotifFile(;
            path           = filepath,
            filename       = filename,
            folder         = folder,
            sequence_type  = sequence_type,
            species        = species,
            is_dna_encoded = is_dna_encoded,
            database_name  = database_name,
            num_motifs     = num_motifs,
        )

        push!(index.files, motif_file)
        file_idx = length(index.files)

        # Build indices
        for sp in species
            sp_lower = lowercase(sp)
            push!(get!(Vector{Int}, index.species_index, sp_lower), file_idx)
        end

        push!(get!(Vector{Int}, index.sequence_type_index, sequence_type), file_idx)
        push!(get!(Vector{Int}, index.folder_index, folder), file_idx)

        for sp in species
            push!(index.available_species, lowercase(sp))
        end
        push!(index.available_sequence_types, sequence_type)
        push!(index.available_databases, database_name)
    end

    if verbose
        @info "Index built: $(length(index.files)) files, " *
              "$(length(index.available_databases)) databases, " *
              "$(length(index.available_species)) species tags"
    end

    return index
end
