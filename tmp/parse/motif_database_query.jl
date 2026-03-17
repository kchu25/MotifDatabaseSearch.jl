#!/usr/bin/env julia
"""
Motif Database Query Tool (Julia)

This module provides functionality to:
1. Discover and classify available PWM/PFM databases by species, sequence type (DNA/RNA/Protein),
   and database category
2. Query the database to find relevant MEME files based on search criteria

MEME file format detection:
- DNA motifs: ALPHABET= ACGT
- RNA motifs: ALPHABET= ACGU
- Protein motifs: ALPHABET= ACDEFGHIKLMNPQRSTVWY (20 amino acids)
"""
module MotifDatabaseQuery

using JSON3
using Dates

export MotifFile, MotifDatabaseIndex
export build_index!, query, get_file_paths
export list_species, list_databases, list_sequence_types
export save_index, load_index!, print_summary

# ============================================================================
# Data Structures
# ============================================================================

"""
Represents a MEME motif file with extracted metadata.
"""
mutable struct MotifFile
    path::String
    filename::String
    folder::String              # Parent folder (e.g., JASPAR, CIS-BP_1.02, RNA, etc.)
    sequence_type::Symbol       # :dna, :rna, :protein, :unknown
    species::Vector{String}     # Extracted species
    is_dna_encoded::Bool        # True if RNA file is DNA-encoded (T instead of U)
    database_name::String       # Database name extracted from folder
    num_motifs::Int             # Number of motifs in file
end

# Constructor with defaults
function MotifFile(;
    path::String = "",
    filename::String = "",
    folder::String = "",
    sequence_type::Symbol = :unknown,
    species::Vector{String} = String[],
    is_dna_encoded::Bool = false,
    database_name::String = "",
    num_motifs::Int = 0
)
    MotifFile(path, filename, folder, sequence_type, species, is_dna_encoded, database_name, num_motifs)
end

"""
Index and query motif databases.
"""
mutable struct MotifDatabaseIndex
    database_root::String
    files::Vector{MotifFile}
    species_index::Dict{String, Vector{Int}}        # species -> file indices
    sequence_type_index::Dict{Symbol, Vector{Int}}  # sequence_type -> file indices
    folder_index::Dict{String, Vector{Int}}         # folder -> file indices
    available_species::Set{String}
    available_sequence_types::Set{Symbol}
    available_databases::Set{String}
end

# Constructor
function MotifDatabaseIndex(database_root::String)
    MotifDatabaseIndex(
        database_root,
        MotifFile[],
        Dict{String, Vector{Int}}(),
        Dict{Symbol, Vector{Int}}(),
        Dict{String, Vector{Int}}(),
        Set{String}(),
        Set{Symbol}(),
        Set{String}()
    )
end

# ============================================================================
# File Detection Functions
# ============================================================================

"""
Peek into a MEME file to detect the alphabet type and count motifs.

Returns: (sequence_type::Symbol, num_motifs::Int)
"""
function detect_alphabet_from_file(filepath::String)::Tuple{Symbol, Int}
    sequence_type = :unknown
    num_motifs = 0
    
    try
        open(filepath, "r") do f
            line_count = 0
            for line in eachline(f)
                line_count += 1
                
                # Detect alphabet in first 100 lines
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
                
                # Count motifs
                if startswith(line, "MOTIF ")
                    num_motifs += 1
                end
            end
        end
    catch e
        @warn "Could not read $filepath" exception=e
    end
    
    return (sequence_type, num_motifs)
end

"""
Extract species names from filename.

Common patterns:
- Genus_species.meme (e.g., Homo_sapiens.meme)
- prefix_Genus_species_suffix.meme
"""
function extract_species_from_filename(filename::String)::Vector{String}
    species = String[]
    
    # Remove .meme extension and .dna_encoded suffix
    base = replace(filename, ".dna_encoded.meme" => "", ".meme" => "")
    
    # Common species patterns (Genus_species format)
    # Match patterns like Homo_sapiens, Mus_musculus, etc.
    species_pattern = r"([A-Z][a-z]+_[a-z]+)"
    
    for m in eachmatch(species_pattern, base)
        match_str = m.match
        # Validate it looks like a species name (not database names)
        if !any(x -> occursin(x, lowercase(match_str)), ["jaspar", "hocomoco", "cisbp", "uniprobe"])
            push!(species, replace(match_str, "_" => " "))
        end
    end
    
    return species
end

"""
Extract species information from folder and filename context.
"""
function extract_species_from_path(folder::String, filename::String)::Vector{String}
    species = String[]
    folder_lower = lowercase(folder)
    filename_lower = lowercase(filename)
    
    # Direct species folder mappings
    folder_species_map = Dict(
        "human" => ["Homo sapiens"],
        "mouse" => ["Mus musculus"],
        "fly" => ["Drosophila melanogaster"],
        "worm" => ["Caenorhabditis elegans"],
        "yeast" => ["Saccharomyces cerevisiae"],
        "ecoli" => ["Escherichia coli"],
        "arabd" => ["Arabidopsis thaliana"],
        "malaria" => ["Plasmodium falciparum"],
    )
    
    # Check folder name
    for (key, sp) in folder_species_map
        if occursin(key, folder_lower)
            append!(species, sp)
            break
        end
    end
    
    # Check filename for species patterns
    append!(species, extract_species_from_filename(filename))
    
    # Check for specific database cues
    cue_map = [
        ("vertebrate", "vertebrates"),
        ("fungi", "fungi"),
        ("insect", "insects"),
        ("plant", "plants"),
        ("nematode", "nematodes"),
        ("urochordate", "urochordates"),
        ("diatom", "diatoms"),
    ]
    
    for (cue, label) in cue_map
        if occursin(cue, filename_lower)
            push!(species, label)
        end
    end
    
    if occursin("prokaryote", folder_lower)
        push!(species, "prokaryotes")
    end
    
    # Remove duplicates while preserving order
    seen = Set{String}()
    unique_species = String[]
    for s in species
        s_lower = lowercase(s)
        if s_lower ∉ seen
            push!(seen, s_lower)
            push!(unique_species, s)
        end
    end
    
    return unique_species
end

# ============================================================================
# Index Building
# ============================================================================

"""
Scan the database folder and build an index of all MEME files.

# Arguments
- `index`: MotifDatabaseIndex to populate
- `peek_files`: If true, read file headers to detect alphabet type
- `verbose`: If true, print progress information
"""
function build_index!(index::MotifDatabaseIndex; peek_files::Bool=true, verbose::Bool=true)
    if verbose
        println("Scanning motif databases in: $(index.database_root)")
    end
    
    # Find all .meme files recursively
    meme_files = String[]
    for (root, dirs, files) in walkdir(index.database_root)
        for file in files
            if endswith(file, ".meme")
                push!(meme_files, joinpath(root, file))
            end
        end
    end
    
    if verbose
        println("Found $(length(meme_files)) MEME files")
    end
    
    for (i, filepath) in enumerate(meme_files)
        if verbose && i % 100 == 0
            println("Processing file $i/$(length(meme_files))...")
        end
        
        # Get folder name (parent directory relative to database_root)
        rel_path = relpath(filepath, index.database_root)
        parts = splitpath(rel_path)
        folder = length(parts) > 1 ? joinpath(parts[1:end-1]...) : ""
        
        filename = basename(filepath)
        
        # Detect if it's a DNA-encoded RNA file
        is_dna_encoded = occursin(".dna_encoded.", filename)
        
        # Get database name (first folder component)
        database_name = length(parts) > 1 ? parts[1] : ""
        
        # Detect sequence type
        if peek_files
            sequence_type, num_motifs = detect_alphabet_from_file(filepath)
        else
            # Infer from folder/filename if not peeking
            sequence_type = :unknown
            num_motifs = 0
            
            if occursin("protein", lowercase(folder))
                sequence_type = :protein
            else
                sequence_type = :dna
            end
        end
        
        # Override: files in RNA-oriented database folders are RNA motifs,
        # even if the on-disk alphabet is ACGT (DNA-encoded).
        # The is_dna_encoded flag tracks the encoding distinction.
        db_lower = lowercase(database_name)
        if db_lower in ("cisbp-rna", "mirbase", "rna")
            sequence_type = :rna
        end
        
        # Extract species
        species = extract_species_from_path(folder, filename)
        
        # Create MotifFile entry
        motif_file = MotifFile(
            path = filepath,
            filename = filename,
            folder = folder,
            sequence_type = sequence_type,
            species = species,
            is_dna_encoded = is_dna_encoded,
            database_name = database_name,
            num_motifs = num_motifs
        )
        
        push!(index.files, motif_file)
        file_idx = length(index.files)
        
        # Build indices
        for sp in species
            sp_lower = lowercase(sp)
            if !haskey(index.species_index, sp_lower)
                index.species_index[sp_lower] = Int[]
            end
            push!(index.species_index[sp_lower], file_idx)
        end
        
        if !haskey(index.sequence_type_index, sequence_type)
            index.sequence_type_index[sequence_type] = Int[]
        end
        push!(index.sequence_type_index[sequence_type], file_idx)
        
        if !haskey(index.folder_index, folder)
            index.folder_index[folder] = Int[]
        end
        push!(index.folder_index[folder], file_idx)
        
        # Track available options
        for sp in species
            push!(index.available_species, lowercase(sp))
        end
        push!(index.available_sequence_types, sequence_type)
        push!(index.available_databases, database_name)
    end
    
    if verbose
        println("\nIndex built successfully!")
        print_summary(index)
    end
    
    return index
end

# ============================================================================
# Summary and Display
# ============================================================================

"""
Print a summary of available databases.
"""
function print_summary(index::MotifDatabaseIndex)
    println("\n" * "─"^60)
    println("MOTIF DATABASE SUMMARY")
    println("─"^60)
    
    println("\nTotal MEME files: $(length(index.files))")
    
    println("\nSequence Types Available:")
    for seq_type in sort(collect(index.available_sequence_types), by=string)
        count = length(get(index.sequence_type_index, seq_type, Int[]))
        println("  - $seq_type: $count files")
    end
    
    println("\nDatabase Folders ($(length(index.available_databases))):")
    for db in sort(collect(index.available_databases))
        if !isempty(db)
            file_count = Base.count(f -> f.database_name == db, index.files)
            println("  - $db: $file_count files")
        end
    end
    
    println("\nSpecies/Groups Available ($(length(index.available_species))):")
    # Get species counts
    species_counts = [(sp, length(index.species_index[sp])) for sp in index.available_species]
    sort!(species_counts, by=x -> -x[2])
    
    println("  Top 30 by file count:")
    for (sp, cnt) in species_counts[1:min(30, length(species_counts))]
        println("    - $sp: $cnt files")
    end
    
    if length(species_counts) > 30
        println("  ... and $(length(species_counts) - 30) more species/groups")
    end
    
    println("─"^60)
end

# ============================================================================
# Query Functions
# ============================================================================

"""
Query the database for matching MEME files.

# Arguments
- `index`: MotifDatabaseIndex to query
- `species`: Species name(s) to search for (case-insensitive, partial match)
- `sequence_type`: Filter by :dna, :rna, or :protein
- `database`: Database folder name(s) to search in
- `include_dna_encoded`: If false, exclude DNA-encoded RNA files
- `filename_pattern`: Regex pattern to match against filenames

# Returns
Vector of matching MotifFile objects
"""
function query(index::MotifDatabaseIndex;
               species::Union{Nothing, String, Vector{String}} = nothing,
               sequence_type::Union{Nothing, Symbol} = nothing,
               database::Union{Nothing, String, Vector{String}} = nothing,
               include_dna_encoded::Bool = true,
               filename_pattern::Union{Nothing, String, Regex} = nothing)::Vector{MotifFile}
    
    # Start with all indices
    results = Set(1:length(index.files))
    
    # Filter by sequence type
    if sequence_type !== nothing
        seq_matches = Set(get(index.sequence_type_index, sequence_type, Int[]))
        results = intersect(results, seq_matches)
    end
    
    # Filter by species
    if species !== nothing
        species_list = species isa String ? [species] : species
        
        species_matches = Set{Int}()
        for sp in species_list
            sp_lower = lowercase(sp)
            for (i, f) in enumerate(index.files)
                # Check if any of the file's species match
                if any(s -> occursin(sp_lower, lowercase(s)), f.species)
                    push!(species_matches, i)
                # Also check filename
                elseif occursin(replace(sp_lower, " " => "_"), lowercase(f.filename))
                    push!(species_matches, i)
                elseif occursin(replace(sp_lower, " " => ""), replace(lowercase(f.filename), "_" => ""))
                    push!(species_matches, i)
                end
            end
        end
        
        results = intersect(results, species_matches)
    end
    
    # Filter by database
    if database !== nothing
        db_list = database isa String ? [database] : database
        
        db_matches = Set{Int}()
        for db in db_list
            db_lower = lowercase(db)
            for (i, f) in enumerate(index.files)
                if occursin(db_lower, lowercase(f.database_name)) || occursin(db_lower, lowercase(f.folder))
                    push!(db_matches, i)
                end
            end
        end
        
        results = intersect(results, db_matches)
    end
    
    # Filter DNA-encoded
    if !include_dna_encoded
        results = Set(i for i in results if !index.files[i].is_dna_encoded)
    end
    
    # Filter by filename pattern
    if filename_pattern !== nothing
        pattern = filename_pattern isa Regex ? filename_pattern : Regex(filename_pattern, "i")
        results = Set(i for i in results if occursin(pattern, index.files[i].filename))
    end
    
    return [index.files[i] for i in sort(collect(results))]
end

"""
Convenience method to get just the file paths from a query.

Returns: Vector of absolute file paths matching the query
"""
function get_file_paths(index::MotifDatabaseIndex; kwargs...)::Vector{String}
    results = query(index; kwargs...)
    return [f.path for f in results]
end

# ============================================================================
# List Functions
# ============================================================================

"""Return sorted list of all available species/groups."""
function list_species(index::MotifDatabaseIndex)::Vector{String}
    return sort(collect(index.available_species))
end

"""Return sorted list of all database folders."""
function list_databases(index::MotifDatabaseIndex)::Vector{String}
    return sort([db for db in index.available_databases if !isempty(db)])
end

"""Return list of available sequence types."""
function list_sequence_types(index::MotifDatabaseIndex)::Vector{Symbol}
    return sort(collect(index.available_sequence_types), by=string)
end

# ============================================================================
# Save/Load Index
# ============================================================================

"""Save the index to a JSON file for faster subsequent loading."""
function save_index(index::MotifDatabaseIndex, output_path::String)
    data = Dict(
        "database_root" => index.database_root,
        "files" => [
            Dict(
                "path" => f.path,
                "filename" => f.filename,
                "folder" => f.folder,
                "sequence_type" => string(f.sequence_type),
                "species" => f.species,
                "is_dna_encoded" => f.is_dna_encoded,
                "database_name" => f.database_name,
                "num_motifs" => f.num_motifs
            )
            for f in index.files
        ]
    )
    
    open(output_path, "w") do io
        JSON3.pretty(io, data)
    end
    
    println("Index saved to: $output_path")
end

"""Load a previously saved index."""
function load_index!(index::MotifDatabaseIndex, index_path::String)
    data = JSON3.read(read(index_path, String))
    
    index.database_root = data.database_root
    empty!(index.files)
    empty!(index.species_index)
    empty!(index.sequence_type_index)
    empty!(index.folder_index)
    empty!(index.available_species)
    empty!(index.available_sequence_types)
    empty!(index.available_databases)
    
    for item in data.files
        motif_file = MotifFile(
            path = item.path,
            filename = item.filename,
            folder = item.folder,
            sequence_type = Symbol(item.sequence_type),
            species = collect(String, item.species),
            is_dna_encoded = item.is_dna_encoded,
            database_name = item.database_name,
            num_motifs = item.num_motifs
        )
        
        push!(index.files, motif_file)
        file_idx = length(index.files)
        
        # Rebuild indices
        for sp in motif_file.species
            sp_lower = lowercase(sp)
            if !haskey(index.species_index, sp_lower)
                index.species_index[sp_lower] = Int[]
            end
            push!(index.species_index[sp_lower], file_idx)
        end
        
        if !haskey(index.sequence_type_index, motif_file.sequence_type)
            index.sequence_type_index[motif_file.sequence_type] = Int[]
        end
        push!(index.sequence_type_index[motif_file.sequence_type], file_idx)
        
        if !haskey(index.folder_index, motif_file.folder)
            index.folder_index[motif_file.folder] = Int[]
        end
        push!(index.folder_index[motif_file.folder], file_idx)
        
        for sp in motif_file.species
            push!(index.available_species, lowercase(sp))
        end
        push!(index.available_sequence_types, motif_file.sequence_type)
        push!(index.available_databases, motif_file.database_name)
    end
    
    println("Index loaded: $(length(index.files)) files")
    return index
end

end # module


# ============================================================================
# Main Script (when run directly)
# ============================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    using .MotifDatabaseQuery
    
    # Default database path
    database_root = get(ARGS, 1, "./motif_databases.12.27/motif_databases")
    
    println("─"^70)
    println("MOTIF DATABASE QUERY TOOL")
    println("─"^70)
    
    # Build index
    index = MotifDatabaseIndex(database_root)
    build_index!(index, peek_files=true, verbose=true)
    
    # Example queries
    println("\n" * "="^70)
    println("EXAMPLE QUERIES")
    println("="^70)
    
    # Query 1: Human DNA motifs
    println("\n1. Human DNA Motifs:")
    println("-"^50)
    results = query(index, species="homo sapiens", sequence_type=:dna)
    for f in results[1:min(5, length(results))]
        println("   $(f.filename) ($(f.num_motifs) motifs)")
    end
    if length(results) > 5
        println("   ... and $(length(results) - 5) more")
    end
    
    # Query 2: RNA motifs
    println("\n2. RNA Motifs (first 5):")
    println("-"^50)
    results = query(index, sequence_type=:rna)
    for f in results[1:min(5, length(results))]
        println("   [$(f.database_name)] $(f.filename)")
    end
    println("   Total: $(length(results)) files")
    
    # Query 3: JASPAR vertebrate motifs
    println("\n3. JASPAR Vertebrate Motifs:")
    println("-"^50)
    results = query(index, database="JASPAR", species="vertebrates")
    for f in results
        println("   $(f.filename)")
    end
    
    # Query 4: Protein motifs
    println("\n4. Protein Motifs:")
    println("-"^50)
    results = query(index, sequence_type=:protein)
    for f in results
        println("   $(f.filename) ($(f.num_motifs) motifs)")
    end
    
    # Query 5: Get just file paths
    println("\n5. File Paths for Mouse DNA (HOCOMOCO):")
    println("-"^50)
    paths = get_file_paths(index, species="mus musculus", database="MOUSE")
    for p in paths
        println("   $p")
    end
    
    # Save index for faster future loading
    println("\n" * "="^70)
    save_index(index, "motif_index.json")
end
