# ============================================================================
# Data Structures
# ============================================================================

# ── Parsing / Indexing types ─────────────────────────────────────────────────

"""
    MotifFile

Represents a MEME motif file with extracted metadata (used during indexing).
"""
mutable struct MotifFile
    path::String
    filename::String
    folder::String              # Parent folder (e.g., JASPAR, CIS-BP_1.02, RNA, …)
    sequence_type::Symbol       # :dna, :rna, :protein, :unknown
    species::Vector{String}     # Extracted species
    is_dna_encoded::Bool        # True if RNA file is DNA-encoded (T instead of U)
    database_name::String       # Database name extracted from folder
    num_motifs::Int             # Number of motifs in file
end

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
    MotifFile(path, filename, folder, sequence_type, species,
              is_dna_encoded, database_name, num_motifs)
end

"""
    MotifDatabaseIndex

Index for scanning and querying raw MEME files on disk.
"""
mutable struct MotifDatabaseIndex
    database_root::String
    files::Vector{MotifFile}
    species_index::Dict{String, Vector{Int}}
    sequence_type_index::Dict{Symbol, Vector{Int}}
    folder_index::Dict{String, Vector{Int}}
    available_species::Set{String}
    available_sequence_types::Set{Symbol}
    available_databases::Set{String}
end

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

"""
    ParsedMotif

A single motif parsed from a MEME file (intermediate representation).
"""
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
    matrix::Vector{Vector{Float64}}   # rows = positions, cols = alphabet letters
    source_file::String
    database_name::String
    species::Vector{String}
    is_dna_encoded::Bool
end

# ── User-facing types ────────────────────────────────────────────────────────

"""
    MotifInfo

A single motif entry returned by queries. Contains the PWM matrix and all
associated metadata.

# Fields
- `id`            — motif identifier (e.g. `"MA0004.1"`)
- `name`          — motif name (e.g. `"Arnt"`)
- `alphabet`      — `"ACGT"`, `"ACGU"`, or amino-acid string
- `sequence_type` — `"dna"`, `"rna"`, `"protein"`, `"unknown"`
- `strand`        — e.g. `"+ -"`
- `width`         — number of positions (rows)
- `nsites`        — number of sites used to build the matrix
- `evalue`        — E-value from the MEME file
- `url`           — source URL
- `matrix`        — `Matrix{Float64}` of size `(width × alphabet_length)`
- `source_file`   — relative path within the database root
- `database_name` — top-level folder name (JASPAR, CIS-BP_1.02, …)
- `species`       — `Vector{String}` of species tags
- `is_dna_encoded`— `true` if the file is a DNA-encoded version of an RNA motif
"""
struct MotifInfo
    id::String
    name::String
    alphabet::String
    sequence_type::String       # "dna" | "rna" | "protein" | "unknown"
    strand::String
    width::Int
    nsites::Int
    evalue::Float64
    url::String
    matrix::Matrix{Float64}     # (width × alphabet_length)
    source_file::String
    database_name::String
    species::Vector{String}
    is_dna_encoded::Bool
end

function Base.show(io::IO, p::MotifInfo)
    print(io, "MotifInfo(\"", p.id, "\", \"", p.name, "\", ",
          p.sequence_type, ", ", p.width, "×", size(p.matrix, 2), ")")
end

function Base.show(io::IO, ::MIME"text/plain", p::MotifInfo)
    println(io, "MotifInfo: ", p.id, " — ", p.name)
    println(io, "  type      = ", p.sequence_type, "  alphabet = ", p.alphabet)
    println(io, "  width     = ", p.width, "  nsites = ", p.nsites, "  E = ", p.evalue)
    println(io, "  database  = ", p.database_name)
    println(io, "  species   = ", join(p.species, ", "))
    println(io, "  source    = ", p.source_file)
    if p.is_dna_encoded
        println(io, "  ⚠  DNA-encoded RNA motif")
    end
    ncols = size(p.matrix, 2)
    header = if p.alphabet == "ACGT"
        ["A", "C", "G", "T"]
    elseif p.alphabet == "ACGU"
        ["A", "C", "G", "U"]
    else
        ["col$i" for i in 1:ncols]
    end
    nshow = min(p.width, 10)
    println(io, "  matrix (", p.width, "×", ncols, "):")
    print(io, "    pos  ")
    for h in header
        print(io, lpad(h, 10))
    end
    println(io)
    for r in 1:nshow
        print(io, "    ", lpad(r, 3), "  ")
        for c in 1:ncols
            print(io, lpad(round(p.matrix[r, c]; digits=5), 10))
        end
        println(io)
    end
    if p.width > nshow
        println(io, "    ... ($(p.width - nshow) more rows)")
    end
end

"""
    MotifDB

Container for all parsed motifs with indices for fast querying.
"""
struct MotifDB
    motifs::Vector{MotifInfo}
    by_species::Dict{String, Vector{Int}}
    by_sequence_type::Dict{String, Vector{Int}}
    by_database::Dict{String, Vector{Int}}
    metadata::Dict{String, Any}
end

function Base.show(io::IO, db::MotifDB)
    print(io, "MotifDB(", length(db.motifs), " motifs, ",
          length(db.by_species), " species, ",
          length(db.by_database), " databases)")
end
