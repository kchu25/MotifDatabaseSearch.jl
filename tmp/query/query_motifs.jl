"""
query_motifs.jl — Load the parsed motifs JSON and query motifs by species,
sequence type, database, etc.

The main export is the `MotifDB` struct, which loads the JSON once, builds
lightweight indices, and provides a `query()` method that returns a
`Vector{MotifInfo}`.

# Quick-start
```julia
include("query_motifs.jl")

db = load_motif_db("motifs.json")          # ~10 s first time
results = query(db; species="Homo_sapiens", sequence_type="dna", database="JASPAR")

for m in results[1:5]
    println(m.id, " — ", m.name, "  (", m.width, " positions)")
    display(m.matrix)
end
```

# MotifInfo struct fields
    id            — motif identifier (e.g. "MA0004.1")
    name          — motif name (e.g. "Arnt")
    alphabet      — "ACGT", "ACGU", or amino-acid string
    sequence_type — "dna", "rna", "protein", "unknown"
    strand        — e.g. "+ -"
    width         — number of positions (rows)
    nsites        — number of sites used to build the matrix
    evalue        — E-value from the MEME file
    url           — source URL
    matrix        — Matrix{Float64} of size (width × alphabet_length)
    source_file   — relative path within the database root
    database_name — top-level folder name (JASPAR, CIS-BP_1.02, …)
    species       — Vector{String} of species tags
    is_dna_encoded — true if the file is a DNA-encoded version of an RNA motif
"""

using JSON3

# ─────────────────────────────────────────────────────────────────────────────
# MotifInfo struct — the main user-facing type
# ─────────────────────────────────────────────────────────────────────────────

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

# ── Pretty-printing ───────────────────────────────────────────────────────

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
    # Show matrix (first 10 rows if large)
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

# ─────────────────────────────────────────────────────────────────────────────
# MotifDB — container with indices for fast querying
# ─────────────────────────────────────────────────────────────────────────────

struct MotifDB
    motifs::Vector{MotifInfo}
    # Indices: key → Set of motif indices (1-based into `motifs`)
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

"""
    load_motif_db(json_path::String) → MotifDB

Load the JSON produced by `parse2json.jl` and build query indices.
"""
function load_motif_db(json_path::String)::MotifDB
    println("Loading $json_path ...")
    t0 = time()

    raw = open(json_path, "r") do io
        JSON3.read(io)
    end

    metadata = Dict{String, Any}(
        "created"      => string(get(raw.metadata, :created, "")),
        "total_motifs" => Int(get(raw.metadata, :total_motifs, 0)),
        "total_files"  => Int(get(raw.metadata, :total_files, 0)),
    )

    motifs = MotifInfo[]
    by_species = Dict{String, Vector{Int}}()
    by_sequence_type = Dict{String, Vector{Int}}()
    by_database = Dict{String, Vector{Int}}()

    sizehint!(motifs, length(raw.motifs))

    for (idx, m) in enumerate(raw.motifs)
        # Convert matrix: Vector{Vector} → Matrix{Float64}
        rows = m.matrix
        nrows = length(rows)
        if nrows == 0
            continue
        end
        ncols = length(rows[1])
        mat = Matrix{Float64}(undef, nrows, ncols)
        for r in 1:nrows
            for c in 1:ncols
                mat[r, c] = Float64(rows[r][c])
            end
        end

        sp = String[string(s) for s in m.species]

        pwm = MotifInfo(
            string(m.motif_id),
            string(m.motif_name),
            string(m.alphabet),
            string(m.sequence_type),
            string(m.strand),
            Int(m.width),
            Int(m.nsites),
            Float64(m.evalue),
            string(m.url),
            mat,
            string(m.source_file),
            string(m.database_name),
            sp,
            Bool(m.is_dna_encoded),
        )
        push!(motifs, pwm)
        i = length(motifs)

        # Index by species
        for s in sp
            sl = lowercase(s)
            v = get!(Vector{Int}, by_species, sl)
            push!(v, i)
        end

        # Index by sequence type
        st = lowercase(pwm.sequence_type)
        v = get!(Vector{Int}, by_sequence_type, st)
        push!(v, i)

        # Index by database
        dl = lowercase(pwm.database_name)
        v = get!(Vector{Int}, by_database, dl)
        push!(v, i)
    end

    elapsed = round(time() - t0; digits=1)
    println("Loaded $(length(motifs)) motifs in $(elapsed)s")

    return MotifDB(motifs, by_species, by_sequence_type, by_database, metadata)
end

# ─────────────────────────────────────────────────────────────────────────────
# Query
# ─────────────────────────────────────────────────────────────────────────────

"""
    query(db::MotifDB; kwargs...) → Vector{MotifInfo}

Return all MotifInfo entries matching the given criteria. All filters are AND-combined.

# Keyword arguments
- `species::String`        — species name (case-insensitive substring match)
- `sequence_type::String`  — "dna", "rna", or "protein"
- `database::String`       — database folder name (case-insensitive substring)
- `motif_id::String`       — substring match on motif ID
- `motif_name::String`     — substring match on motif name
- `min_width::Int`         — minimum number of positions
- `max_width::Int`         — maximum number of positions
- `include_dna_encoded::Bool` — include DNA-encoded RNA motifs? (default: true)

# Examples
```julia
# All human DNA motifs from JASPAR
query(db; species="Homo_sapiens", sequence_type="dna", database="JASPAR")

# All RNA motifs with "RBFOX" in the name
query(db; sequence_type="rna", motif_name="RBFOX")

# All mouse motifs with width 8–12
query(db; species="Mus_musculus", min_width=8, max_width=12)
```
"""
function query(db::MotifDB;
               species::String = "",
               sequence_type::String = "",
               database::String = "",
               motif_id::String = "",
               motif_name::String = "",
               min_width::Int = 0,
               max_width::Int = typemax(Int),
               include_dna_encoded::Bool = true)::Vector{MotifInfo}

    # Start with candidate indices — intersect index lookups for fast filtering
    candidates = nothing  # nothing = "all"

    # ── Use indices for exact-ish lookups ──────────────────────
    if !isempty(sequence_type)
        key = lowercase(sequence_type)
        idxs = get(db.by_sequence_type, key, Int[])
        candidates = _intersect(candidates, Set(idxs))
    end

    if !isempty(species)
        spl = lowercase(species)
        spl_alt = replace(spl, "_" => " ")
        matched = Int[]
        for (k, v) in db.by_species
            if occursin(spl, k) || occursin(spl_alt, k)
                append!(matched, v)
            end
        end
        if isempty(matched)
            @warn "species=\"$species\" matched no motifs. Note: some databases (e.g. JASPAR) do not tag species per-motif — try omitting species= and filtering by database= instead. Use list_species(db) to see available species."
        end
        candidates = _intersect(candidates, Set(matched))
    end

    if !isempty(database)
        dbl = lowercase(database)
        dbl_alt = replace(dbl, "_" => " ")
        matched = Int[]
        for (k, v) in db.by_database
            if occursin(dbl, k) || occursin(dbl_alt, k)
                append!(matched, v)
            end
        end
        if isempty(matched)
            @warn "database=\"$database\" matched no known database. Use list_databases(db) to see available databases."
        end
        candidates = _intersect(candidates, Set(matched))
    end

    # ── Linear scan on remaining candidates ────────────────────
    result = MotifInfo[]
    iter = candidates === nothing ? (1:length(db.motifs)) : sort(collect(candidates))

    id_lower = lowercase(motif_id)
    name_lower = lowercase(motif_name)

    for i in iter
        pwm = db.motifs[i]

        if !include_dna_encoded && pwm.is_dna_encoded
            continue
        end
        if !isempty(motif_id) && !occursin(id_lower, lowercase(pwm.id))
            continue
        end
        if !isempty(motif_name) && !occursin(name_lower, lowercase(pwm.name))
            continue
        end
        if pwm.width < min_width || pwm.width > max_width
            continue
        end

        push!(result, pwm)
    end

    return result
end

function _intersect(a::Nothing, b::Set{Int})
    return b
end
function _intersect(a::Set{Int}, b::Set{Int})
    return intersect(a, b)
end

# ─────────────────────────────────────────────────────────────────────────────
# Convenience helpers
# ─────────────────────────────────────────────────────────────────────────────

"""List all distinct species in the database."""
list_species(db::MotifDB) = sort(collect(keys(db.by_species)))

"""List all distinct databases in the database."""
list_databases(db::MotifDB) = sort(collect(keys(db.by_database)))

"""List all distinct sequence types in the database."""
list_sequence_types(db::MotifDB) = sort(collect(keys(db.by_sequence_type)))

"""Get the matrix as a named tuple of column vectors for easier access."""
function matrix_columns(pwm::MotifInfo)
    if pwm.alphabet == "ACGT"
        return (A=pwm.matrix[:, 1], C=pwm.matrix[:, 2],
                G=pwm.matrix[:, 3], T=pwm.matrix[:, 4])
    elseif pwm.alphabet == "ACGU"
        return (A=pwm.matrix[:, 1], C=pwm.matrix[:, 2],
                G=pwm.matrix[:, 3], U=pwm.matrix[:, 4])
    else
        return Tuple(pwm.matrix[:, c] for c in 1:size(pwm.matrix, 2))
    end
end

"""Compute the information content (bits) at each position."""
function information_content(pwm::MotifInfo)
    npos, nalph = size(pwm.matrix)
    max_bits = log2(nalph)
    ic = Vector{Float64}(undef, npos)
    for r in 1:npos
        H = 0.0
        for c in 1:nalph
            p = pwm.matrix[r, c]
            if p > 0
                H -= p * log2(p)
            end
        end
        ic[r] = max_bits - H
    end
    return ic
end

"""Return the consensus sequence string."""
function consensus(pwm::MotifInfo)
    letters = collect(pwm.alphabet)
    npos = size(pwm.matrix, 1)
    buf = IOBuffer()
    for r in 1:npos
        _, best = findmax(pwm.matrix[r, :])
        if best <= length(letters)
            write(buf, letters[best])
        else
            write(buf, '?')
        end
    end
    return String(take!(buf))
end

# ─────────────────────────────────────────────────────────────────────────────
# Main — demo usage when run directly
# ─────────────────────────────────────────────────────────────────────────────

function main()
    json_path = length(ARGS) >= 1 ? ARGS[1] : "motifs.json"

    db = load_motif_db(json_path)
    println()
    println(repeat('─', 70))
    println(db)
    println(repeat('─', 70))

    # Example queries
    println("\n▶ Human DNA motifs (from HUMAN database):")
    results = query(db; species="homo_sapiens", sequence_type="dna")
    println("  Found $(length(results)) motifs")
    for pwm in results[1:min(3, length(results))]
        println("  • $(pwm.id) — $(pwm.name)  [$(pwm.width) pos]  consensus: $(consensus(pwm))")
    end

    println("\n▶ All JASPAR DNA motifs:")
    results = query(db; sequence_type="dna", database="jaspar")
    println("  Found $(length(results)) motifs")
    for pwm in results[1:min(3, length(results))]
        println("  • $(pwm.id) — $(pwm.name)  [$(pwm.width) pos]  consensus: $(consensus(pwm))")
    end

    println("\n▶ All RNA motifs for Homo sapiens:")
    results = query(db; species="homo_sapiens", sequence_type="rna")
    println("  Found $(length(results)) motifs")
    for pwm in results[1:min(3, length(results))]
        println("  • $(pwm.id) — $(pwm.name)  [$(pwm.width) pos]")
    end

    println("\n▶ All E. coli motifs:")
    results = query(db; species="escherichia_coli")
    println("  Found $(length(results)) motifs")

    println("\n▶ Motifs with 'Sox' in name (DNA, width 8–12):")
    results = query(db; motif_name="sox", sequence_type="dna", min_width=8, max_width=12)
    println("  Found $(length(results)) motifs")
    for pwm in results[1:min(3, length(results))]
        println("  • $(pwm.id) — $(pwm.name)  [$(pwm.width) pos]  consensus: $(consensus(pwm))")
    end

    # Show a detailed view of one motif
    if !isempty(results)
        println("\n", repeat('─', 70))
        println("Detailed view of first Sox match:")
        println(repeat('─', 70))
        display(results[1])
        println()
        ic = information_content(results[1])
        println("  IC per position: ", round.(ic; digits=3))
        cols = matrix_columns(results[1])
        println("  Column A: ", round.(cols.A; digits=4))
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
