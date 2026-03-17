# ============================================================================
# Query — load the JSON cache and query motifs
# ============================================================================

# Global database handle (loaded lazily)
const _MOTIF_DB = Ref{Union{Nothing, MotifDB}}(nothing)

"""
    load_motif_db(json_path::String) → MotifDB

Load the JSON produced by the parsing step and build query indices.
"""
function load_motif_db(json_path::String)::MotifDB
    @info "Loading $json_path …"
    t0 = time()

    raw = open(json_path, "r") do io
        JSON3.read(io)
    end

    metadata = Dict{String, Any}(
        "created"      => string(get(raw.metadata, :created, "")),
        "total_motifs" => Int(get(raw.metadata, :total_motifs, 0)),
        "total_files"  => Int(get(raw.metadata, :total_files, 0)),
    )

    motifs           = MotifInfo[]
    by_species       = Dict{String, Vector{Int}}()
    by_sequence_type = Dict{String, Vector{Int}}()
    by_database      = Dict{String, Vector{Int}}()

    sizehint!(motifs, length(raw.motifs))

    for m in raw.motifs
        rows  = m.matrix
        nrows = length(rows)
        nrows == 0 && continue
        ncols = length(rows[1])

        mat = Matrix{Float64}(undef, nrows, ncols)
        for r in 1:nrows, c in 1:ncols
            mat[r, c] = Float64(rows[r][c])
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
            push!(get!(Vector{Int}, by_species, lowercase(s)), i)
        end
        # Index by sequence type
        push!(get!(Vector{Int}, by_sequence_type, lowercase(pwm.sequence_type)), i)
        # Index by database
        push!(get!(Vector{Int}, by_database, lowercase(pwm.database_name)), i)
    end

    elapsed = round(time() - t0; digits=1)
    @info "Loaded $(length(motifs)) motifs in $(elapsed)s"

    return MotifDB(motifs, by_species, by_sequence_type, by_database, metadata)
end

# ── Lazy-loading accessor ───────────────────────────────────────────────────

"""
    get_db(; force_reload=false) → MotifDB

Return the global `MotifDB` handle, downloading / parsing the database if
needed and loading it on first access.
"""
function get_db(; force_reload::Bool = false)
    if force_reload || _MOTIF_DB[] === nothing
        json_path = ensure_database()
        _MOTIF_DB[] = load_motif_db(json_path)
    end
    return _MOTIF_DB[]
end

# ── Query ────────────────────────────────────────────────────────────────────

"""
    query(db::MotifDB; kwargs...) → Vector{MotifInfo}

Return all `MotifInfo` entries matching the given criteria (AND-combined).

# Keyword arguments
| argument              | description                                        |
|:-----------------------|:---------------------------------------------------|
| `species`             | species name (case-insensitive substring match)    |
| `sequence_type`       | `"dna"`, `"rna"`, or `"protein"`                   |
| `database`            | database folder name (case-insensitive substring)  |
| `motif_id`            | substring match on motif ID                        |
| `motif_name`          | substring match on motif name                      |
| `min_width`           | minimum number of positions                        |
| `max_width`           | maximum number of positions                        |
| `include_dna_encoded` | include DNA-encoded RNA motifs? (default `true`)   |

# Examples
```julia
db = get_db()
query(db; species="Homo_sapiens", sequence_type="dna", database="JASPAR")
query(db; sequence_type="rna", motif_name="RBFOX")
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

    candidates = nothing  # nothing means "all"

    # Use indices for fast set lookups
    if !isempty(sequence_type)
        key  = lowercase(sequence_type)
        idxs = get(db.by_sequence_type, key, Int[])
        candidates = _intersect(candidates, Set(idxs))
    end

    if !isempty(species)
        spl     = lowercase(species)
        spl_alt = replace(spl, "_" => " ")
        matched = Int[]
        for (k, v) in db.by_species
            if occursin(spl, k) || occursin(spl_alt, k)
                append!(matched, v)
            end
        end
        candidates = _intersect(candidates, Set(matched))
    end

    if !isempty(database)
        dbl     = lowercase(database)
        dbl_alt = replace(dbl, "_" => " ")
        matched = Int[]
        for (k, v) in db.by_database
            if occursin(dbl, k) || occursin(dbl_alt, k)
                append!(matched, v)
            end
        end
        candidates = _intersect(candidates, Set(matched))
    end

    # Linear scan on remaining candidates
    result    = MotifInfo[]
    iter      = candidates === nothing ? (1:length(db.motifs)) : sort(collect(candidates))
    id_lower  = lowercase(motif_id)
    name_lower = lowercase(motif_name)

    for i in iter
        pwm = db.motifs[i]

        !include_dna_encoded && pwm.is_dna_encoded && continue
        !isempty(motif_id)   && !occursin(id_lower, lowercase(pwm.id)) && continue
        !isempty(motif_name) && !occursin(name_lower, lowercase(pwm.name)) && continue
        (pwm.width < min_width || pwm.width > max_width) && continue

        push!(result, pwm)
    end

    return result
end

_intersect(::Nothing, b::Set{Int}) = b
_intersect(a::Set{Int}, b::Set{Int}) = intersect(a, b)

# ── Convenience helpers ──────────────────────────────────────────────────────

"""List all distinct species tags in the database."""
list_species(db::MotifDB) = sort(collect(keys(db.by_species)))

"""List all distinct database names in the database."""
list_databases(db::MotifDB) = sort(collect(keys(db.by_database)))

"""List all distinct sequence types in the database."""
list_sequence_types(db::MotifDB) = sort(collect(keys(db.by_sequence_type)))

"""
    matrix_columns(pwm::MotifInfo)

Get the matrix as a named tuple of column vectors for easier access.
"""
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

"""
    information_content(pwm::MotifInfo) → Vector{Float64}

Compute the information content (bits) at each position.
"""
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

"""
    consensus(pwm::MotifInfo) → String

Return the consensus sequence string.
"""
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
