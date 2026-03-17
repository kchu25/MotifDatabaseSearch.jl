# ============================================================================
# Data Management — download, decompress, parse, and cache
# ============================================================================
#
# Uses Scratch.jl to store the parsed motifs.json in a package-owned scratch
# space inside the Julia depot.  The scratch dir is:
#   ~/.julia/scratchspaces/<pkg-uuid>/motif_data/
#
# Lifecycle:
#   1. `ensure_database()` — called lazily on first query
#   2. If motifs.json exists in the scratch dir → done
#   3. Otherwise: download .tgz → extract → parse all MEME → write JSON → cleanup

using Downloads: download

# The default URL for the motif databases tarball
const DEFAULT_MOTIF_DB_URL =
    "https://meme-suite.org/meme/meme-software/Databases/motifs/motif_databases.12.27.tgz"

# Name of the cached JSON inside the scratch space
const MOTIFS_JSON = "motifs.json"

# ── Scratch space helpers ────────────────────────────────────────────────────

"""
    data_dir() → String

Return the path to the package scratch space for storing the motif database.
Creates the directory if it does not exist.
"""
function data_dir()
    d = @get_scratch!("motif_data")
    return d
end

"""
    motifs_json_path() → String

Return the full path to the cached `motifs.json`.
"""
motifs_json_path() = joinpath(data_dir(), MOTIFS_JSON)

"""
    database_exists() → Bool

Return `true` if the parsed JSON cache already exists.
"""
database_exists() = isfile(motifs_json_path())

# ── Download and extract ─────────────────────────────────────────────────────

"""
    download_and_extract(url::String, dest_dir::String) → String

Download the tarball from `url`, extract it into `dest_dir`, and return the
path to the extracted motif_databases folder.
"""
function download_and_extract(url::String, dest_dir::String)
    tgz_path = joinpath(dest_dir, "motif_databases.tgz")

    @info "Downloading motif databases from $url …"
    download(url, tgz_path)
    @info "Download complete ($(round(filesize(tgz_path) / 1024 / 1024; digits=1)) MB)"

    @info "Extracting archive …"
    run(`tar -xzf $tgz_path -C $dest_dir`)

    # Remove the tarball to save space
    rm(tgz_path; force=true)

    # Locate the extracted folder — typically motif_databases.X.Y/motif_databases
    # Search for the first directory that contains .meme files
    extracted_root = ""
    for entry in readdir(dest_dir; join=true)
        if isdir(entry) && startswith(basename(entry), "motif_databases")
            # Could be motif_databases.12.27 containing motif_databases/
            inner = joinpath(entry, "motif_databases")
            if isdir(inner)
                extracted_root = inner
            else
                extracted_root = entry
            end
            break
        end
    end

    if isempty(extracted_root) || !isdir(extracted_root)
        error("Could not find extracted motif_databases folder in $dest_dir")
    end

    @info "Extracted to $extracted_root"
    return extracted_root
end

# ── Build JSON from extracted database ───────────────────────────────────────

"""
    build_json(database_root::String, output_path::String)

Scan all MEME files under `database_root`, parse every motif, and write a
single JSON file to `output_path`.
"""
function build_json(database_root::String, output_path::String)
    @info "Building file index …"
    index = MotifDatabaseIndex(database_root)
    build_index!(index; peek_files=true, verbose=true)

    n_files = length(index.files)
    @info "Parsing $n_files MEME files …"

    all_motifs = Dict{String, Any}[]
    t0 = time()

    for (i, mf) in enumerate(index.files)
        if i % 200 == 0
            @info "  $i / $n_files files …"
        end

        rel = relpath(mf.path, database_root)

        parsed = parse_meme_file(mf.path;
            source_file            = rel,
            database_name          = mf.database_name,
            species                = mf.species,
            sequence_type_override = string(mf.sequence_type),
            is_dna_encoded         = mf.is_dna_encoded,
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
                "is_dna_encoded" => pm.is_dna_encoded,
            ))
        end
    end

    elapsed = round(time() - t0; digits=1)
    @info "Parsed $(length(all_motifs)) motifs from $n_files files in $(elapsed)s"

    output = Dict{String, Any}(
        "metadata" => Dict{String, Any}(
            "created"       => string(Dates.now()),
            "total_motifs"  => length(all_motifs),
            "total_files"   => n_files,
            "database_root" => database_root,
        ),
        "motifs" => all_motifs,
    )

    @info "Writing $output_path …"
    open(output_path, "w") do io
        JSON3.write(io, output)
    end
    fsize = round(filesize(output_path) / 1024 / 1024; digits=1)
    @info "Done — $(output_path) ($(fsize) MB, $(length(all_motifs)) motifs)"
end

# ── Top-level entry point ────────────────────────────────────────────────────

"""
    ensure_database(; url = DEFAULT_MOTIF_DB_URL, force = false)

Make sure the parsed `motifs.json` cache exists. If it doesn't (or `force` is
true), download the motif databases tarball, extract it, parse all MEME files,
and write the JSON.

This is called automatically by [`load_motif_db()`](@ref) if no JSON is found.
"""
function ensure_database(; url::String = DEFAULT_MOTIF_DB_URL, force::Bool = false)
    json_path = motifs_json_path()

    if !force && isfile(json_path)
        @info "Motif database already cached at $json_path"
        return json_path
    end

    dd = data_dir()

    # Download and extract into the scratch space
    database_root = download_and_extract(url, dd)

    # Parse everything into a single JSON
    build_json(database_root, json_path)

    # Clean up the raw extracted files to save disk space
    @info "Cleaning up raw MEME files …"
    for entry in readdir(dd; join=true)
        if isdir(entry) && startswith(basename(entry), "motif_databases")
            rm(entry; recursive=true, force=true)
        end
    end

    @info "Motif database ready at $json_path"
    return json_path
end

"""
    clear_cache!()

Delete the cached motif database. The next call to `ensure_database()` or
`load_motif_db()` will re-download and re-parse everything.
"""
function clear_cache!()
    dd = data_dir()
    for f in readdir(dd; join=true)
        rm(f; recursive=true, force=true)
    end
    @info "Cache cleared."
end
