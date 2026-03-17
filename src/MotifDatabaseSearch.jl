module MotifDatabaseSearch

using JSON3
using Dates
using Scratch

# ── Types ────────────────────────────────────────────────────────────────────
include("types.jl")

# ── MEME file parsing ────────────────────────────────────────────────────────
include("parse.jl")

# ── Index building (scan raw MEME files) ─────────────────────────────────────
include("index.jl")

# ── Data management (download / decompress / cache) ──────────────────────────
include("data_management.jl")

# ── Query (load JSON, search, helpers) ───────────────────────────────────────
include("query.jl")

# ── Public API ───────────────────────────────────────────────────────────────
export
    # Types
    MotifInfo, MotifDB,
    # Database management
    ensure_database, clear_cache!, database_exists, motifs_json_path,
    # Loading & querying
    get_db, load_motif_db, query,
    # Listing helpers
    list_species, list_databases, list_sequence_types,
    # Motif utilities
    matrix_columns, information_content, consensus

end
