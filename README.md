# MotifDatabaseSearch

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/MotifDatabaseSearch.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/MotifDatabaseSearch.jl/dev/)
[![Build Status](https://github.com/kchu25/MotifDatabaseSearch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/MotifDatabaseSearch.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/MotifDatabaseSearch.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/MotifDatabaseSearch.jl)

A Julia package for searching and querying motif databases (JASPAR, HOCOMOCO,
CIS-BP, CISBP-RNA, and more) from the
[MEME Suite motif database](https://meme-suite.org/meme/doc/download.html).
This covers DNA motifs for transcription factor binding sites (TFBS), RNA
motifs for RNA-binding proteins (RBPs), and protein motifs. PWMs are stored
locally as a single JSON cache so that subsequent queries are fast without
re-downloading.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/kchu25/MotifDatabaseSearch.jl")
```

---

## Quick Start

```julia
using MotifDatabaseSearch

# Load the database.
# On the FIRST call the motif databases tarball (~500 MB) is downloaded,
# extracted, parsed into a single JSON cache, and the raw files are removed.
# Every subsequent call simply loads the JSON — no download needed.
db = get_db()
# MotifDB(100000 motifs, 312 species, 25 databases)
```

---

## Use Cases

### 1. Browse what's available

```julia
db = get_db()

list_databases(db)       # ["cis-bp_1.02", "eukaryote_pssms", "jaspar", ...]
list_sequence_types(db)  # ["dna", "protein", "rna"]
list_species(db)         # ["arabidopsis thaliana", "homo sapiens", ...]
```

### 2. Get all human DNA motifs from JASPAR

```julia
results = query(db; species="Homo_sapiens", sequence_type="dna", database="JASPAR")
println("Found $(length(results)) motifs")

for pwm in results[1:5]
    println(pwm.id, "  ", pwm.name, "  (width=$(pwm.width))")
end
```

### 3. Search by motif name

```julia
# Find all Sox-family DNA motifs between 8 and 12 positions wide
sox = query(db; motif_name="sox", sequence_type="dna", min_width=8, max_width=12)

for pwm in sox
    println(pwm.id, "  ", pwm.name, "  consensus: ", consensus(pwm))
end
```

### 4. Inspect a PWM

```julia
pwm = results[1]
display(pwm)
# MotifInfo: MA0004.1 — Arnt
#   type      = dna  alphabet = ACGT
#   width     = 6  nsites = 20  E = 0.0
#   database  = JASPAR
#   ...

# Access the probability matrix (width × alphabet_length)
pwm.matrix          # Matrix{Float64}

# Named column vectors for DNA (A, C, G, T)
cols = matrix_columns(pwm)
cols.A              # Vector of A probabilities at each position

# Information content (bits) per position
ic = information_content(pwm)

# Consensus sequence string
consensus(pwm)      # e.e.g "CACGTG"
```

### 5. Get all RNA motifs

```julia
rna = query(db; sequence_type="rna")
println("$(length(rna)) RNA motifs across $(length(unique(m.database_name for m in rna))) databases")
```

### 6. Query mouse motifs from a specific database

```julia
mouse = query(db; species="Mus_musculus", database="HOCOMOCO")
```

### 7. Look up a motif by ID

```julia
hits = query(db; motif_id="MA0002")
```

---

## Cache Management

The motif database is stored in Julia's
[scratch space](https://github.com/JuliaPackaging/Scratch.jl) for this package
(`~/.julia/scratchspaces/<uuid>/motif_data/motifs.json`).

```julia
# Check whether the cache already exists
database_exists()        # true / false

# Show the path to the cached JSON
motifs_json_path()

# Force a re-download and re-parse (e.g. to update to a newer database version)
ensure_database(; force=true)

# Or use a custom URL
ensure_database(; url="https://example.com/custom_motif_databases.tgz")

# Delete the cache entirely (next get_db() call will re-download)
clear_cache!()
```
