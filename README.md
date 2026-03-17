# MotifDatabaseSearch.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/MotifDatabaseSearch.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/MotifDatabaseSearch.jl/dev/)
[![Build Status](https://github.com/kchu25/MotifDatabaseSearch.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/MotifDatabaseSearch.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/MotifDatabaseSearch.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/MotifDatabaseSearch.jl)

Query PWM/PFM motif databases — JASPAR, HOCOMOCO, CIS-BP, CISBP-RNA, and more — directly from Julia. Covers DNA motifs (TFBS), RNA motifs (RBPs), and protein motifs sourced from the [MEME Suite motif collection](https://meme-suite.org/meme/doc/download.html).

On first use the database is downloaded, parsed into a local JSON cache, and the raw files are removed. All subsequent queries read from the cache — no repeated downloads.

---

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/kchu25/MotifDatabaseSearch.jl")
```

---

## Getting Started

```julia
using MotifDatabaseSearch

# First call: downloads & parses the database (~500 MB tarball, one-time).
# Later calls: loads the local JSON cache in seconds.
db = get_db()
```

---

## Examples

#### Explore what's in the database

```julia
list_databases(db)       # ["cis-bp_1.02", "jaspar", "hocomoco", ...]
list_sequence_types(db)  # ["dna", "rna", "protein"]
list_species(db)         # ["homo sapiens", "mus musculus", ...]
```

#### Query human TF motifs from JASPAR

```julia
results = query(db; species="Homo_sapiens", sequence_type="dna", database="JASPAR")

for pwm in results[1:5]
    println(pwm.id, "  ", pwm.name, "  (width=$(pwm.width))  ", consensus(pwm))
end
```

#### Find Sox-family motifs by name

```julia
sox = query(db; motif_name="sox", sequence_type="dna", min_width=8, max_width=12)
```

#### Query RBP motifs (RNA)

```julia
rbp = query(db; sequence_type="rna")
```

#### Look up a motif by ID

```julia
hits = query(db; motif_id="MA0002")
```

#### Inspect a PWM

```julia
pwm = results[1]
display(pwm)

pwm.matrix                 # Matrix{Float64} — (width × alphabet_length)
cols = matrix_columns(pwm) # named tuple: cols.A, cols.C, cols.G, cols.T
information_content(pwm)   # IC in bits per position
consensus(pwm)             # e.g. "CACGTG"
```

---

## Cache Management

The JSON cache lives in Julia's [scratch space](https://github.com/JuliaPackaging/Scratch.jl) for this package — typically `~/.julia/scratchspaces/<uuid>/motif_data/motifs.json`. It persists across sessions and is removed automatically if the package is uninstalled.

```julia
database_exists()           # check whether the cache is present
motifs_json_path()          # show the cache path

ensure_database(force=true) # force a re-download and re-parse
clear_cache!()              # delete the cache entirely
```

