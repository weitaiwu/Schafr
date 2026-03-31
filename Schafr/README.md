# Schafr

**S**ingle-**C**ell **H**igh-throughput **A**nalysis **F**ramework for **R**NA 5' End

An integrated R package for single-cell 5' RNA-seq analysis. Replaces the multi-step shell/Python pipeline with a pure R workflow. Supports mouse, human, and custom species.

## Pipeline

```
BAM (per cell type)
  │
  ▼  filter_tso()          — Filter reads containing TSO sequence
  │
  ▼  bam_to_ctss()         — Extract 5' end positions + count
  │
  ▼  paraclu()             — Parametric density clustering (pure R)
  │
  ▼  paraclu_cut()         — Filter by length, density ratio, containment
  │
  ▼  overlap_with_peaks()  — Retain clusters overlapping ATAC peaks
  │
  ▼  annotate_intergenic() — (Optional) ChIPseeker annotation
  │
  ▼  Final filtered clusters
```

## Installation

```r
# Install Bioconductor dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "GenomicRanges", "GenomicAlignments", "Rsamtools",
  "IRanges", "Biostrings", "rtracklayer"
))

# Optional: for annotation (install based on your species)
# Mouse
BiocManager::install(c("ChIPseeker", "TxDb.Mmusculus.UCSC.mm10.knownGene", "org.Mm.eg.db"))
# Human
BiocManager::install(c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"))

# Install Schafr from GitHub
devtools::install_github("weitaiwu/Schafr", subdir = "Schafr")
```

## Quick Start

### One-line pipeline

```r
library(Schafr)

results <- run_pipeline(
  bam_files  = c("cluster_0.bam", "cluster_1.bam", "cluster_2.bam"),
  atac_peaks = "ATAC_peaks.narrowPeak",
  output_dir = "schafr_output"
)
```

### With ChIPseeker annotation

```r
results <- run_pipeline(
  bam_files         = list.files("bam_dir", pattern = "\\.bam$", full.names = TRUE),
  atac_peaks        = "ATAC_peaks.narrowPeak",
  annotate          = TRUE,
  filter_intergenic = TRUE,
  output_dir        = "schafr_output"
)
```

### Human data (hg38)

```r
results <- run_pipeline(
  bam_files         = list.files("bam_dir", pattern = "\\.bam$", full.names = TRUE),
  atac_peaks        = "ATAC_peaks.narrowPeak",
  genome            = "hg38",
  annotate          = TRUE,
  filter_intergenic = TRUE,
  output_dir        = "schafr_output"
)
```

### Custom TSO sequence and position

```r
results <- run_pipeline(
  bam_files    = c("sample.bam"),
  atac_peaks   = "peaks.narrowPeak",
  tso_sequence = "AAGCAGTGGTATCAACGCAGAGT",
  tso_start    = 1,
  tso_end      = 23,
  max_mismatch = 3,
  output_dir   = "output"
)
```

### Custom species (e.g., rat)

```r
library(TxDb.Rnorvegicus.UCSC.rn6.refGene)

results <- run_pipeline(
  bam_files         = c("rat_cluster_0.bam"),
  atac_peaks        = "rat_ATAC.narrowPeak",
  annotate          = TRUE,
  txdb              = TxDb.Rnorvegicus.UCSC.rn6.refGene,
  anno_db           = "org.Rn.eg.db",
  output_dir        = "rat_output"
)
```

### Step-by-step usage

```r
library(Schafr)

# 1. TSO filter (custom position)
ga <- filter_tso("cluster_0.bam",
  tso_sequence = "CTAACGGG",
  tso_start = 25, tso_end = 32,
  max_mismatch = 2)

# 2. CTSS counting
ctss <- bam_to_ctss(ga)

# 3. Paraclu clustering
clusters <- paraclu(ctss, min_value = 1)

# 4. Filter clusters
clusters_cut <- paraclu_cut(clusters, max_length = 200, min_density_ratio = 2)

# 5. ATAC overlap
final <- overlap_with_peaks(clusters_cut, "ATAC_peaks.narrowPeak")

# 6. (Optional) Annotation — human hg38
annotated <- annotate_intergenic(final, genome = "hg38", filter_intergenic = TRUE)
```

## Functions

| Function | Description |
|---|---|
| `filter_tso()` | Filter BAM reads by TSO sequence (customizable position & sequence) |
| `bam_to_ctss()` | Convert BAM/GAlignments to 5' end counts |
| `paraclu()` | Parametric clustering (pure R implementation) |
| `paraclu_cut()` | Filter clusters by length, density, containment |
| `overlap_with_peaks()` | Retain clusters overlapping peak regions |
| `annotate_intergenic()` | ChIPseeker annotation with species & Distal Intergenic filter |
| `run_pipeline()` | Run the complete pipeline in one call |

## Parameters

### TSO filtering (`filter_tso`)
| Parameter | Default | Description |
|---|---|---|
| `tso_sequence` | `"CTAACGGG"` | TSO oligo sequence |
| `tso_start` | `25` | Check window start position (1-based) |
| `tso_end` | `32` | Check window end position |
| `max_mismatch` | `2` | Max Hamming distance |
| `check_reverse` | `FALSE` | Also check reverse complement at symmetric end position |

### Paraclu (`paraclu` / `paraclu_cut`)
| Parameter | Default | Description |
|---|---|---|
| `min_value` | `1` | Minimum total count per cluster |
| `max_length` | `200` | Max cluster length (bp) |
| `min_density_ratio` | `2` | Min max_d/min_d ratio |

### Annotation (`annotate_intergenic`)
| Parameter | Default | Description |
|---|---|---|
| `genome` | `"mm10"` | Built-in: `"mm10"`, `"mm39"`, `"hg19"`, `"hg38"` |
| `txdb` | `NULL` | Custom TxDb object (overrides genome) |
| `anno_db` | `NULL` | Custom annotation DB (overrides genome) |
| `filter_intergenic` | `FALSE` | Only return Distal Intergenic clusters |

## Species Support

| Genome | Species | TxDb | annoDb |
|---|---|---|---|
| `mm10` | Mouse | TxDb.Mmusculus.UCSC.mm10.knownGene | org.Mm.eg.db |
| `mm39` | Mouse | TxDb.Mmusculus.UCSC.mm39.knownGene | org.Mm.eg.db |
| `hg19` | Human | TxDb.Hsapiens.UCSC.hg19.knownGene | org.Hs.eg.db |
| `hg38` | Human | TxDb.Hsapiens.UCSC.hg38.knownGene | org.Hs.eg.db |
| Custom | Any | User-provided TxDb | User-provided annoDb |

## License

MIT
