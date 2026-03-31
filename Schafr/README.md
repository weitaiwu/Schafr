# Schafr

**S**ingle-**C**ell **H**igh-throughput **A**nalysis **F**ramework for **R**NA 5' End

An integrated R package for single-cell 5' RNA-seq analysis. Replaces the multi-step shell/Python pipeline with a pure R workflow.

## Pipeline

```
BAM (per cell type)
  │
  ▼  filter_tso()        — Filter reads containing TSO sequence
  │
  ▼  bam_to_ctss()       — Extract 5' end positions + count
  │
  ▼  paraclu()           — Parametric density clustering (pure R)
  │
  ▼  paraclu_cut()       — Filter by length, density ratio, containment
  │
  ▼  overlap_with_peaks() — Retain clusters overlapping ATAC peaks
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

# Optional: for annotation
BiocManager::install(c(
  "ChIPseeker",
  "TxDb.Mmusculus.UCSC.mm10.knownGene",
  "org.Mm.eg.db"
))

# Install Schafr
# devtools::install("path/to/Schafr")
# or
# devtools::install_github("your_username/Schafr")
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

### Step-by-step usage

```r
library(Schafr)

# 1. TSO filter
ga <- filter_tso("cluster_0.bam", tso_sequence = "CTAACGGG", max_mismatch = 2)

# 2. CTSS counting
ctss <- bam_to_ctss(ga)

# 3. Paraclu clustering
clusters <- paraclu(ctss, min_value = 1)

# 4. Filter clusters
clusters_cut <- paraclu_cut(clusters, max_length = 200, min_density_ratio = 2)

# 5. ATAC overlap
final <- overlap_with_peaks(clusters_cut, "ATAC_peaks.narrowPeak")

# 6. (Optional) Annotation
annotated <- annotate_intergenic(final, filter_intergenic = TRUE)
```

## Functions

| Function | Description |
|---|---|
| `filter_tso()` | Filter BAM reads by TSO sequence (Hamming distance) |
| `bam_to_ctss()` | Convert BAM/GAlignments to 5' end counts |
| `paraclu()` | Parametric clustering (pure R implementation) |
| `paraclu_cut()` | Filter clusters by length, density, containment |
| `overlap_with_peaks()` | Retain clusters overlapping peak regions |
| `annotate_intergenic()` | ChIPseeker annotation with Distal Intergenic filter |
| `run_pipeline()` | Run the complete pipeline in one call |

## Parameters

### TSO filtering
- `tso_sequence`: TSO oligo sequence (default: "CTAACGGG")
- `max_mismatch`: Max Hamming distance (default: 2)

### Paraclu
- `min_value`: Minimum total count (default: 1)
- `max_length`: Max cluster length in bp (default: 200)
- `min_density_ratio`: Min max_d/min_d ratio (default: 2)

## License

MIT
