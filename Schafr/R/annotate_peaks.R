#' Annotate Clusters and Filter Distal Intergenic Regions
#'
#' Uses ChIPseeker to annotate paraclu clusters and optionally filters to
#' retain only those in Distal Intergenic regions. Supports mouse (mm10/mm39),
#' human (hg19/hg38), and custom species via user-provided TxDb/annoDb.
#'
#' @param clusters A data.frame with at least columns: chr, start, end, strand.
#' @param genome Genome identifier. Built-in: "mm10", "mm39", "hg19", "hg38".
#'   Default: "mm10". Ignored if txdb and anno_db are both provided.
#' @param txdb A TxDb object. If NULL, auto-resolved from \code{genome}.
#' @param anno_db Annotation DB name. If NULL, auto-resolved from \code{genome}.
#' @param tss_region Numeric vector of length 2. Default: c(-3000, 3000).
#' @param filter_intergenic Logical; only return Distal Intergenic? Default: FALSE.
#'
#' @details
#' Built-in genome support:
#' \itemize{
#'   \item \strong{mm10}: TxDb.Mmusculus.UCSC.mm10.knownGene + org.Mm.eg.db
#'   \item \strong{mm39}: TxDb.Mmusculus.UCSC.mm39.knownGene + org.Mm.eg.db
#'   \item \strong{hg19}: TxDb.Hsapiens.UCSC.hg19.knownGene + org.Hs.eg.db
#'   \item \strong{hg38}: TxDb.Hsapiens.UCSC.hg38.knownGene + org.Hs.eg.db
#' }
#' For other species, supply your own \code{txdb} and \code{anno_db}.
#'
#' @return A data.frame with annotation columns. If \code{filter_intergenic = TRUE},
#'   only Distal Intergenic clusters.
#'
#' @examples
#' \dontrun{
#' # Mouse mm10 (default)
#' annotate_intergenic(clusters)
#'
#' # Human hg38
#' annotate_intergenic(clusters, genome = "hg38")
#'
#' # Custom species (e.g., rat)
#' annotate_intergenic(clusters,
#'   txdb = TxDb.Rnorvegicus.UCSC.rn6.refGene,
#'   anno_db = "org.Rn.eg.db")
#' }
#'
#' @export
annotate_intergenic <- function(clusters,
                                genome = "mm10",
                                txdb = NULL,
                                anno_db = NULL,
                                tss_region = c(-3000, 3000),
                                filter_intergenic = FALSE) {

  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    stop("Package 'ChIPseeker' is required. Install with:\n",
         "  BiocManager::install('ChIPseeker')")
  }

  # --- Resolve TxDb and annoDb from genome ---
  genome_db <- .resolve_genome_db(genome, txdb, anno_db)
  txdb    <- genome_db$txdb
  anno_db <- genome_db$anno_db

  # Build GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = clusters$chr,
    ranges   = IRanges::IRanges(start = clusters$start, end = clusters$end),
    strand   = if ("strand" %in% names(clusters)) clusters$strand else "*"
  )

  # Annotate
  message(">> Running ChIPseeker annotatePeak (", genome, ") ...")
  peakAnno <- ChIPseeker::annotatePeak(
    gr,
    TxDb      = txdb,
    tssRegion = tss_region,
    annoDb    = anno_db,
    verbose   = FALSE
  )

  anno_df <- as.data.frame(peakAnno)

  extra_cols <- setdiff(names(clusters), c("chr", "start", "end", "strand"))
  if (length(extra_cols) > 0) {
    for (col in extra_cols) {
      anno_df[[col]] <- clusters[[col]]
    }
  }

  if (filter_intergenic) {
    anno_df <- anno_df[anno_df$annotation == "Distal Intergenic", ]
    message(sprintf("Distal Intergenic: %d / %d clusters",
                    nrow(anno_df), nrow(clusters)))
  }

  rownames(anno_df) <- NULL
  anno_df
}


# ---------------------------------------------------------------------------
# Internal: resolve genome to TxDb + annoDb
# ---------------------------------------------------------------------------

.resolve_genome_db <- function(genome, txdb = NULL, anno_db = NULL) {

  # If user provides both, use them directly
  if (!is.null(txdb) && !is.null(anno_db)) {
    return(list(txdb = txdb, anno_db = anno_db))
  }

  # Built-in genome registry
  registry <- list(
    mm10 = list(
      txdb_pkg = "TxDb.Mmusculus.UCSC.mm10.knownGene",
      anno_db  = "org.Mm.eg.db"
    ),
    mm39 = list(
      txdb_pkg = "TxDb.Mmusculus.UCSC.mm39.knownGene",
      anno_db  = "org.Mm.eg.db"
    ),
    hg19 = list(
      txdb_pkg = "TxDb.Hsapiens.UCSC.hg19.knownGene",
      anno_db  = "org.Hs.eg.db"
    ),
    hg38 = list(
      txdb_pkg = "TxDb.Hsapiens.UCSC.hg38.knownGene",
      anno_db  = "org.Hs.eg.db"
    )
  )

  if (!genome %in% names(registry)) {
    stop("Genome '", genome, "' is not built-in. Supported: ",
         paste(names(registry), collapse = ", "), ".\n",
         "For other species, provide txdb and anno_db arguments directly.")
  }

  info <- registry[[genome]]

  # Resolve TxDb
  if (is.null(txdb)) {
    if (!requireNamespace(info$txdb_pkg, quietly = TRUE)) {
      stop("TxDb package '", info$txdb_pkg, "' not installed. Install with:\n",
           "  BiocManager::install('", info$txdb_pkg, "')")
    }
    txdb <- getExportedValue(info$txdb_pkg, info$txdb_pkg)
  }

  # Resolve annoDb
  if (is.null(anno_db)) {
    anno_db <- info$anno_db
    if (!requireNamespace(anno_db, quietly = TRUE)) {
      stop("Annotation package '", anno_db, "' not installed. Install with:\n",
           "  BiocManager::install('", anno_db, "')")
    }
  }

  list(txdb = txdb, anno_db = anno_db)
}
