#' Annotate Clusters and Filter Distal Intergenic Regions
#'
#' Uses ChIPseeker to annotate paraclu clusters and optionally filters to
#' retain only those in Distal Intergenic regions.
#'
#' @param clusters A data.frame with at least columns: chr, start, end, strand.
#' @param txdb A TxDb object for annotation. Default uses
#'   TxDb.Mmusculus.UCSC.mm10.knownGene if available.
#' @param anno_db Annotation database name. Default: "org.Mm.eg.db".
#' @param tss_region Numeric vector of length 2 defining TSS region.
#'   Default: c(-3000, 3000).
#' @param filter_intergenic Logical; if TRUE, only return clusters annotated as
#'   "Distal Intergenic". Default: FALSE.
#'
#' @return A data.frame with annotation columns appended. If
#'   \code{filter_intergenic = TRUE}, only Distal Intergenic clusters are returned.
#'
#' @export
annotate_intergenic <- function(clusters,
                                txdb = NULL,
                                anno_db = "org.Mm.eg.db",
                                tss_region = c(-3000, 3000),
                                filter_intergenic = FALSE) {

  if (!requireNamespace("ChIPseeker", quietly = TRUE)) {
    stop("Package 'ChIPseeker' is required. Install with:\n",
         "  BiocManager::install('ChIPseeker')")
  }

  if (is.null(txdb)) {
    if (!requireNamespace("TxDb.Mmusculus.UCSC.mm10.knownGene", quietly = TRUE)) {
      stop("Default TxDb not available. Install with:\n",
           "  BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')")
    }
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene
  }

  # 构建 GRanges
  gr <- GenomicRanges::GRanges(
    seqnames = clusters$chr,
    ranges   = IRanges::IRanges(start = clusters$start, end = clusters$end),
    strand   = if ("strand" %in% names(clusters)) clusters$strand else "*"
  )

  # 注释
  message(">> Running ChIPseeker annotatePeak ...")
  peakAnno <- ChIPseeker::annotatePeak(
    gr,
    TxDb      = txdb,
    tssRegion = tss_region,
    annoDb    = anno_db,
    verbose   = FALSE
  )

  anno_df <- as.data.frame(peakAnno)

  # 合并回原始 clusters 的额外列（如 sites, total, min_d, max_d）
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
