#' Overlap Paraclu Clusters with Peak Regions
#'
#' Retains paraclu clusters that overlap with a set of peak regions
#' (e.g., ATAC-seq peaks). Equivalent to \code{bedtools intersect -u}.
#'
#' @param clusters A data.frame from \code{paraclu_cut} with at least columns:
#'   chr, strand, start, end.
#' @param peaks Either a path to a BED/narrowPeak file (character) or a
#'   GRanges object.
#' @param ignore_strand Logical; ignore strand when finding overlaps?
#'   Default: TRUE (ATAC-seq peaks are typically unstranded).
#'
#' @return A data.frame of clusters that overlap with at least one peak.
#'
#' @export
overlap_with_peaks <- function(clusters, peaks, ignore_strand = TRUE) {

  if (!is.data.frame(clusters) || nrow(clusters) == 0) {
    message("No clusters to overlap.")
    return(clusters)
  }

  # 构建 cluster GRanges
  cluster_gr <- GenomicRanges::GRanges(
    seqnames = clusters$chr,
    ranges   = IRanges::IRanges(start = clusters$start, end = clusters$end),
    strand   = if ("strand" %in% names(clusters)) clusters$strand else "*"
  )

  # 读取 peaks
  if (is.character(peaks)) {
    stopifnot(file.exists(peaks))
    peak_df <- read.table(peaks, header = FALSE, sep = "\t",
                          stringsAsFactors = FALSE, fill = TRUE)
    # BED 格式至少 3 列
    peak_gr <- GenomicRanges::GRanges(
      seqnames = peak_df[[1]],
      ranges   = IRanges::IRanges(start = as.integer(peak_df[[2]]) + 1L,
                                  end   = as.integer(peak_df[[3]])),
      strand   = if (ncol(peak_df) >= 6) peak_df[[6]] else "*"
    )
  } else if (methods::is(peaks, "GRanges")) {
    peak_gr <- peaks
  } else {
    stop("'peaks' must be a file path or GRanges object.")
  }

  # findOverlaps — 等同于 bedtools intersect -u
  hits <- GenomicRanges::findOverlaps(cluster_gr, peak_gr,
                                       ignore.strand = ignore_strand)
  overlap_idx <- unique(S4Vectors::queryHits(hits))

  filtered <- clusters[overlap_idx, ]
  rownames(filtered) <- NULL

  message(sprintf("Overlap: %d / %d clusters overlap with peaks",
                  nrow(filtered), nrow(clusters)))
  filtered
}
