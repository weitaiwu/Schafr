#' Filter BAM reads containing TSO sequence
#'
#' Identifies and retains reads that contain a Template Switching Oligo (TSO)
#' sequence at specific positions (bp 25-32 from either end), allowing for
#' mismatches via Hamming distance.
#'
#' @param bam_file Path to the input BAM file.
#' @param tso_sequence TSO sequence string. Default: "CTAACGGG".
#' @param max_mismatch Maximum allowed Hamming distance. Default: 2.
#' @param output_bam Optional output BAM path. If NULL, returns a GAlignments object.
#'
#' @return If \code{output_bam} is NULL, returns a GAlignments object of filtered reads.
#'   Otherwise writes a BAM file and returns the output path invisibly.
#'
#' @export
filter_tso <- function(bam_file,
                       tso_sequence = "CTAACGGG",
                       max_mismatch = 2L,
                       output_bam = NULL) {

  stopifnot(file.exists(bam_file))

  # TSO 及其反向互补
  tso <- Biostrings::DNAString(tso_sequence)
  tso_rc <- as.character(Biostrings::reverseComplement(tso))
  tso <- as.character(tso)
  tso_len <- nchar(tso)

  # 读取 BAM
  param <- Rsamtools::ScanBamParam(
    what = c("seq", "flag", "qname"),
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
  )
  ga <- GenomicAlignments::readGAlignments(bam_file, param = param, use.names = TRUE)

  if (length(ga) == 0) {
    message("No aligned reads found in ", bam_file)
    return(ga)
  }

  # 提取序列
  seqs <- as.character(S4Vectors::mcols(ga)$seq)
  seq_lens <- nchar(seqs)

  # Hamming 距离（向量化）
  .hamming <- function(s1_vec, s2) {
    s1_split <- strsplit(s1_vec, "")
    s2_split <- strsplit(s2, "")[[1]]
    vapply(s1_split, function(s1) sum(s1 != s2_split), integer(1))
  }

  keep <- logical(length(ga))

  # 检查正向 25-32bp 区域
  long_enough <- seq_lens >= 32L
  if (any(long_enough)) {
    forward_sub <- substr(seqs[long_enough], 25L, 24L + tso_len)
    valid_len <- nchar(forward_sub) == tso_len
    if (any(valid_len)) {
      idx <- which(long_enough)[valid_len]
      hd <- .hamming(forward_sub[valid_len], tso)
      keep[idx] <- keep[idx] | (hd <= max_mismatch)
    }
  }

  # 检查倒数 25-32bp 区域
  not_yet <- long_enough & !keep
  if (any(not_yet)) {
    rev_start <- seq_lens[not_yet] - 32L + 1L
    rev_end <- seq_lens[not_yet] - 24L
    rev_sub <- substr(seqs[not_yet], rev_start, rev_end)
    valid_len <- nchar(rev_sub) == tso_len
    if (any(valid_len)) {
      idx <- which(not_yet)[valid_len]
      hd <- .hamming(rev_sub[valid_len], tso_rc)
      keep[idx] <- keep[idx] | (hd <= max_mismatch)
    }
  }

  ga_filtered <- ga[keep]
  message(sprintf("TSO filter: %d / %d reads retained (%.1f%%)",
                  length(ga_filtered), length(ga),
                  100 * length(ga_filtered) / length(ga)))

  if (!is.null(output_bam)) {
    rtracklayer::export(ga_filtered, output_bam, format = "BAM")
    return(invisible(output_bam))
  }

  ga_filtered
}
