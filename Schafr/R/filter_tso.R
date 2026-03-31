#' Filter BAM reads containing TSO sequence
#'
#' Identifies and retains reads that contain a Template Switching Oligo (TSO)
#' sequence at customizable positions, allowing for mismatches via Hamming distance.
#' Checks both the forward strand position and the reverse complement at the
#' symmetric position from the read end.
#'
#' @param bam_file Path to the input BAM file.
#' @param tso_sequence TSO sequence string. Default: "CTAACGGG".
#' @param max_mismatch Maximum allowed Hamming distance. Default: 2.
#' @param tso_start Start position (1-based) to check for TSO in the read.
#'   Default: 25.
#' @param tso_end End position (1-based) to check for TSO in the read.
#'   Default: 32. The window length (tso_end - tso_start + 1) should equal
#'   \code{nchar(tso_sequence)}.
#' @param check_reverse Logical; also check the symmetric position from the
#'   read end for the reverse complement? Default: TRUE.
#' @param output_bam Optional output BAM path. If NULL, returns a GAlignments object.
#'
#' @details
#' The forward check extracts \code{seq[tso_start:tso_end]} and compares to
#' \code{tso_sequence}. The reverse check extracts
#' \code{seq[(len-tso_end+1):(len-tso_start+1)]} and compares to the reverse
#' complement of \code{tso_sequence}.
#'
#' @return If \code{output_bam} is NULL, returns a GAlignments object of filtered reads.
#'   Otherwise writes a BAM file and returns the output path invisibly.
#'
#' @examples
#' \dontrun{
#' # Default: check position 25-32
#' ga <- filter_tso("cluster_0.bam")
#'
#' # Custom TSO and position
#' ga <- filter_tso("cluster_0.bam",
#'   tso_sequence = "AAGCAGTGGTATCAACGCAGAGT",
#'   tso_start = 1, tso_end = 23, max_mismatch = 3)
#'
#' # Only check forward, not reverse
#' ga <- filter_tso("cluster_0.bam", check_reverse = FALSE)
#' }
#'
#' @export
filter_tso <- function(bam_file,
                       tso_sequence = "CTAACGGG",
                       max_mismatch = 2L,
                       tso_start = 25L,
                       tso_end = 32L,
                       check_reverse = FALSE,
                       output_bam = NULL) {

  stopifnot(file.exists(bam_file))

  # Validate position window vs TSO length
  tso <- Biostrings::DNAString(tso_sequence)
  tso_rc <- as.character(Biostrings::reverseComplement(tso))
  tso <- as.character(tso)
  tso_len <- nchar(tso)

  window_len <- tso_end - tso_start + 1L
  if (window_len != tso_len) {
    stop(sprintf(
      "Window size (tso_start=%d to tso_end=%d = %d bp) does not match TSO length (%d bp).",
      tso_start, tso_end, window_len, tso_len))
  }

  # Read BAM
  param <- Rsamtools::ScanBamParam(
    what = c("seq", "flag", "qname"),
    flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
  )
  ga <- GenomicAlignments::readGAlignments(bam_file, param = param, use.names = TRUE)

  if (length(ga) == 0) {
    message("No aligned reads found in ", bam_file)
    return(ga)
  }

  # Extract sequences
  seqs <- as.character(S4Vectors::mcols(ga)$seq)
  seq_lens <- nchar(seqs)

  # Vectorized Hamming distance
  .hamming <- function(s1_vec, s2) {
    s1_split <- strsplit(s1_vec, "")
    s2_split <- strsplit(s2, "")[[1]]
    vapply(s1_split, function(s1) sum(s1 != s2_split), integer(1))
  }

  keep <- logical(length(ga))

  # Forward check: position tso_start to tso_end
  min_len_fwd <- tso_end
  long_enough <- seq_lens >= min_len_fwd
  if (any(long_enough)) {
    forward_sub <- substr(seqs[long_enough], tso_start, tso_end)
    valid_len <- nchar(forward_sub) == tso_len
    if (any(valid_len)) {
      idx <- which(long_enough)[valid_len]
      hd <- .hamming(forward_sub[valid_len], tso)
      keep[idx] <- keep[idx] | (hd <= max_mismatch)
    }
  }

  # Reverse check: symmetric position from end
  if (check_reverse) {
    not_yet <- long_enough & !keep
    if (any(not_yet)) {
      rev_start <- seq_lens[not_yet] - tso_end + 1L
      rev_end   <- seq_lens[not_yet] - tso_start + 1L
      rev_sub <- substr(seqs[not_yet], rev_start, rev_end)
      valid_len <- nchar(rev_sub) == tso_len
      if (any(valid_len)) {
        idx <- which(not_yet)[valid_len]
        hd <- .hamming(rev_sub[valid_len], tso_rc)
        keep[idx] <- keep[idx] | (hd <= max_mismatch)
      }
    }
  }

  ga_filtered <- ga[keep]
  message(sprintf("TSO filter (pos %d-%d): %d / %d reads retained (%.1f%%)",
                  tso_start, tso_end,
                  length(ga_filtered), length(ga),
                  100 * length(ga_filtered) / length(ga)))

  if (!is.null(output_bam)) {
    rtracklayer::export(ga_filtered, output_bam, format = "BAM")
    return(invisible(output_bam))
  }

  ga_filtered
}
