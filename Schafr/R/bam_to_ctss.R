#' Convert BAM/GAlignments to CTSS (Cage Tag Starting Site) counts
#'
#' Extracts the 5' end position of each read (considering strand) and
#' aggregates into counts per position, producing the input format for paraclu.
#'
#' @param input Either a path to a BAM file (character) or a GAlignments object
#'   (e.g., output from \code{filter_tso}).
#' @param output_file Optional path to write the CTSS counts file
#'   (tab-separated: chr, strand, pos, count). If NULL, returns a data.frame.
#'
#' @return A data.frame with columns: chr, strand, pos, count.
#'   If \code{output_file} is specified, also writes to file.
#'
#' @export
bam_to_ctss <- function(input, output_file = NULL) {

  if (is.character(input)) {
    stopifnot(file.exists(input))
    param <- Rsamtools::ScanBamParam(
      flag = Rsamtools::scanBamFlag(isUnmappedQuery = FALSE)
    )
    ga <- GenomicAlignments::readGAlignments(input, param = param)
  } else if (methods::is(input, "GAlignments")) {
    ga <- input
  } else {
    stop("'input' must be a BAM file path or a GAlignments object.")
  }

  if (length(ga) == 0) {
    message("No reads to process.")
    empty_df <- data.frame(chr = character(), strand = character(),
                           pos = integer(), count = integer(),
                           stringsAsFactors = FALSE)
    if (!is.null(output_file)) {
      write.table(empty_df, output_file, sep = "\t",
                  row.names = FALSE, col.names = FALSE, quote = FALSE)
    }
    return(empty_df)
  }

  # 提取 5' 端位置
  gr <- GenomicRanges::granges(ga)
  strands <- as.character(BiocGenerics::strand(gr))

  # 正链: start 即 5' 端; 负链: end 即 5' 端
  pos <- ifelse(strands == "-",
                GenomicRanges::end(gr),
                GenomicRanges::start(gr))

  dt <- data.table::data.table(
    chr    = as.character(GenomicRanges::seqnames(gr)),
    strand = strands,
    pos    = pos
  )

  # 按位置聚合计数
  ctss <- dt[, .(count = .N), by = .(chr, strand, pos)]
  data.table::setkey(ctss, chr, strand, pos)

  ctss_df <- as.data.frame(ctss)

  message(sprintf("CTSS: %d unique sites from %d reads", nrow(ctss_df), length(ga)))

  if (!is.null(output_file)) {
    write.table(ctss_df, output_file, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    message("Written to: ", output_file)
  }

  ctss_df
}
