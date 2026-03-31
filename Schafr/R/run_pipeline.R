#' Run the Complete Schafr Pipeline
#'
#' One-function workflow: TSO filtering → CTSS counting → Paraclu clustering →
#' Paraclu-cut filtering → ATAC peak overlap → (optional) ChIPseeker annotation.
#'
#' @param bam_files Character vector of BAM file paths (one per cell type).
#' @param atac_peaks Path to ATAC-seq peak file (BED/narrowPeak) or GRanges.
#' @param tso_sequence TSO sequence. Default: "CTAACGGG".
#' @param max_mismatch Max Hamming distance for TSO filter. Default: 2.
#' @param min_paraclu_value Min count for paraclu. Default: 1.
#' @param max_cluster_length Max cluster length for paraclu-cut. Default: 200.
#' @param min_density_ratio Min density ratio for paraclu-cut. Default: 2.
#' @param annotate Logical; run ChIPseeker annotation? Default: FALSE.
#' @param filter_intergenic Logical; if annotating, filter to Distal Intergenic
#'   only? Default: FALSE.
#' @param genome Genome identifier: "mm10", "mm39", "hg19", "hg38".
#'   Default: "mm10". Ignored if txdb and anno_db are both provided.
#' @param txdb TxDb object. If NULL, auto-resolved from \code{genome}.
#' @param anno_db Annotation DB. If NULL, auto-resolved from \code{genome}.
#' @param output_dir Optional directory to write intermediate and final files.
#'
#' @return A named list, one element per BAM file, each containing:
#'   \describe{
#'     \item{ctss}{CTSS count data.frame}
#'     \item{clusters_raw}{Raw paraclu clusters}
#'     \item{clusters_cut}{Filtered paraclu clusters}
#'     \item{clusters_overlap}{Clusters overlapping ATAC peaks}
#'     \item{clusters_annotated}{(if annotate=TRUE) Annotated clusters}
#'   }
#'
#' @examples
#' \dontrun{
#' results <- run_pipeline(
#'   bam_files  = c("cluster_0.bam", "cluster_1.bam"),
#'   atac_peaks = "ATAC_peaks.narrowPeak",
#'   annotate   = TRUE,
#'   filter_intergenic = TRUE,
#'   output_dir = "schafr_output"
#' )
#' }
#'
#' @export
run_pipeline <- function(bam_files,
                         atac_peaks,
                         tso_sequence = "CTAACGGG",
                         max_mismatch = 2L,
                         min_paraclu_value = 1L,
                         max_cluster_length = 200L,
                         min_density_ratio = 2,
                         annotate = FALSE,
                         filter_intergenic = FALSE,
                         genome = "mm10",
                         txdb = NULL,
                         anno_db = NULL,
                         output_dir = NULL) {

  stopifnot(length(bam_files) > 0)

  if (!is.null(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  results <- list()

  for (i in seq_along(bam_files)) {
    bam <- bam_files[i]
    sample_name <- tools::file_path_sans_ext(basename(bam))

    message("\n", strrep("=", 60))
    message(sprintf(" [%d/%d] Processing: %s", i, length(bam_files), sample_name))
    message(strrep("=", 60))

    res <- list()

    # Step 1: TSO filter
    message("\n>> Step 1/5: TSO filtering ...")
    ga_filtered <- filter_tso(
      bam_file      = bam,
      tso_sequence  = tso_sequence,
      max_mismatch  = max_mismatch
    )

    # Step 2: BAM → CTSS
    message("\n>> Step 2/5: Generating CTSS ...")
    ctss_file <- if (!is.null(output_dir)) {
      file.path(output_dir, paste0(sample_name, "_ctss.txt"))
    } else NULL

    res$ctss <- bam_to_ctss(ga_filtered, output_file = ctss_file)

    # Step 3: Paraclu
    message("\n>> Step 3/5: Running Paraclu ...")
    res$clusters_raw <- paraclu(res$ctss, min_value = min_paraclu_value)

    # Step 4: Paraclu-cut
    message("\n>> Step 4/5: Paraclu-cut filtering ...")
    res$clusters_cut <- paraclu_cut(
      res$clusters_raw,
      max_length        = max_cluster_length,
      min_density_ratio = min_density_ratio
    )

    # Step 5: ATAC overlap
    message("\n>> Step 5/5: Overlapping with ATAC peaks ...")
    res$clusters_overlap <- overlap_with_peaks(
      res$clusters_cut,
      peaks         = atac_peaks,
      ignore_strand = TRUE
    )

    # Optional: ChIPseeker annotation
    if (annotate) {
      message("\n>> Bonus: ChIPseeker annotation ...")
      res$clusters_annotated <- annotate_intergenic(
        res$clusters_overlap,
        genome             = genome,
        txdb               = txdb,
        anno_db            = anno_db,
        filter_intergenic  = filter_intergenic
      )
    }

    # 写出最终结果
    if (!is.null(output_dir)) {
      final <- if (annotate && filter_intergenic) {
        res$clusters_annotated
      } else {
        res$clusters_overlap
      }

      out_file <- file.path(output_dir, paste0(sample_name, "_final.bed"))
      if (nrow(final) > 0) {
        bed_cols <- intersect(c("chr", "start", "end", "sites", "total", "strand"),
                              names(final))
        write.table(final[, bed_cols, drop = FALSE], out_file, sep = "\t",
                    row.names = FALSE, col.names = FALSE, quote = FALSE)
        message(">> Final result: ", out_file)
      }
    }

    results[[sample_name]] <- res
  }

  message("\n", strrep("=", 60))
  message(" Pipeline complete! Processed ", length(bam_files), " sample(s).")
  message(strrep("=", 60))

  invisible(results)
}
