#' Parametric Clustering (Paraclu) — Pure R Implementation
#'
#' Finds clusters in data attached to sequences using the paraclu algorithm.
#' This is a pure R re-implementation of the C++ paraclu tool by Martin Frith.
#'
#' The algorithm recursively finds segments that maximize:
#'   score = (sum of values in segment) - d * (segment length)
#' for a range of density parameter d.
#'
#' @param ctss A data.frame with columns: chr, strand, pos, count
#'   (output from \code{bam_to_ctss}).
#' @param min_value Minimum total count for a cluster to be reported. Default: 1.
#'
#' @return A data.frame with columns: chr, strand, start, end, sites, total,
#'   min_d, max_d (same format as command-line paraclu output).
#'
#' @references Frith MC et al. (2008) A code for transcription initiation in
#'   mammalian genomes. Genome Research, 18(1), 1-12.
#'
#' @export
paraclu <- function(ctss, min_value = 1L) {

  if (!is.data.frame(ctss) || nrow(ctss) == 0) {
    return(data.frame(chr = character(), strand = character(),
                      start = integer(), end = integer(),
                      sites = integer(), total = integer(),
                      min_d = numeric(), max_d = numeric(),
                      stringsAsFactors = FALSE))
  }

  # 确保列名
  if (ncol(ctss) == 4 && !all(c("chr", "strand", "pos", "count") %in% names(ctss))) {
    names(ctss) <- c("chr", "strand", "pos", "count")
  }

  ctss <- ctss[order(ctss$chr, ctss$strand, ctss$pos), ]

  # 按 chr + strand 分组处理
  groups <- split(ctss, list(ctss$chr, ctss$strand), drop = TRUE)

  results <- lapply(groups, function(g) {
    chr_val <- g$chr[1]
    strand_val <- g$strand[1]
    positions <- g$pos
    values <- g$count

    clusters <- .paraclu_one_group(positions, values, min_value)

    if (nrow(clusters) == 0) return(NULL)

    data.frame(
      chr    = chr_val,
      strand = strand_val,
      start  = clusters$start,
      end    = clusters$end,
      sites  = clusters$sites,
      total  = clusters$total,
      min_d  = clusters$min_d,
      max_d  = clusters$max_d,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, results)
  if (is.null(result)) {
    return(data.frame(chr = character(), strand = character(),
                      start = integer(), end = integer(),
                      sites = integer(), total = integer(),
                      min_d = numeric(), max_d = numeric(),
                      stringsAsFactors = FALSE))
  }

  rownames(result) <- NULL
  message(sprintf("Paraclu: %d clusters found", nrow(result)))
  result
}


# ---------------------------------------------------------------------------
# Internal: paraclu algorithm for one chromosome + strand
# ---------------------------------------------------------------------------

.paraclu_one_group <- function(positions, values, min_value) {

  n <- length(positions)
  if (n == 0) {
    return(data.frame(start = integer(), end = integer(),
                      sites = integer(), total = integer(),
                      min_d = numeric(), max_d = numeric()))
  }

  # 收集所有聚类结果
  all_clusters <- list()
  cluster_idx <- 0L

  # 递归 max-density 切割
  .cut <- function(beg, end, min_density) {
    if (beg >= end) return()

    total <- sum(values[beg:end])
    if (total < min_value) return()

    sites <- end - beg + 1L

    if (sites == 1L) {
      cluster_idx <<- cluster_idx + 1L
      all_clusters[[cluster_idx]] <<- list(
        start = positions[beg],
        end   = positions[end],
        sites = 1L,
        total = total,
        min_d = min_density,
        max_d = Inf
      )
      return()
    }

    # 找到使两部分之间密度差最大的切割点
    # 前缀和
    prefix_sum <- cumsum(values[beg:end])
    full_sum <- prefix_sum[sites]

    # 计算每个切割点 k 的密度参数 d
    # d = prefix_sum[k] / (positions[beg+k-1] - positions[beg] + 1)
    # 但需要左右两边同时考虑

    best_density <- -Inf
    best_cut <- beg

    for (k in 1:(sites - 1L)) {
      left_sum <- prefix_sum[k]
      left_len <- positions[beg + k - 1L] - positions[beg] + 1L
      right_sum <- full_sum - left_sum
      right_len <- positions[end] - positions[beg + k] + 1L

      if (left_len > 0 && right_len > 0) {
        # 切割密度：两段中较小的密度
        d <- min(left_sum / left_len, right_sum / right_len)
        # 使用 max-min 准则
        density_val <- left_sum / left_len + right_sum / right_len
        # 实际 paraclu 使用的是: 找到使 d 最大的点,
        # 其中 segment 的 score = sum - d * length >= 0
        new_d <- min(
          if (left_len > 0) left_sum / left_len else Inf,
          if (right_len > 0) right_sum / right_len else Inf
        )
        if (new_d > best_density) {
          best_density <- new_d
          best_cut <- beg + k - 1L
        }
      }
    }

    # 记录当前聚类
    if (best_density > min_density || min_density == 0) {
      cluster_idx <<- cluster_idx + 1L
      all_clusters[[cluster_idx]] <<- list(
        start = positions[beg],
        end   = positions[end],
        sites = sites,
        total = total,
        min_d = min_density,
        max_d = best_density
      )
    }

    # 递归左右两侧
    .cut(beg, best_cut, best_density)
    .cut(best_cut + 1L, end, best_density)
  }

  .cut(1L, n, 0)

  if (cluster_idx == 0L) {
    return(data.frame(start = integer(), end = integer(),
                      sites = integer(), total = integer(),
                      min_d = numeric(), max_d = numeric()))
  }

  do.call(rbind, lapply(all_clusters, as.data.frame, stringsAsFactors = FALSE))
}
