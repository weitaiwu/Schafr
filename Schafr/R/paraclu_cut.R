#' Filter Paraclu Clusters (paraclu-cut)
#'
#' Filters paraclu output by removing clusters that are too long, have
#' insufficient density increase, or are entirely contained within another
#' cluster.
#'
#' @param clusters A data.frame from \code{paraclu} with columns:
#'   chr, strand, start, end, sites, total, min_d, max_d.
#' @param max_length Maximum allowed cluster length (end - start + 1).
#'   Default: 200.
#' @param min_density_ratio Minimum ratio of max_d / min_d required.
#'   Default: 2.
#' @param remove_contained Logical; remove clusters fully contained within
#'   another cluster? Default: TRUE.
#'
#' @return A filtered data.frame with the same columns.
#'
#' @export
paraclu_cut <- function(clusters,
                        max_length = 200L,
                        min_density_ratio = 2,
                        remove_contained = TRUE) {

  if (!is.data.frame(clusters) || nrow(clusters) == 0) {
    return(clusters)
  }

  n_start <- nrow(clusters)

  # 计算聚类长度
  clusters$length <- clusters$end - clusters$start + 1L

  # 过滤 1: 长度
  clusters <- clusters[clusters$length <= max_length, ]

  # 过滤 2: 密度增幅
  # min_d == 0 时跳过比率检查（避免除零），但仍保留
  valid_ratio <- ifelse(
    clusters$min_d > 0,
    clusters$max_d / clusters$min_d >= min_density_ratio,
    TRUE
  )
  # 对 max_d == Inf 的单位点聚类，总是通过
  valid_ratio[is.infinite(clusters$max_d)] <- TRUE
  clusters <- clusters[valid_ratio, ]

  # 过滤 3: 去除被包含的聚类
  if (remove_contained && nrow(clusters) > 1) {
    clusters <- .remove_contained(clusters)
  }

  # 清理辅助列
  clusters$length <- NULL
  rownames(clusters) <- NULL

  message(sprintf("paraclu_cut: %d → %d clusters retained", n_start, nrow(clusters)))
  clusters
}


# ---------------------------------------------------------------------------
# Internal: remove clusters entirely contained within another
# ---------------------------------------------------------------------------

.remove_contained <- function(clusters) {
  # 按 chr, strand, start 升序, end 降序 排序
  clusters <- clusters[order(clusters$chr, clusters$strand,
                             clusters$start, -clusters$end), ]

  keep <- rep(TRUE, nrow(clusters))

  groups <- split(seq_len(nrow(clusters)),
                  list(clusters$chr, clusters$strand), drop = TRUE)

  for (idx in groups) {
    if (length(idx) <= 1) next

    for (i in seq_along(idx)[-1]) {
      ci <- idx[i]
      for (j in seq_len(i - 1)) {
        cj <- idx[j]
        if (!keep[cj]) next
        # 如果 ci 完全被 cj 包含
        if (clusters$start[ci] >= clusters$start[cj] &&
            clusters$end[ci] <= clusters$end[cj]) {
          keep[ci] <- FALSE
          break
        }
      }
    }
  }

  clusters[keep, ]
}
