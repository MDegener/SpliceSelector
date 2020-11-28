#' Generating a correlation plot
#'
#' @param x Matrix or dataframe containing only numeric data
#' @param y An optional second matrix or dataframe with the same number of rows as x
#' @param method Options are "pearson", "spearman" or "kendall"
#' @param adjust Method for p-value correction; Options are "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr" and "none"
#'
#' @export

makeCorrplot <- function(x, y=NULL, method="pearson", adjust="fdr"){

  assertthat::assert_that(is.numeric(x))
  assertthat::assert_that(is.numeric(y) || is.null(y))
  assertthat::assert_that(method %in% c("pearson", "spearman", "kendall"),
                          msg="'method' argument must be either \"pearson\", \"spearman\", \"kendall\"")
  assertthat::assert_that(adjust %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                          msg="'adjust' argument must be either \"holm\", \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\", \"none\"")

  corr <- psych::corr.test(x, y, method = method, adjust = adjust)

  # plot above diagonal correlations
  corrplot::corrplot(corr$r, p.mat = corr$p,
                     order = "original",
                     method = "color",
                     col = colorRampPalette(RColorBrewer::brewer.pal(n = 8, name = "PuOr"))(100),
                     tl.pos = "td",
                     tl.col="black",
                     tl.srt=60, # label rotation
                     type = ifelse(is.null(y), "upper", "full"),
                     insig = "label_sig", pch = "*", pch.cex = 0.7, sig.level = 0.05,
                     diag = TRUE)
}

