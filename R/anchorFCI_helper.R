removeForbiddenCI <- function(citestResults, anchors_ids, m.max=Inf) {
  forb_tests <- c()
  if (!is.null(citestResults) && !is.null(anchors_ids) && length(anchors_ids) > 0) {
    anchor_pairs <- which(citestResults$X %in% anchors_ids & citestResults$Y %in% anchors_ids)

    Slist <- sapply(citestResults$S, getSepVector)
    forb_anchor_seps <- which(sapply(Slist, function(x) { !all(x %in% anchors_ids)} ))
    forb_tests <- intersect(anchor_pairs, forb_anchor_seps)
  }

  if (!is.null(citestResults)) {
    forb_tests <- unique(c(which(citestResults$ord > m.max), forb_tests))
    if (length(forb_tests) > 0) {
      citestResults[forb_tests, c("pvalue")] <-
        matrix(c(rep(0, length(forb_tests))), nrow=length(forb_tests), ncol=1)
    }
  }
  return(citestResults)
}
