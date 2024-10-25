rm(list=ls())

library(FCI.Utils)
library(lavaan)

library(doFuture)
library(future.apply)
library(rje)
library(dplyr)
library(RUcausal)
library(pcalg)
source("anchorRFCI.R")

n_cores <- 8 # 30
plan("multicore", workers = n_cores, gc=TRUE)
#plan("cluster", workers = n_cores)


f.args <- list()
#node_labels <- colnames(truePAGs[[1]])[1:n_nodes] # all pags have the same labels
#anchor_labels <- colnames(truePAGs[[1]])[(n_nodes+1):(n_nodes+n_anchors)] # all pags have the same labels
node_labels <- c("A", "B", "C", "D", "E") # all pags have the same labels
anchor_labels <- c("G1", "G2", "G3")
for (lab in node_labels) {
  f.args[[lab]] <- list(levels = 1) # simulating continuous phenotpes
}
for (lab in anchor_labels) {
  f.args[[lab]] <- list(levels = 3) # simulating genotypes with three levels
}

# verbose = TRUE
# n_graphs = 17
# n_anchors=3; n_nodes = 5; dir_edges_prob = 0.2; bidir_edges_prob = 0.3
generateUniqueAnchoredRandomPAGs <- function(n_graphs = 20, n_nodes = 5, n_anchors = 3,
                                     dir_edges_prob = 0.2, bidir_edges_prob = 0.3,
                                     n_circ_anchors = NULL,
                                     verbose=FALSE) {

  if (is.null(n_circ_anchors)) {
    n_circ_anchors <- sample(1:n_anchors, 1)
  }

  truePAGs <- list()
  trueMAGs <- list()
  stats <- c()

  while (length(truePAGs) < n_graphs) {
    amat.mag <- getRandomMAG(n_nodes, dir_edges_prob = dir_edges_prob, bidir_edges_prob = bidir_edges_prob)
    labels <- colnames(amat.mag)
    # renderAG(amat.mag)

    amag <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
    truePAG <- getTruePAG(amag)
    amat.pag <- truePAG@amat
    #renderAG(amat.pag)

    circles <- which(amat.pag == 1, arr.ind = T)
    circle_vars <- unique(circles[,2])
    non_circle_vars <- setdiff(1:n_nodes, circle_vars)
    #print(length(circle_vars) > n_anchors)

    if (!any(colSums(amat.pag) == 0) &&
        length(circle_vars) >= n_circ_anchors &&
        length(non_circle_vars) >= (n_anchors - n_circ_anchors)) {
      #renderAG(amat.pag)
      for (circi in 1:n_circ_anchors) {
        amat.mag <- cbind(amat.mag, rep(0, dim(amat.mag)[1]))
        amat.mag <- rbind(amat.mag, rep(0, dim(amat.mag)[2]))
        amat.mag[circle_vars[circi], n_nodes + circi] <- 3
        amat.mag[n_nodes + circi, circle_vars[circi]] <- 2
      }

      for (ncirci in 1:(n_anchors - n_circ_anchors)) {
        amat.mag <- cbind(amat.mag, rep(0, dim(amat.mag)[1]))
        amat.mag <- rbind(amat.mag, rep(0, dim(amat.mag)[2]))
        amat.mag[non_circle_vars[ncirci], n_nodes + n_circ_anchors + ncirci] <- 3
        amat.mag[n_nodes + n_circ_anchors + ncirci, non_circle_vars[ncirci]] <- 2
      }

      colnames(amat.mag) <- c(labels, paste0("G", seq(1,n_anchors)))
      rownames(amat.mag) <- colnames(amat.mag)
      #renderAG(amat.mag)

      amag <- pcalg::pcalg2dagitty(amat.mag, colnames(amat.mag), type="mag")
      truePAG <- getTruePAG(amag)
      amat.pag <- truePAG@amat
      #renderAG(amat.pag)

      circ_anchors <- which(amat.pag[(n_nodes+1):(n_nodes+n_anchors), ] == 1, arr.ind = T)
      if (n_anchors > 0 && nrow(circ_anchors) < 1) {
        next()
      }

      test_dat <- generateDatasetFromPAG(apag = amat.pag, N=1000,
                                       type="mixed", f.args=f.args)
      if (is.null(test_dat$dat)) {
        next()
      }

      mec <- MAGtoMEC(amat.mag, verbose=FALSE)
      nCK0 <- length(which(mec$CK$ord == 0))
      nCK1 <- length(which(mec$CK$ord > 0))
      nNCK0 <- length(which(mec$NCK$ord == 0))
      nNCK1 <- length(which(mec$NCK$ord > 0))



      sum_amat.pag <- amat.pag + t(amat.pag)
      edge_counts <- table(factor(
        sum_amat.pag[upper.tri(sum_amat.pag)], levels = c(0, 1, 2, 3, 4, 5, 6)))

      if (verbose) {
        cat("PAG", length(truePAGs),
            "with nCK0:", nCK0, "; nCK1+:", nCK1,
            "; nNCK0", nNCK0, "; and nNCK1", nNCK1, "\n")
      }

      cur_stats <- c(nCK0=nCK0, nCK1=nCK1, nNCK0=nNCK0, nNCK1=nNCK1, edge_counts)
      stats <- rbind(stats, cur_stats)

      truePAGs[[length(truePAGs) + 1]] <- amat.pag
      trueMAGs[[length(trueMAGs) + 1]] <- amat.mag

      #print(length(truePAGs))


      #dupl_ids <- which(duplicated(stats))
      dupl_ids <- which(duplicated(truePAGs))

      if (length(dupl_ids) > 0) {
        truePAGs <- truePAGs[-dupl_ids]
        trueMAGs <- trueMAGs[-dupl_ids]
        stats <- stats[-dupl_ids, ]
      }
    }
  }
  stats <- as.data.frame(stats)

  return(list(pags=truePAGs, mags=trueMAGs, stats=stats))
}


sims_output_folder <- "../Simulation_Results/"
if (!file.exists(sims_output_folder)) {
  dir.create(sims_output_folder, recursive = T)
}

n_nodes = 5
n_anchors = 3
n_graphs = 50
out_pags_mags_file <- paste0(sims_output_folder, "out.RData")
if (!file.exists(out_pags_mags_file)) {
  out_pags_mags <- generateUniqueAnchoredRandomPAGs(n_graphs = n_graphs,
                                          n_nodes = n_nodes, n_anchors = n_anchors,
                                          dir_edges_prob = 0.2, bidir_edges_prob = 0.3,
                                          n_circ_anchors = 2,
                                          verbose=TRUE)
  save(out_pags_mags, file=out_pags_mags_file)
} else {
  load(file=out_pags_mags_file)
}

truePAGs <- out_pags_mags$pags
trueMAGs <- out_pags_mags$mags
stats <- out_pags_mags$stats

#lapply(truePAGs, renderAG)


pag_i <- 1
n_dats = 30
m.max = 2
data_type <- "mixed"
alpha <- 0.05

# Example shown in the paper
# pag_i = 25
# dat_i = 26
# N=1000


for (N in c(500, 1000, 5000, 10000)) {
  metrics <- list()
  metrics_phens <- list()

  for (pag_i in 1:n_graphs) {

    cur_pag <- truePAGs[[pag_i]]
    cur_mag <- trueMAGs[[pag_i]]
    cur_labels <- colnames(cur_pag)
    vars_names <- cur_labels # only observed variables
    covs_names = c()
    indepTest <- mixedCITest

    dat_i = 1
    while(dat_i <= n_dats) {
      cur_sim_output_folder <- paste0(sims_output_folder, pag_i, "/", N, "/", dat_i, "/")
      if (!file.exists(cur_sim_output_folder)) {
        dir.create(cur_sim_output_folder, recursive = T)
      }

      cur_fileid <- paste0("pag", pag_i, "_N", N, "_dat", dat_i)

      cur_dat_file <- paste0(cur_sim_output_folder, "cur_dat.RData")
      if (!file.exists(cur_dat_file)) {
        adat_out <- generateDatasetFromPAG(apag = cur_pag, N=N,
                                         type=data_type, f.args=f.args)
        cur_dat <- adat_out$dat
        save(cur_dat, file = cur_dat_file)
      } else {
        cat("Loading ", cur_dat_file, "...\n")
        load(file = cur_dat_file)
      }

      if (is.null(cur_dat)) {
        dat_i = n_dats + 1 # forces going to the next pag_i
        next()
      }

      suffStat <- getMixedCISuffStat(cur_dat, vars_names, covs_names)

      citestResults_file <- paste0(cur_sim_output_folder, "citestResults.RData")
      if (!file.exists(citestResults_file)) {
        citestResults <- tryCatch({
            citestResults <- getAllCITestResults(cur_dat, indepTest, suffStat,
                                             m.max=m.max, computeProbs = FALSE,
                                             saveFiles = TRUE, fileid = cur_fileid,
                                             citestResults_folder=cur_sim_output_folder)
          },
          error = function(cond) {
            message(conditionMessage(cond))
            NULL
          })
        if (is.null(citestResults)) {
          file.remove(cur_dat_file)
          next()
        }
        save(citestResults, file = citestResults_file)
      } else {
        cat("Loading ", citestResults_file, "...\n")
        load(citestResults_file)
      }

      suffStat$citestResults <- citestResults

      rfci_file <- paste0(cur_sim_output_folder, "fit_rfci.RData")
      if (!file.exists(rfci_file)) {
        fit_rfci <- pcalg::rfci(suffStat, indepTest = indepTest,
                                skel.method = "stable", labels = cur_labels, m.max=m.max,
                                NAdelete = FALSE, #type = "normal",
                                alpha = alpha,
                                verbose = TRUE, conservative = FALSE,
                                maj.rule = TRUE)
        #renderAG(fit_rfci@amat)
        save(fit_rfci, file = rfci_file)
      } else {
        cat("Loading ", rfci_file, "...\n")
        load(rfci_file)
      }

      anchor_rfci_file <- paste0(cur_sim_output_folder, "fit_anchor_rfci.RData")
      if (!file.exists(anchor_rfci_file)) {
        fit_anchor_rfci <- anchorRFCI(suffStat, indepTest, vars_names, cur_labels,
                                    fileid=cur_fileid,
                                    anchor_names = anchor_labels, m.max=m.max,
                                    alpha = alpha,
                                    conservative = FALSE, maj.rule = TRUE,
                                    renderAll = FALSE, width=800, height = 450,
                                    output_folder=cur_sim_output_folder)
        #renderAG(fit_anchor_rfci$fit_rfci@amat)
        save(fit_anchor_rfci, file = anchor_rfci_file)
      } else {
        cat("Loading ", anchor_rfci_file, "...\n")
        load(anchor_rfci_file)
      }


      faith_anchors <- c()
      if (!is.null(fit_anchor_rfci$amb_triplets_df)) {
        faith_anchors <- fit_anchor_rfci$amb_triplets_df[which(!fit_anchor_rfci$amb_triplets_df$Unf), "Vk"]
        faith_anchors <- unique(faith_anchors[which(faith_anchors > n_nodes)])
        faith_anchors <- sort(faith_anchors)
      }
      n_faith_anchors <- length(faith_anchors)

      anchor_rfci_file_faith <- paste0(cur_sim_output_folder, "fit_anchor_rfci_faith.RData")
      if (!file.exists(anchor_rfci_file_faith)) {
        all_citestResults <- suffStat$citestResults
        faith_vars <- vars_names[c(1:n_nodes, faith_anchors)]
        suffStat$citestResults  <- extractValidCITestResults(all_citestResults, vars_names, faith_vars)

        fit_anchor_rfci_faith <- anchorRFCI(suffStat, indepTest, faith_vars, faith_vars,
                                      fileid=paste0(cur_fileid, "_faith"),
                                      anchor_names = vars_names[faith_anchors], m.max=m.max,
                                      alpha = alpha,
                                      conservative = FALSE, maj.rule = TRUE,
                                      renderAll = FALSE, width=800, height = 450,
                                      output_folder=cur_sim_output_folder)
        #renderAG(fit_anchor_rfci_faith$fit_rfci@amat)
        save(fit_anchor_rfci_faith, file = anchor_rfci_file_faith)
      } else {
        cat("Loading ", anchor_rfci_file_faith, "...\n")
        load(anchor_rfci_file_faith)
      }


      # renderAG(cur_pag, fileid = paste0("truePAG_", cur_fileid), add_index = F)
      # renderAG(cur_mag, fileid = paste0("trueMAG_", cur_fileid), add_index = F)
      # renderAG(fit_rfci@amat, fileid = paste0("rfciPAG_", cur_fileid), add_index = F)
      # renderAG(fit_anchor_rfci$fit_rfci@amat, fileid = paste0("anchorFCI_PAG_", cur_fileid), add_index = F)
      # renderAG(fit_anchor_rfci_faith$fit_rfci@amat, fileid = paste0("faith_anchorFCI_PAG_", cur_fileid), add_index = F)

      # comparing true pag with true mag using all vars
      cur_true_shd <- shd_PAG(cur_mag, cur_pag)
      cur_rfci_shd <- shd_PAG(cur_mag, fit_rfci@amat)
      #cur_rfci_shd <- shd_PAG(cur_pag, fit_rfci@amat)
      cur_anch_rfci_shd <- shd_PAG(cur_mag, fit_anchor_rfci$fit_rfci@amat)

      cur_diff_rfci_shd <- cur_rfci_shd - cur_true_shd
      cur_diff_anch_sdh <- cur_anch_rfci_shd - cur_true_shd

      metrics <- rbind.data.frame(metrics,
                       data.frame(pag = pag_i,
                                  N = N,
                                  dat = dat_i,
                                  true_shd = cur_true_shd,
                                  rfci_shd = cur_rfci_shd,
                                  anch_rfci_shd = cur_anch_rfci_shd,
                                  diff_rfci_shd = cur_diff_rfci_shd,
                                  diff_anch_sdh = cur_diff_anch_sdh))

     #print(tail(metrics))


      # comparing true pag with true mag using only phens
      cur_true_phens_shd <- shd_PAG(cur_mag[1:n_nodes, 1:n_nodes], cur_pag[1:n_nodes, 1:n_nodes])
      cur_rfci_phens_shd <- shd_PAG(cur_mag[1:n_nodes, 1:n_nodes], fit_rfci@amat[1:n_nodes, 1:n_nodes])
      #cur_rfci_phens_shd <- shd_PAG(cur_pag[1:n_nodes, 1:n_nodes], fit_rfci@amat[1:n_nodes, 1:n_nodes])

      cur_diff_rfci_phens_shd <- cur_rfci_phens_shd - cur_true_phens_shd

      # comparing true pag with true mag using only phens from fit with all anchors
      cur_anch_rfci_phens_shd <- shd_PAG(cur_mag[1:n_nodes, 1:n_nodes], fit_anchor_rfci$fit_rfci@amat[1:n_nodes, 1:n_nodes])
      cur_diff_anch_phens_shd <- cur_anch_rfci_phens_shd - cur_true_phens_shd

      # comparing true pag with true mag using only phens from fit with faithful anchors
      cur_anch_rfci_phens_faith_shd <- shd_PAG(cur_mag[1:n_nodes, 1:n_nodes], fit_anchor_rfci_faith$fit_rfci@amat[1:n_nodes, 1:n_nodes])
      cur_diff_anch_phens_faith_shd <- cur_anch_rfci_phens_faith_shd - cur_true_phens_shd

      metrics_phens <- rbind.data.frame(metrics_phens,
                                  data.frame(pag = pag_i,
                                             N = N,
                                             dat = dat_i,
                                             true_shd = cur_true_phens_shd,
                                             rfci_shd = cur_rfci_phens_shd,
                                             diff_rfci_shd = cur_diff_rfci_phens_shd,
                                             anch_rfci_shd = cur_anch_rfci_phens_shd,
                                             diff_anch_sdh = cur_diff_anch_phens_shd,
                                             n_fanchors = n_faith_anchors,
                                             anch_rfci_faith_shd = cur_anch_rfci_phens_faith_shd,
                                             diff_anch_faith_sdh = cur_diff_anch_phens_faith_shd))

      #print(tail(metrics_phens))


      dat_i = dat_i + 1
    }
  }
  save(metrics_phens, file = paste0(sims_output_folder, "metrics_phens", "_N", N, ".RData"))
  save(metrics, file = paste0(sims_output_folder, "metrics", "_N", N, ".RData"))
}



# test for H1: mu_A > mu_B
wilcox_summary <- function(grpA, grpB) {
  wtest <- wilcox.test(grpA, # group A
                       grpB, # group B
                       paired=TRUE,  #  paired Wilcoxon signed-rank test
                       # alternative = "less", # H1: mu_A < mu_B
                       alternative = "greater", # H0: mu_B < mu_A; H1: mu_A > mu_B
                       conf.int=FALSE, conf.level=0.95) # confidence interval for the
  return(wtest)
}


# Overall Metrics:

overall_metrics <- c()
overall_metrics_phen <- c()

for (N in c(500, 1000, 5000, 10000)) {
  cur_metrics_phens_file <- paste0(sims_output_folder, "metrics_phens", "_N", N, ".RData")
  load(cur_metrics_phens_file)

  cur_metrics_file <- paste0(sims_output_folder, "metrics", "_N", N, ".RData")
  load(cur_metrics_file)

  agg_metrics_df <- aggregate(metrics[, c("diff_rfci_shd", "diff_anch_sdh")],
                              by = as.list(metrics[, c("pag","N")]), FUN=mean)

  # All together, over all variables
  grpA <- subset(metrics, N == N)[, "diff_rfci_shd"]
  grpB <- subset(metrics, N == N)[, "diff_anch_sdh"]
  mins <- sapply(agg_metrics_df[, c(3,4)], min)
  maxs <- sapply(agg_metrics_df[, c(3,4)], max)
  intervs <- apply(cbind(mins, maxs), 1, function(x) {
    paste0("[", paste0(signif(x,3), collapse = ", "), "]")})

  means <- colMeans(agg_metrics_df[, c(3,4)])
  sds <- sapply(agg_metrics_df[, c(3,4)], sd)
  mean_sds <- apply(cbind(means, sds), 1, function(x) {
    paste0(signif(x,3), collapse = " \\pm ")})

  ave_diff <- signif(means[[1]] - means[[2]],3)
  pval <- signif(wilcox_summary(grpA, grpB)$p.value,3)

  overall_metrics <- rbind.data.frame(overall_metrics,
                       data.frame(N=N,
                                  minmaxRFCI=intervs[[1]],
                                  meansdRFCI=mean_sds[[1]],
                                  minmaxAnchFCI=intervs[[2]],
                                  meansdAnchFCI=mean_sds[[2]],
                                  ave_diff=ave_diff, pval=pval))

  # All together, over phenotyps
  agg_metrics_phens_df <- aggregate(metrics_phens[, c("diff_rfci_shd", #"diff_anch_sdh",
                                                      "diff_anch_faith_sdh", "n_fanchors")],
                                    by = as.list(metrics_phens[, c("pag","N")]), FUN=mean)

  grpA <- subset(metrics_phens, N == N)[, "diff_rfci_shd"]
  #grpB <- subset(metrics_phens, N == N)[, "diff_anch_sdh"]
  grpB <- subset(metrics_phens, N == N)[, "diff_anch_faith_sdh"]
  mins <- sapply(agg_metrics_phens_df[, c(3,4)], min)
  maxs <- sapply(agg_metrics_phens_df[, c(3,4)], max)
  intervs <- apply(cbind(mins, maxs), 1, function(x) {
    paste0("[", paste0(signif(x,3), collapse = ", "), "]")})

  means <- colMeans(agg_metrics_phens_df[, c(3,4)])
  sds <- sapply(agg_metrics_phens_df[, c(3,4)], sd)
  mean_sds <- apply(cbind(means, sds), 1, function(x) {
    paste0(signif(x,3), collapse = " \\pm ")})

  ave_diff <- signif(means[[1]] - means[[2]],3)
  ave_nfanchor <- signif(mean(agg_metrics_phens_df$n_fanchors),3)
  pval <- signif(wilcox_summary(grpA, grpB)$p.value,3)
  overall_metrics_phen <- rbind.data.frame(overall_metrics_phen,
                                data.frame(N=N,
                                  minmaxRFCI=intervs[[1]],
                                  meansdRFCI=mean_sds[[1]],
                                  minmaxAnchFCI=intervs[[2]],
                                  meansdAnchFCI=mean_sds[[2]],
                                  ave_anch=ave_nfanchor,
                                  ave_diff=ave_diff, pval=pval))
}

library(xtable)

overall_metrics
print(xtable(overall_metrics, floating=FALSE, latex.environments=NULL,
             display=c("d", "d", "s", "s", "s", "s", "g", "g"),
             digits=rep(3, ncol(overall_metrics) + 1)),
      math.style.exponents = TRUE,  include.rownames=FALSE)

overall_metrics_phen
print(xtable(overall_metrics_phen, #floating=FALSE, latex.environments=NULL,
             display=c("d", "d", "s", "s", "s", "s", "g", "g", "g"),
             digits=rep(3, ncol(overall_metrics_phen) + 1)),
      math.style.exponents =TRUE,  include.rownames=FALSE)



# Metrics by PAG :

agg_metrics_df <- aggregate(metrics[, c("diff_rfci_shd", "diff_anch_sdh")],
                            by = as.list(metrics[, c("pag","N")]), FUN=mean)
agg_metrics_df
mean_diff <- agg_metrics_df$diff_rfci_shd - agg_metrics_df$diff_anch_sdh

pvalues <- c()
for (pag_i in 1:50) {
  grpA <- subset(metrics, N == N & pag == pag_i)[, "diff_rfci_shd"]
  grpB <- subset(metrics, N == N & pag == pag_i)[, "diff_anch_sdh"]
  pv1 <- wilcox_summary(grpA, grpB)$p.value
  pvalues <- c(pvalues, pv1)
  cat(paste0("N: ", N, "; pag: ", pag_i, " - p:", pv1), "\n")
}

agg_metrics_df <- cbind(agg_metrics_df, diff=mean_diff, pvalue=pvalues)
agg_metrics_df

agg_metrics_df[, 3:5] <- lapply(agg_metrics_df[, 3:5], sprintf, fmt = "%0.2f")

print(xtable(agg_metrics_df, floating=FALSE, latex.environments=NULL,
             display=c("s", "d", "d", rep("s", 3), "g"),
             digits=c(0,0,0,3,3,3,3)),
      math.style.exponents = TRUE,  include.rownames=FALSE)



###################################
# Consider only phenotype network #
###################################


agg_metrics_phens_df <- aggregate(metrics_phens[, c("diff_rfci_shd", "diff_anch_sdh", "diff_anch_faith_sdh", "n_fanchors")],
                            by = as.list(metrics_phens[, c("pag","N")]), FUN=mean)

agg_metrics_phens_df
mean_diff <- agg_metrics_phens_df$diff_rfci_shd - agg_metrics_phens_df$diff_anch_sdh

pvalues <- c()
for (pag_i in 1:50) {
  grpA <- subset(metrics_phens, N == N & pag == pag_i)[, "diff_rfci_shd"]
  grpB <- subset(metrics_phens, N == N & pag == pag_i)[, "diff_anch_faith_sdh"]
  pv1 <- wilcox_summary(grpA, grpB)$p.value
  pvalues <- c(pvalues, pv1)
  cat(paste0("N: ", N, "; pag: ", pag_i, " - p:", pv1), "\n")
}

agg_metrics_phens_df <- cbind(agg_metrics_phens_df, diff=mean_diff, pvalue=pvalues)
agg_metrics_phens_df

agg_metrics_phens_df[, c(3:7)] <- lapply(agg_metrics_phens_df[, c(3:7)], sprintf, fmt = "%0.2f")

print(xtable(agg_metrics_phens_df, floating=FALSE, latex.environments=NULL,
             display=c("s", "d", "d", rep("s", 5), "g"),
             digits=c(0,0,0,3,3,3,3,3,3)),
      math.style.exponents = TRUE,  include.rownames=FALSE)









