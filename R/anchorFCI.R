#' @title anchorFCI: Causal Discovery using Anchors
#'  
#' @param suffStat A list of parameters used by the `indepTest` function.
#' @param indepTest A function that performs conditional independence tests and returns a p-value.
#' @param var_names Names of all variables to be represented as nodes in the graph, corresponding to dataset columns.
#' @param var_labels Labels used for `var_names` when generating plots.
#' @param prec_var_names Names of additional variables, where `prec_var_names` are prioritized before `var_names`.
#' @param prec_var_labels Labels used for `prec_var_names` when generating plots.
#' @param select_anchors A Boolean indicating whether only a subset of `prec_var_labels` deemed reliable anchors should be used in the learning process.
#' @param m.max Maximum number of variables to condition on in conditional independence tests.
#' @param alpha Significance level for the conditional independence tests.
#' @param fileid A unique identifier string added to all saved objects.
#' @param conservative A Boolean indicating if the conservative approach should be applied.
#' @param maj.rule A Boolean indicating if the majority rule should be applied.
#' @param renderAll A Boolean indicating if all relevant plots should be rendered and saved within a subfolder in `output_folder`.
#' @param width Width of the rendered figures.
#' @param height Height of the rendered figures.
#' @param output_folder Path to the folder where output files should be saved.
#' 
#' @import FCI.Utils  
#' @import pcalg  
#' @export anchorFCI
anchorFCI <- function(suffStat, indepTest, var_names, var_labels, 
                      prec_var_names = c(), prec_var_labels = c(), 
                      select_anchors = TRUE,
                      m.max=2, alpha = 0.05, fileid="fit", 
                      conservative = FALSE, maj.rule = TRUE,
                      renderAll = FALSE, width=800, height = 450,
                      output_folder="./") {

  plots_folder <- paste0(output_folder, "plots/")
  if (!file.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
  }

  all_vars <- c(var_names, prec_var_names)
  
  if (!select_anchors) {
    anchor_names = prec_var_names
  } else {
    anchor_names = selectReliableAnchors(
      suffStat, indepTest, var_names, prec_var_names,
      m.max, alpha, fileid, conservative, maj.rule)
  }
  

  mjr_rfci_fileid <- paste0(fileid, "_mjr_anchorfci_", alpha)

  p <- length(var_labels)
  fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  colnames(fixedEdges) <- rownames(fixedEdges) <- all_vars

  fixedGaps = NULL
  NAdelete = FALSE
  verbose = TRUE
  rules <- rep(TRUE, 10)

  skel <- skeleton(suffStat, indepTest, alpha, labels = all_vars,
                   method = "stable",
                   fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                   NAdelete = NAdelete, m.max = m.max, numCores = 1,
                   verbose = verbose)


  skel_amat <- as(skel@graph, "matrix")
  sepset <- fixSepsetList(skel@sepset)

  if (renderAll) {
    renderAG(skel_amat, output_folder=plots_folder,
           fileid = paste0(mjr_rfci_fileid, "_skel"), type = "pdf",
           width = width, height = height,
           labels=var_labels, add_index = FALSE)
  }

  unsh_triplets <- find.unsh.triple(skel_amat, check = FALSE)
  rfci_vstr <- rfci.vStruc(suffStat, indepTest, alpha, sepset,
                           skel_amat, unshTripl = unsh_triplets$unshTripl,
                           unshVect = unsh_triplets$unshVect,
                           conservative = (conservative || maj.rule),
                           version.unf = c(1,1), maj.rule = maj.rule,
                           verbose = verbose)

  pag_vstr <- rfci_vstr$amat
  if (renderAll) {
    renderAG(pag_vstr, output_folder=plots_folder,
           fileid = paste0(mjr_rfci_fileid, "_vstr1"), type = "pdf",
           width = width, height = height,
           labels=var_labels, add_index = FALSE)
  }

  # Ambiguous Triplets #
  ######################

  amb_triplets_df <- c()
  if (length(unsh_triplets$unshVect) > 0) {
    amb_triplets_df <- cbind.data.frame(t(unsh_triplets$unshTripl), unsh_triplets$unshVect, FALSE)
    colnames(amb_triplets_df) <- c("Vi", "Vj", "Vk", "ID", "Unf")
    amb_triplets_df$Unf[which(unsh_triplets$unshVect %in% rfci_vstr$unfTripl)] <- TRUE
    sink(paste0(output_folder, fileid, "_unf_triplets.txt"))
    print(amb_triplets_df)
    sink()
  }

  # Adding non-ancestralities -- SNP *-> Phen #
  #############################################

  edges_ids <- which(pag_vstr >=1, arr.ind = T)
  anchor_edges_ids <- edges_ids[
    (colnames(pag_vstr)[edges_ids[,1]] %in% anchor_names &
       colnames(pag_vstr)[edges_ids[,2]] %in% var_names),, drop=FALSE]
  pag_vstr[anchor_edges_ids] <- 2

  if (renderAll) {
    renderAG(pag_vstr, output_folder=plots_folder,
           fileid = paste0(mjr_rfci_fileid, "_vstr2"), type = "pdf",
           width = width, height = height,
           labels=var_labels, add_index = FALSE)
  }

  # Triplets Anchor *-> Var <-* are always unambiguous #
  ####################################################

  # triplets i,j,k, such that k is phenotype and i, j are SNPs are unambiguous colliders
  phen_ids <- which(all_vars %in% var_names)
  snps_ids <- which(all_vars %in% anchor_names)
  unamb_triplets <- expand.grid(snps_ids, snps_ids, phen_ids)
  colnames(unamb_triplets) <- c("i", "j", "k")
  unamb_triplets <- unamb_triplets[-which(unamb_triplets$i == unamb_triplets$j),]
  p <- length(all_vars)
  unamb_triplets$numb <- apply(unamb_triplets, 1, function(x) { triple2numb(p, x[1], x[2], x[3]) })

  #which(rfci_vstr$unfTripl %in% unamb_triplets$numb)

  rfci_vstr$unfTripl <- setdiff(rfci_vstr$unfTripl, unamb_triplets$numb)

  if (length(unamb_triplets$numb) > 0 && length(unsh_triplets$unshVect) > 0) {
    amb_triplets_df <- cbind.data.frame(t(unsh_triplets$unshTripl), unsh_triplets$unshVect, FALSE)
    colnames(amb_triplets_df) <- c("Vi", "Vj", "Vk", "ID", "Unf")
    amb_triplets_df$Unf[which(unsh_triplets$unshVect %in% rfci_vstr$unfTripl)] <- TRUE
    sink(paste0(output_folder, fileid, "_unf_triplets2.txt"))
    print(amb_triplets_df)
    sink()
  }


  ####################################################

  sepset <- fixSepsetList(rfci_vstr$sepset) # fit_rfci@sepset)
  #formatSepset(sepset)

  sink(paste0(output_folder, mjr_rfci_fileid, "_output.txt"))
  res <- udag2apag(pag_vstr, suffStat, indepTest, alpha, sepset,
                   rules = rules, unfVect = rfci_vstr$unfTripl, verbose = verbose)
  sink()

  pag.amat <- res$graph
  labels <- all_vars
  colnames(pag.amat) <- rownames(pag.amat) <- labels

  fit_rfci <- new("fciAlgo", amat = pag.amat,  n = integer(0), max.ord = as.integer(skel@max.ord),
                  max.ordPDSEP = 0L, n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
                  sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
  save(fit_rfci, file=paste0(output_folder, mjr_rfci_fileid, ".RData"))

  renderAG(pag.amat, output_folder=plots_folder,
           fileid = mjr_rfci_fileid, type = "pdf",
           width = width, height = height,
           labels=var_labels, add_index = FALSE)

  return(list(fit_rfci=fit_rfci, amb_triplets_df=amb_triplets_df))
}
