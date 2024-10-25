anchorFCI <- function(suffStat, indepTest, cur_vars, cur_labels, fileid="fit",
                       anchor_names = c(), m.max=2,
                       alpha = 0.05,
                       conservative = FALSE, maj.rule = TRUE,
                       renderAll = FALSE, width=800, height = 450,
                       output_folder="./") {

  plots_folder <- paste0(output_folder, "plots/")
  if (!file.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
  }

  noanchor_names = setdiff(cur_vars, anchor_names)

  mjr_rfci_fileid <- paste0(fileid, "_mjr_anchorfci_", alpha)

  p <- length(cur_labels)
  fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
  colnames(fixedEdges) <- rownames(fixedEdges) <- cur_vars

  fixedGaps = NULL
  NAdelete = FALSE
  verbose = TRUE
  rules <- rep(TRUE, 10)

  skel <- skeleton(suffStat, indepTest, alpha, labels = cur_vars,
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
           labels=cur_labels, add_index = FALSE)
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
           labels=cur_labels, add_index = FALSE)
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
       colnames(pag_vstr)[edges_ids[,2]] %in% noanchor_names),, drop=FALSE]
  pag_vstr[anchor_edges_ids] <- 2

  if (renderAll) {
    renderAG(pag_vstr, output_folder=plots_folder,
           fileid = paste0(mjr_rfci_fileid, "_vstr2"), type = "pdf",
           width = width, height = height,
           labels=cur_labels, add_index = FALSE)
  }

  # Triplets Anchor *-> Var <-* are always unambiguous #
  ####################################################

  # triplets i,j,k, such that k is phenotype and i, j are SNPs are unambiguous colliders
  phen_ids <- which(cur_vars %in% noanchor_names)
  snps_ids <- which(cur_vars %in% anchor_names)
  unamb_triplets <- expand.grid(snps_ids, snps_ids, phen_ids)
  colnames(unamb_triplets) <- c("i", "j", "k")
  unamb_triplets <- unamb_triplets[-which(unamb_triplets$i == unamb_triplets$j),]
  p <- length(cur_vars)
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
  labels <- cur_vars
  colnames(pag.amat) <- rownames(pag.amat) <- labels

  fit_rfci <- new("fciAlgo", amat = pag.amat,  n = integer(0), max.ord = as.integer(skel@max.ord),
                  max.ordPDSEP = 0L, n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
                  sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
  save(fit_rfci, file=paste0(output_folder, mjr_rfci_fileid, ".RData"))

  renderAG(pag.amat, output_folder=plots_folder,
           fileid = mjr_rfci_fileid, type = "pdf",
           width = width, height = height,
           labels=cur_labels, add_index = FALSE)

  return(list(fit_rfci=fit_rfci, amb_triplets_df=amb_triplets_df))
}
