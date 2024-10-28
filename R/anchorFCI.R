#' @title anchorFCI: Causal Discovery using Anchors
#'  
#' @param suffStat A list of parameters used by the `indepTest` function.
#' @param indepTest A function that performs conditional independence tests and returns a p-value.
#' @param var_names Names of all variables to be represented as nodes in the graph, corresponding to dataset columns.
#' @param var_labels Labels used for `var_names` when generating plots.
#' @param prec_var_names Names of additional variables, where `prec_var_names` precede `var_names` in the partial ordering.
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
                      output_folder="./", verbose = TRUE) {

  if (is.null(var_labels)) {
    var_labels = var_names
  }
  
  if (is.null(prec_var_labels)) {
    prec_var_labels = prec_var_names
  }
  
  anchor_fci_fileid <- paste0(fileid, "_anchorfci_", alpha)
  
  
  # suffStat is a list with original suffStat and indepTest, and anchor_ids
  anchorCITest <- function(x, y, S, suffStat) {
    if (x %in% suffStat$anchor_ids && y %in% suffStat$anchor_ids) {
      if (!all(S %in% suffStat$anchor_ids)) {
        cat("Forbidden test between ", x, "and", y, 
            "given", paste0(S, collapse=","), "\n")
        return(0) # indicating dependence
      }
    }
    return(suffStat$indepTest(x,y,S,suffStat$suffStat))
  }

  getR0Skeleton <- function(all_vars, all_var_labels, anchor_ids) {
  
    p <- length(all_vars)
    fixedEdges <- matrix(rep(FALSE, p * p), nrow = p, ncol = p)
    colnames(fixedEdges) <- rownames(fixedEdges) <- all_vars
  
    fixedGaps = NULL
    NAdelete = FALSE
    
    if (!(length(orig_all_vars) == length(all_vars) && 
          all(orig_all_vars == all_vars))) {
      if (!is.null(all_citestResults)) {
        suffStat$citestResults <- extractValidCITestResults(
          all_citestResults, orig_all_vars, all_vars)
      }
    }
    
    anchor_suffStat <- list(suffStat=suffStat, indepTest=indepTest, 
                            anchor_ids = anchor_ids)
  
    skel <- skeleton(anchor_suffStat, anchorCITest, alpha, 
                     labels = all_vars, method = "stable",
                     fixedGaps = fixedGaps, fixedEdges = fixedEdges,
                     NAdelete = NAdelete, m.max = m.max, numCores = 1,
                     verbose = verbose)
  
  
    skel_amat <- as(skel@graph, "matrix")
    sepset <- fixSepsetList(skel@sepset)
  
    if (renderAll) {
      renderAG(skel_amat, output_folder=plots_folder,
             fileid = paste0(anchor_fci_fileid, "_skel"), type = "pdf",
             width = width, height = height,
             labels=all_var_labels, add_index = FALSE)
    }
  
    unsh_triplets <- find.unsh.triple(skel_amat, check = FALSE)
    rfci_vstr <- rfci.vStruc(anchor_suffStat, anchorCITest, alpha, sepset,
                             skel_amat, unshTripl = unsh_triplets$unshTripl,
                             unshVect = unsh_triplets$unshVect,
                             conservative = (conservative || maj.rule),
                             version.unf = c(1,1), maj.rule = maj.rule,
                             verbose = verbose)

    pag_vstr <- rfci_vstr$amat
    if (renderAll) {
      renderAG(pag_vstr, output_folder=plots_folder,
             fileid = paste0(anchor_fci_fileid, "_vstr1"), type = "pdf",
             width = width, height = height,
             labels=all_var_labels, add_index = FALSE)
    }

    ##################################
    # Identifying Ambiguous Triplets # 
    ##################################
  
    amb_triplets_df <- c()
    if (length(unsh_triplets$unshVect) > 0) {
      amb_triplets_df <- cbind.data.frame(t(unsh_triplets$unshTripl), 
                                          unsh_triplets$unshVect, FALSE)
      colnames(amb_triplets_df) <- c("Vi", "Vj", "Vk", "ID", "Unf")
      amb_triplets_df$Unf[
        which(unsh_triplets$unshVect %in% rfci_vstr$unfTripl)] <- TRUE
      sink(paste0(output_folder, fileid, "_unf_triplets.txt"))
      print(amb_triplets_df)
      sink()
    }
    
    return(list(rfci_vstr=rfci_vstr, skel=skel, amb_triplets_df=amb_triplets_df))
  }
  
  
  
  plots_folder <- paste0(output_folder, "plots/")
  if (!file.exists(plots_folder)) {
    dir.create(plots_folder, recursive = TRUE)
  }
  
  
  orig_all_vars <- c(var_names, prec_var_names)
  orig_all_var_labels <-  c(var_labels, prec_var_labels)
  all_citestResults <- NULL
  if (!is.null(suffStat$citestResults)) {
    all_citestResults <- suffStat$citestResults
  }
  
  
  all_vars <- orig_all_vars
  all_var_labels <- orig_all_var_labels

  prec_vars_ids <- which(orig_all_vars %in% prec_var_names)
  
  outR0Skel <- getR0Skeleton(orig_all_vars, orig_all_var_labels, anchor_ids = prec_vars_ids)
  amb_triplets_df <- outR0Skel$amb_triplets_df
  rfci_vstr <- outR0Skel$rfci_vstr
  skel <- outR0Skel$skel
  
  ##############################
  # Selecting Reliable Anchors #
  ##############################
  
  if (!select_anchors) {
    anchor_ids <- prec_vars_ids
  } else {
    anchor_ids <- c()
    if (!is.null(amb_triplets_df)) {
      anch_triplets_df <- subset(amb_triplets_df, 
                                 Vi %in% prec_vars_ids | Vj %in% prec_vars_ids | 
                                   Vk %in% prec_vars_ids)
      anch_triplets_df <- subset(anch_triplets_df, Unf == FALSE, drop = FALSE)

      anchor_ids <- sort(unique(unlist(c(anch_triplets_df[, c("Vi", "Vj", "Vk")]))))
      anchor_ids <- anchor_ids[which(anchor_ids %in% prec_vars_ids)]
    }
  }
  
  n_anchors <- length(anchor_ids)
   
  #anchor_suffStat <- list(suffStat=suffStat, indepTest=indepTest, 
  #                        anchor_ids = anchor_ids)
  
  unf_prec_var_ids <- setdiff(prec_vars_ids, anchor_ids)
  if (length(unf_prec_var_ids) > 0) {
    all_vars <- orig_all_vars[-c(unf_prec_var_ids)]
    all_var_labels <- orig_all_var_labels[-c(unf_prec_var_ids)]
    
    anchor_ids = which(all_vars %in% prec_var_names)
    outR0Skel <- getR0Skeleton(all_vars, all_var_labels, anchor_ids = anchor_ids)
    amb_triplets_df <- outR0Skel$amb_triplets_df
    rfci_vstr <- outR0Skel$rfci_vstr
    skel <- outR0Skel$skel
  }  
  
  anchor_names = all_vars[which(all_vars %in% prec_var_names)]
  cat("Selected", n_anchors, "reliable anchors:", 
      paste0(anchor_names, collapse = ", "), "\n")
  
  
  ###############################################
  # Adding non-ancestralities -- Anchor *-> Var #
  ###############################################

  pag_vstr <- rfci_vstr$amat
  
  edges_ids <- which(pag_vstr >=1, arr.ind = T)
  anchor_edges_ids <- edges_ids[
    (colnames(pag_vstr)[edges_ids[,1]] %in% anchor_names &
       colnames(pag_vstr)[edges_ids[,2]] %in% var_names),, drop=FALSE]
  pag_vstr[anchor_edges_ids] <- 2

  if (renderAll) {
    renderAG(pag_vstr, output_folder=plots_folder,
           fileid = paste0(anchor_fci_fileid, "_vstr2"), type = "pdf",
           width = width, height = height,
           labels=all_var_labels, add_index = FALSE)
  }

  # Triplets Anchor *-> Var <-* Anchor are always unambiguous #
  #############################################################

  # triplets i,j,k, such that k is phenotype and i, j are SNPs are unambiguous colliders
  var_ids <- which(all_vars %in% var_names)
  anchor_ids <- which(all_vars %in% anchor_names)
  unamb_triplets <- expand.grid(anchor_ids, anchor_ids, var_ids)
  colnames(unamb_triplets) <- c("i", "j", "k")
  unamb_triplets <- unamb_triplets[-which(unamb_triplets$i == unamb_triplets$j),]
  p <- length(all_vars)
  unamb_triplets$numb <- apply(unamb_triplets, 1, function(x) { triple2numb(p, x[1], x[2], x[3]) })


  rfci_vstr$unfTripl <- setdiff(rfci_vstr$unfTripl, unamb_triplets$numb)

  if (length(unamb_triplets$numb) > 0 && !is.null(amb_triplets_df)) {
    amb_triplets_df$Unf[which(amb_triplets_df$ID %in% unamb_triplets$numb)] <- FALSE
    sink(paste0(output_folder, fileid, "_unf_triplets2.txt"))
    print(amb_triplets_df)
    sink()
  }

  ####################################################

  sepset <- fixSepsetList(rfci_vstr$sepset) 
  #formatSepset(sepset)

  rules <- rep(TRUE, 10)
  sink(paste0(output_folder, anchor_fci_fileid, "_output.txt"))
  res <- udag2apag(pag_vstr, suffStat, indepTest, alpha, sepset,
                   rules = rules, unfVect = rfci_vstr$unfTripl, verbose = verbose)
  sink()

  pag.amat <- res$graph
  colnames(pag.amat) <- rownames(pag.amat) <- all_vars

  fit_anchor_fci <- new("fciAlgo", amat = pag.amat,  n = integer(0), max.ord = as.integer(skel@max.ord),
                  max.ordPDSEP = 0L, n.edgetests = skel@n.edgetests, n.edgetestsPDSEP = 0,
                  sepset = res$sepset, pMax = skel@pMax, allPdsep = vector("list", p))
  save(fit_anchor_fci, file=paste0(output_folder, anchor_fci_fileid, ".RData"))

  renderAG(pag.amat, output_folder=plots_folder,
           fileid = anchor_fci_fileid, type = "pdf",
           width = width, height = height,
           labels=all_var_labels, add_index = FALSE)

  return(list(fit_anchor_fci=fit_anchor_fci, amb_triplets_df=amb_triplets_df))
}
