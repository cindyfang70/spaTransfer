#' Run NMF on a dataset
#'
#' @param data SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param k integer indicating the number of factors for NMF
#' @param seed a random seed
#' @param ... additional arguments passed to singlet::run_nmf or singlet::RunNMF
#'
#' @return NMF model object
#'
#' @import SingleCellExperiment
#' @import RcppML
#' @export
run_nmf<- function(data, assay, k=NULL, seed=1237,...){

  #add in checks for k
  message("Running NMF")

  if(is.null(k)){
    warning("Number of factors for NMF not specified. Using cross-validation to idenitfy optimal number of factors.", immediate. = TRUE)
      #k <- find_num_factors(A)
    set.seed(seed)
    model <- run_rank_determination_nmf(data, assay,...)
  }else{
    A <- as.matrix(assay(data, assay))
    set.seed(seed)
    model <- singlet::run_nmf(A, rank=k,...)
  }

  return(model)
}

find_num_factors <- function(A, ranks = c(50,100,200)){
  if(any(ranks >= ncol(A))){
    stop("ranks must be less than the number of columns in A")
    }
  cv <- singlet::cross_validate_nmf(A, ranks = ranks,
                                    n_replicates = 3,
                                    verbose=3)
  num_factors <- singlet::GetBestRank(cv)

  return(num_factors)
}

#' run_rank_determination_nmf
#'
#' @param data A SingleCellExperiment or SpatialExperiment object
#' @param assay string indicating the assay to run NMF on
#' @param ... additional arguments passed to singlet::RunNMF
#'
#' @return a NMF model object
#' @import singlet SingleCellExperiment
run_rank_determination_nmf <- function(data, assay,...){
  data_nmf <- RunNMF(data, assay=assay,...)
  nmf_mod <- S4Vectors::metadata(data_nmf)$nmf_model
  return(nmf_mod)
}

#' Project a target dataset onto the source NMF factors.
#'
#' The projected weights are returned on the same scale as the source factors
#' `t(nmf_model$h)`, so that a model fitted on the source can be applied to them.
#'
#' `singlet` factorises the source as `A ~ w %*% diag(d) %*% h` with the rows of
#' `h` summing to 1, so `d` carries the per-factor scale. `RcppML::project()`
#' solves `A_target ~ w %*% h_new`, meaning `h_new` absorbs that scale and must
#' be divided by a per-factor quantity to return it to the `h` scale.
#'
#' Which quantity is appropriate depends on how much of the source's gene space
#' the target actually measures:
#'
#' \itemize{
#'   \item When the target shares most of the source's genes, the source `d` is
#'     a good estimate of the target's per-factor scale and is used directly.
#'   \item When the target measures only a small panel, the source `d` no longer
#'     describes the target: on a 266-gene Xenium panel against a 28,916-gene
#'     Visium factorisation, `cor(log(d_source), log(d_target))` is about 0.2
#'     with per-factor ratios spanning several hundred fold, against about 0.66
#'     and roughly ten fold for a Visium target. In that regime the scale is
#'     estimated from the target itself as `rowMeans(h_new) * n_spots`, which is
#'     what `d` measures on the source expressed as a mean per observation. The
#'     mean matters: a total would make the scale depend on the number of cells
#'     in the section, which collapses the predictions for large sections.
#' }
#'
#' @param source A SingleCellExperiment or SpatialExperiment used to fit `nmf_model`.
#' @param target The object to project.
#' @param assay Assay to project.
#' @param nmf_model An NMF model from `singlet` (components `w`, `d`, `h`).
#' @param d_scale Optional length-k vector of per-factor scales to divide by. When
#'   supplied it overrides the rule above. `transfer_labels()` does not set it:
#'   pooling the target-side estimate across targets shrinks the projections
#'   below the scale of the source factors the model was fitted on, which
#'   collapses the predictions (see the note in `transfer_labels.list`).
#' @param overlap_threshold Fraction of the source's genes that the target must
#'   measure for the source `d` to be used. Defaults to 0.5.
#'
#' @return A cells x factors matrix on the source factor scale.
project_factors <- function(source, target, assay, nmf_model,
                            d_scale = NULL, overlap_threshold = 0.5){
  pr <- project_raw(source, target, assay, nmf_model)
  d_use <- if (!is.null(d_scale)) d_scale else
    target_factor_scale(pr$proj, nmf_model, pr$n_shared, overlap_threshold)
  factors <- t(pr$proj / d_use)
  colnames(factors) <- paste0("NMF", 1:ncol(factors))
  factors
}

#' Project a target onto the source loadings without rescaling.
#'
#' Returns the raw `k x n` output of `RcppML::project` together with the number
#' of genes shared with the source, so that a caller can choose the per-factor
#' scale itself (for example by pooling across several targets).
#'
#' @inheritParams project_factors
#' @return A list with `proj` (k x n) and `n_shared`.
project_raw <- function(source, target, assay, nmf_model){
  if(is(target, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(target)))){
      stop("Please provide gene symbols in your target dataset as a column named 'gene_name' in rowData.")
    }
  }

  if(is(source, "SpatialExperiment")){
    if (!(assay %in% assayNames(source))){
      stop(sprintf("Assay %s not found in the source dataset. %s needs to be available in both the source and target datasets.", assay, assay))
    }

    if(!("gene_name" %in% colnames(rowData(source)))){
      stop("Please provide gene symbols in your source dataset as a column named 'gene_name' in rowData.")
    }
  }



  loadings <- nmf_model$w
  #subset to the genes shared between the source and the target
  i<-intersect(rowData(target)$gene_name, rowData(source)$gene_name) # need to check for gene names in source object later

  if(length(i) == 0){
    stop("No intersecting genes between target and source dataset.")
  }

  rownames(loadings) <- rowData(source)$gene_name
  loadings<-loadings[rownames(loadings) %in% i,]
  loadings <- loadings[unique(rownames(loadings)),] # genes may get duplicated

  target <- target[rowData(target)$gene_name %in% i, ]
  loadings<-loadings[match(rowData(target)$gene_name,rownames(loadings)),]
  #print(any(is.na(loadings)))

  A <- assay(target, assay)
  options(RcppML.threads = 0) #line below doesn't work otherwise
  proj<-RcppML::project(data=as.matrix(A), w=loadings, threads=0, L1=0, mask=NULL)

  # proj is k x n, factors in ROWS. Any rescaling by a length-k vector must be
  # applied here, before the transpose: t(proj)/d divides an n x k matrix by a
  # length-k vector, which recycles down columns and gives all but 1 in k
  # entries the wrong divisor.
  list(proj = proj, n_shared = length(i))
}

#' Per-factor scale to return a projection to the source factor scale.
#'
#' Returns the source `d` when the target measures at least `overlap_threshold`
#' of the source's genes, and otherwise a scale estimated from the target as
#' `rowMeans(h_new) * n_spots`, which puts the target's mean weight per cell on
#' the source's mean weight per spot. See [project_factors()] for why the two
#' regimes differ, and why the mean rather than the total is used.
#'
#' @param proj A k x n projection from `RcppML::project`.
#' @param nmf_model The source NMF model.
#' @param n_shared Number of genes shared between source and target.
#' @param overlap_threshold Fraction of source genes the target must measure.
#' @param pooled_d_target Optional pre-pooled `rowSums(h_new)` summed over several
#'   targets, used instead of deriving it from a single `proj`.
#'
#' @return A length-k vector of per-factor scales.
target_factor_scale <- function(proj, nmf_model, n_shared, overlap_threshold = 0.5,
                                pooled_d_target = NULL){
  overlap <- n_shared / nrow(nmf_model$w)
  if (overlap >= overlap_threshold) {
    return(nmf_model$d)
  }
  # Normalise the MEAN weight per cell, not the TOTAL over the target.
  #
  # `rowSums(proj)` sets each factor's total mass over the target to 1, so every
  # cell gets ~1/n_cells and the scale depends on how many cells the section has.
  # The source factors `t(h)` the model was fitted on average 1/n_spots per spot,
  # so a target with more cells than the source has spots lands below the scale
  # the model expects, the intercepts dominate, and the predictions collapse onto
  # the most abundant labels. Using the mean removes n_cells from the expression
  # (rowMeans * n_spots == rowSums * n_spots / n_cells), giving a target scale
  # that does not depend on section size.
  d_target <- if (!is.null(pooled_d_target)) pooled_d_target else
    rowMeans(proj) * ncol(nmf_model$h)
  # A factor with essentially no support on the target panel would otherwise be
  # divided by ~0 and blown up; fall back to the source d for those.
  tiny <- d_target < 0.01 * stats::median(d_target)
  if (any(tiny)) {
    warning(sprintf("%d factor(s) have almost no mass in the target; using the source d for them.",
                    sum(tiny)), immediate. = TRUE)
    d_target[tiny] <- nmf_model$d[tiny]
  }
  d_target
}
