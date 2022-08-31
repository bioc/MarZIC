##' Marginal Mediation Model for Zero-Inflated Compositional Mediators
##'
##' @description
##' \loadmathjax
##' MarZICM is used for calculating mediation effects specifically for zero-inflated compositional
##' mediators.
##' The marginal outcome model for taxon \mjeqn{j}{} is:
##' \mjdeqn{Y=\beta_0+\beta_1M_j+\beta_21_{M_j>0}+\beta_3X+\beta_4X1_{M_j>0}+\beta_5XM_j+\epsilon}{}
##' where \mjeqn{1_{()}}{} is indicator function
##' X is the covariate of interest
##' The probability of \mjeqn{M_j}{} being structure zero (\mjeqn{\Delta}{}) is:
##' \mjdeqn{\log(\frac{\Delta_j}{1-\Delta_j})=\gamma_0 + \gamma_1X}{}
##' The mean of \mjeqn{M_j}{} in compositional structure is
##' \mjdeqn{\log(\frac{\mu_j}{1-\mu_j})=\alpha_0 + \alpha_1X}{}
##'
##'
##' @param Experiment_dat A SummarizedExperiment object containing microbiome data as assay and
##' covariates, outcome and library size as colData. The microbiome data could be relative abundance or absolute
##' abundance. Missing value should be imputed in advance.
##' @param lib_name      Name of library size variable within colData.
##' @param y_name   Name of outcome variable within colData.
##' @param x_name Name of covariate of interest within colData.
##' @param conf_name Name of confounders within colData. Defaule is NULL, meaning no confounder.
##' @param mediator_mix_range Number of mixtures in mediator. Default is 1, meaning no mixture.
##' @param transfer_to_RA Logical variable indicating whether the microbiome data should be
##' transferred to relative abundance. Default is FALSE. If TRUE, microbiome data will be rescaled
##' by its row sum.
##' @param num_cores Number of CPU cores to be used in parallelization task.
##' @param adjust_method P value adjustment method. Same as p.adjust. Default is "fdr".
##' @param fdr_rate FDR cutoff for significance. Default is 0.2.
##' @param taxDropThresh The threshold of dropping taxon due to high zero percentage. Default is
##' 0.9, meaning taxon will be dropped for analysis if zero percentage is higher than 90\%.
##' @param taxDropCount The threshold of dropping taxon due to not enough non-zero observation counts.
##' Default is 20, meaning taxon will be dropped if non-zero observation is less than 20.
##' @param SDThresh The threshold of dropping taxon due to low coefficient of variation (CV)
##' to avoid constant taxon.
##' Default is 0.05, meaning any taxon has CV less than 0.05 will be dropped.
##' @param SDx The threshold of stopping analysis due to low CV of covariate of interest.
##' Default is 0.05, meaning when CV of covariate of interest is less than 0.05, the analysis
##' will be stopped.
##' @param SDy The threshold of stopping analysis due to low CV of outcome.
##' Default is 0.05, meaning when CV of outcome. is less than 0.05, the analysis
##' will be stopped.
##'
##' @return
##' A `list` of `4` datasets containing the results for `NIE1`, `NIE2`, `NDE`, and `NIE`.
##' Each dataset has row representing each taxon, 6 columns for `Estimates`, `Standard Error`,
##' `Lower bound for 95% Confidence Interval`, `Upper bound for 95% Confidence Interval`,
##' `Adjusted p value`, `Significance indicator`.
##'
##' @examples {
##' library(MarZICM)
##' library(SummarizedExperiment)
##' ## A make up example with 1 taxon and 100 subjects.
##' set.seed(1)
##' nSub<-100
##' ## generate covariate of interest X
##' X <- rbinom(nSub,1,0.5)
##' ## generate mean of mediator M
##' mu <- exp(X)/(1+exp(X))
##' phi <- 10
##' ## generate mediator M
##' M<-rbeta(nSub,mu*phi,(1-mu)*phi)
##' ## library size set to 10000
##' libsize <- 10000
##' ## generate observed zero, both structure zero and zero due to LOD mechanism
##' non_zero_ind <- 1 - rbinom(nSub,1,exp(0.3*X)/(1+exp(0.3*X)))
##' obs_m<-M * ((M * libsize)>1) * non_zero_ind
##' ## generate outcome Y
##' Y <- 1 + 100 * M + (M > 0) + X
##' ## Construct SummerizedExperiment object
##' CovData <- cbind(Y=Y,X=X,libsize=libsize)
##' test_dat <- SummarizedExperiment(assays=list(MicrobData=t(M)),colData=CovData)
##' res<-MarZICM(Experiment_dat = test_dat,
##' lib_name = "libsize",
##' y_name = "Y",
##' x_name = "X",
##' num_cores = 1,
##' mediator_mix_range = 1)
##'
##' ## Pull out significant NIE1
##' NIE1 <- res$NIE1
##' subset(NIE1,significance == TRUE)
##' }
##'
##' @importFrom foreach foreach %dopar% registerDoSEQ
##' @importFrom parallel makeCluster clusterExport stopCluster clusterSetRNGStream detectCores
##' @importFrom doParallel registerDoParallel
##' @importFrom NlcOptim solnl
##' @importFrom betareg betareg
##' @importFrom SummarizedExperiment assays colData SummarizedExperiment
##' @importFrom S4Vectors DataFrame
##' @importFrom pracma jacobian grad hessian
##' @import Rcpp
##' @import stats
##' @import mathjaxr
##'
##'
##'
##' @export


MarZICM <- function(Experiment_dat,
                    lib_name,
                    y_name,
                    x_name,
                    conf_name = NULL,
                    mediator_mix_range = 1,
                    transfer_to_RA = FALSE,
                    num_cores = detectCores() - 2,
                    adjust_method = "fdr",
                    fdr_rate = 0.2,
                    taxDropThresh = 0.8,
                    taxDropCount = 20,
                    SDThresh = 0.05,
                    SDx = 0.05,
                    SDy = 0.05) {
  assay_name <- names(assays(Experiment_dat))
  MicrobData <- t(assays(Experiment_dat)[[assay_name]])
  CovData <- colData(Experiment_dat)

  clean_dat <- data_clean(
    MicrobData = MicrobData,
    CovData = CovData,
    lib_name = lib_name,
    y_name = y_name,
    x_name = x_name,
    conf_name = conf_name,
    taxDropThresh = taxDropThresh,
    taxDropCount = taxDropCount,
    SDThresh = SDThresh,
    SDx = SDx,
    SDy = SDy,
    transfer_to_RA = transfer_to_RA
  )

  MicrobData_clean <- clean_dat$MicrobData_clean
  conf_name_remain <- clean_dat$conf_name_remain
  res_list <- suppressWarnings(
    apply_real_data_func(
      MicrobData = MicrobData_clean,
      CovData = CovData,
      lib_name = lib_name,
      y_name = y_name,
      x_name = x_name,
      conf_name = conf_name_remain,
      k_range = mediator_mix_range,
      num_cores = num_cores
    )
  )

  nTaxa <- res_list$nTaxa
  nSub <- res_list$nSub

  NIE1_save <- DataFrame(matrix(nrow = nTaxa, ncol = 6))
  NIE2_save <- DataFrame(matrix(nrow = nTaxa, ncol = 6))
  NDE_save <- DataFrame(matrix(nrow = nTaxa, ncol = 6))
  NIE_save <- DataFrame(matrix(nrow = nTaxa, ncol = 6))

  rownames(NIE1_save) <-
    rownames(NIE2_save) <-
    rownames(NDE_save) <- rownames(NIE_save) <- res_list$taxon_ori_name
  colnames(NIE1_save) <- colnames(NIE2_save) <- colnames(NDE_save) <- colnames(NIE_save) <-
    c("est", "se", "CI low", "CI up", "p value adj", "significance")
  for (i in seq_len(length(res_list$list_save))) {
    if (is.na(res_list$list_save[[i]]$res_fin_med)[1]) {
      NIE1_save[i, ] <- NIE2_save[i, ] <- NDE_save[i, ] <- NIE_save[i, ] <- NA
    } else {
      res_temp <- res_list$list_save[[i]]$res_fin_med
      NIE1_save[i, 1] <- res_temp$mediation_effect[1]
      NIE2_save[i, 1] <- res_temp$mediation_effect[2]
      NDE_save[i, 1] <- res_temp$mediation_effect[3]
      NIE_save[i, 1] <- res_temp$mediation_effect[4]

      NIE1_save[i, 2] <- res_temp$NIE_sd[1]
      NIE2_save[i, 2] <- res_temp$NIE_sd[2]
      NDE_save[i, 2] <- res_temp$NIE_sd[3]
      NIE_save[i, 2] <- res_temp$NIE_sd[4]

      NIE1_save[i, 3] <- res_temp$mediation_effect[1] - 1.96 * res_temp$NIE_sd[1]
      NIE2_save[i, 3] <- res_temp$mediation_effect[2] - 1.96 * res_temp$NIE_sd[2]
      NDE_save[i, 3] <- res_temp$mediation_effect[3] - 1.96 * res_temp$NIE_sd[3]
      NIE_save[i, 3] <- res_temp$mediation_effect[4] - 1.96 * res_temp$NIE_sd[4]

      NIE1_save[i, 4] <- res_temp$mediation_effect[1] + 1.96 * res_temp$NIE_sd[1]
      NIE2_save[i, 4] <- res_temp$mediation_effect[2] + 1.96 * res_temp$NIE_sd[2]
      NDE_save[i, 4] <- res_temp$mediation_effect[3] + 1.96 * res_temp$NIE_sd[3]
      NIE_save[i, 4] <- res_temp$mediation_effect[4] + 1.96 * res_temp$NIE_sd[4]
    }
  }

  NIE1_save[, 5] <- p.adjust((1 - pnorm(abs(NIE1_save[, 1] / NIE1_save[, 2]))) * 2, adjust_method)
  NIE2_save[, 5] <- p.adjust((1 - pnorm(abs(NIE2_save[, 1] / NIE2_save[, 2]))) * 2, adjust_method)
  NDE_save[, 5] <- p.adjust((1 - pnorm(abs(NDE_save[, 1] / NDE_save[, 2]))) * 2, adjust_method)
  NIE_save[, 5] <- p.adjust((1 - pnorm(abs(NIE_save[, 1] / NIE_save[, 2]))) * 2, adjust_method)

  NIE1_save[, 6] <- NIE1_save[, 5] < fdr_rate
  NIE2_save[, 6] <- NIE2_save[, 5] < fdr_rate
  NDE_save[, 6] <- NDE_save[, 5] < fdr_rate
  NIE_save[, 6] <- NIE_save[, 5] < fdr_rate

  output_list <- list(
    NIE1_save = NIE1_save,
    NIE2_save = NIE2_save,
    NDE_save = NDE_save,
    NIE_save = NIE_save
  )

  return(output_list)
}
