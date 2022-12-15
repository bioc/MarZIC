##' Marginal Mediation Model for Zero-Inflated Compositional Mediators
##'
##' @description
##' \loadmathjax
##' MarZIC is used for calculating marginal mediation effects for zero-inflated compositional
##' mediators. For microbiome data, the marginal outcome model for the \mjeqn{j}{}th taxon (or OTU, ASV) is:
##' \mjdeqn{Y=\beta_0+\beta_1M_j+\beta_21_{M_j>0}+\beta_3X+\beta_4X1_{M_j>0}+\beta_5XM_j+\epsilon}{}
##' where \mjeqn{1_{()}}{} is indicator function, 
##' X is the covariate of interest 
##' and \mjeqn{M_j}{} is the relative abundance of the \mjeqn{j}{}th taxon. 
##' The probability of \mjeqn{M_j}{} being structure zero (ie, true zeros) is:
##' \mjdeqn{\log(\frac{\Delta_j}{1-\Delta_j})=\gamma_0 + \gamma_1X}{}
##' The mean of \mjeqn{M_j}{} in compositional structure is modeled as:
##' \mjdeqn{\log(\frac{\mu_j}{1-\mu_j})=\alpha_0 + \alpha_1X}{}
##' Typically, users just need to feed the first six inputs to the function: 
##'`Experiment_dat`, `lib_name`, `y_name`, `x_name`, `conf_name` and `taxa_of_interest`.
##'
##' @param Experiment_dat A SummarizedExperiment object containing microbiome data as assay and
##' covariates, outcome and library size as colData. The microbiome data could be relative abundance or absolute
##' abundance. Subjects with missing value will be removed during analysis.
##' @param lib_name      Name of library size variable within colData.
##' @param y_name   Name of outcome variable within colData.
##' @param x_name Name of covariate of interest within colData.
##' @param conf_name Name of confounders within colData. Defaule is NULL, meaning no confounder.
##' @param x4_inter Whether to include the interaction term \mjeqn{\beta_4}{}. Default is TRUE.
##' @param x5_inter Whether to include the interaction term \mjeqn{\beta_5}{}. Default is TRUE.
##' @param taxa_of_interest A character vector for taxa names indicating taxa that should be analyzed. Default
##' is "all", meaning all taxa should be included into analysis.
##' @param mediator_mix_range Number of mixtures in mediator. Default is 1, meaning no mixture.
##' @param transfer_to_RA Logical variable indicating whether the microbiome data should be
##' transferred to relative abundance. Default is TRUE. If TRUE, microbiome data will be rescaled
##' by its row sum.
##' @param num_cores Number of CPU cores to be used in parallelization task.
##' @param adjust_method P value adjustment method. Same as p.adjust. Default is "fdr".
##' @param fdr_rate FDR cutoff for significance. Default is 0.05.
##' @param taxDropThresh The threshold of dropping taxon due to high zero percentage. Default is
##' 0.9, meaning taxon will be dropped for analysis if zero percentage is higher than 90\%.
##' @param taxDropCount The threshold of dropping taxon due to not enough non-zero observation counts.
##' Default is 4 * (length(conf_name)+2), meaning taxon will be dropped if non-zero observation is less than four times
##' of number of covariates plus 1.
##' @param zero_prop_NIE2 The threshold of zero percentage for calculating NIE2. Default is 0.1,
##' meaning NIE2 will be calculated for taxon with zero percentage greater than 10\%.
##' @param zero_count_NIE2 The threshold of zero counts for calculating NIE2.
##' Default is 4 * (length(conf_name)+2), meaning NIE2 will be calculated for taxon with zero counts
##' greater than four times of number of covariates plus 1.
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
##' library(MarZIC)
##'
##' ## A make up example with 10 taxon and 200 subjects.
##' set.seed(1)
##' nSub <- 200
##' nTaxa <- 10
##' ## generate covariate of interest X
##' X <- rbinom(nSub, 1, 0.5)
##' ## generate confounders
##' conf1<-rnorm(nSub)
##' conf2<-rbinom(nSub,1,0.5)
##' ## generate mean of each taxon. All taxon are having the same mean for simplicity.
##' mu <- exp(-5 + X + 0.1 * conf1 + 0.1 * conf2) /
##'  (1 + exp(-5 + X + 0.1 * conf1 + 0.1 * conf2))
##' phi <- 10
##'
##' ## generate true RA
##' M_taxon<-t(sapply(mu,function(x) dirmult::rdirichlet(n=1,rep(x*phi,nTaxa))))
##'
##' P_zero <- exp(-3 + 0.3 * X + 0.1 * conf1 + 0.1 * conf2) /
##'  (1 + exp(-3 + 0.3 * X + 0.1 * conf1 + 0.1 * conf2))
##'
##' non_zero_ind <- t(sapply(P_zero,function(x) 1-rbinom(nTaxa,1,rep(x,nTaxa))))
##'
##' True_RA<-t(apply(M_taxon*non_zero_ind,1,function(x) x/sum(x)))
##'
##' ## generate outcome Y based on true RA
##' Y <- 1 + 100 * True_RA[,1] + 5 * (True_RA[,1] > 0) + X + conf1 + conf2 + rnorm(nSub)
##'
##' ## library size was set to 10,000 for all subjects for simplicity.
##' libsize <- 10000
##'
##' ## generate observed RA
##' observed_AA <- floor(M_taxon*libsize*non_zero_ind)
##'
##' observed_RA <- t(apply(observed_AA,1,function(x) x/sum(x)))
##' colnames(observed_RA)<-paste0("rawCount",seq_len(nTaxa))
##'
##'
##'
##' ## Construct SummerizedExperiment object
##' CovData <- cbind(Y = Y, X = X, libsize = libsize, conf1 = conf1, conf2 = conf2)
##' test_dat <-
##'   SummarizedExperiment::SummarizedExperiment(assays = list(MicrobData = t(observed_RA)), 
##'   colData = CovData)
##' ## run the analysis
##' res <- MarZIC(
##'   Experiment_dat = test_dat,
##'   lib_name = "libsize",
##'   y_name = "Y",
##'   x_name = "X",
##'   conf_name = c("conf1","conf2"),
##'   taxa_of_interest = c("rawCount1","rawCount2","rawCount3"),
##'   num_cores = 1,
##'   mediator_mix_range = 1
##' )
##' }
##' @references Wu et al.(2022) MarZIC: A Marginal Mediation Model for Zero-Inflated Compositional Mediators with Applications to Microbiome Data. Genes 2022, 13, 1049.
##'
##'
##' @importFrom foreach foreach %dopar% registerDoSEQ
##' @importFrom parallel makeCluster clusterExport stopCluster clusterSetRNGStream detectCores
##' @importFrom doParallel registerDoParallel
##' @importFrom NlcOptim solnl
##' @importFrom betareg betareg
##' @importFrom SummarizedExperiment assays colData SummarizedExperiment
##' @importFrom S4Vectors DataFrame
##' @importFrom pracma jacobian grad hessian
##' @importFrom dirmult rdirichlet
##' @import Rcpp
##' @import stats
##' @import mathjaxr
##'
##'
##'
##' @export
##' @useDynLib MarZIC, .registration=TRUE


MarZIC <- function(Experiment_dat,
                   lib_name,
                   y_name,
                   x_name,
                   conf_name = NULL,
                   x4_inter = TRUE,
                   x5_inter = TRUE,
                   taxa_of_interest = "all",
                   mediator_mix_range = 1,
                   transfer_to_RA = TRUE,
                   num_cores = max(detectCores() - 2,1),
                   adjust_method = "fdr",
                   fdr_rate = 0.05,
                   taxDropThresh = 0.8,
                   taxDropCount = 4 * (length(conf_name) + 2),
                   zero_prop_NIE2 = 0.1,
                   zero_count_NIE2 = 4 * (length(conf_name) + 2),
                   SDThresh = 0.05,
                   SDx = 0.05,
                   SDy = 0.05) {
  assay_name <- names(assays(Experiment_dat))
  MicrobData <- t(assays(Experiment_dat)[[assay_name]])
  CovData <- as.data.frame(colData(Experiment_dat))
  
  clean_dat <- data_clean(
    MicrobData = MicrobData,
    CovData = CovData,
    lib_name = lib_name,
    y_name = y_name,
    x_name = x_name,
    conf_name = conf_name,
    taxa_of_interest = taxa_of_interest,
    taxDropThresh = taxDropThresh,
    taxDropCount = taxDropCount,
    SDThresh = SDThresh,
    SDx = SDx,
    SDy = SDy,
    transfer_to_RA = transfer_to_RA
  )
  
  MicrobData_clean <- clean_dat$MicrobData_clean
  conf_name_remain <- clean_dat$conf_name_remain
  CovData <- clean_dat$CovData
  res_list <- suppressWarnings(
    apply_real_data_func(
      MicrobData = MicrobData_clean,
      CovData = CovData,
      lib_name = lib_name,
      y_name = y_name,
      x_name = x_name,
      conf_name = conf_name_remain,
      k_range = mediator_mix_range,
      num_cores = num_cores,
      zero_prop_NIE2 = zero_prop_NIE2,
      zero_count_NIE2 = zero_count_NIE2,
      x4_inter = x4_inter,
      x5_inter = x5_inter
    )
  )
  
  nTaxa <- res_list$nTaxa
  nSub <- res_list$nSub
  
  NIE1_save <- DataFrame(matrix(nrow = nTaxa, ncol = 7))
  NIE2_save <- DataFrame(matrix(nrow = nTaxa, ncol = 7))
  NDE_save <- DataFrame(matrix(nrow = nTaxa, ncol = 7))
  NIE_save <- DataFrame(matrix(nrow = nTaxa, ncol = 7))
  
  rownames(NIE1_save) <-
    rownames(NIE2_save) <-
    rownames(NDE_save) <-
    rownames(NIE_save) <- res_list$taxon_ori_name
  colnames(NIE1_save) <-
    colnames(NIE2_save) <- colnames(NDE_save) <- colnames(NIE_save) <-
    c("est",
      "se",
      "CI low",
      "CI up",
      "p value unadj" ,
      "p value adj",
      "significance")
  for (i in seq_len(length(res_list$list_save))) {
    if (is.na(res_list$list_save[[i]]$res_fin_med)[1]) {
      NIE1_save[i,] <-
        NIE2_save[i,] <- NDE_save[i,] <- NIE_save[i,] <- NA
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
      
      NIE1_save[i, 3] <-
        res_temp$mediation_effect[1] - 1.96 * res_temp$NIE_sd[1]
      NIE2_save[i, 3] <-
        res_temp$mediation_effect[2] - 1.96 * res_temp$NIE_sd[2]
      NDE_save[i, 3] <-
        res_temp$mediation_effect[3] - 1.96 * res_temp$NIE_sd[3]
      NIE_save[i, 3] <-
        res_temp$mediation_effect[4] - 1.96 * res_temp$NIE_sd[4]
      
      NIE1_save[i, 4] <-
        res_temp$mediation_effect[1] + 1.96 * res_temp$NIE_sd[1]
      NIE2_save[i, 4] <-
        res_temp$mediation_effect[2] + 1.96 * res_temp$NIE_sd[2]
      NDE_save[i, 4] <-
        res_temp$mediation_effect[3] + 1.96 * res_temp$NIE_sd[3]
      NIE_save[i, 4] <-
        res_temp$mediation_effect[4] + 1.96 * res_temp$NIE_sd[4]
    }
  }
  
  NIE1_save[, 5] <-
    (1 - pnorm(abs(NIE1_save[, 1] / NIE1_save[, 2]))) * 2
  NIE2_save[, 5] <-
    (1 - pnorm(abs(NIE2_save[, 1] / NIE2_save[, 2]))) * 2
  NDE_save[, 5] <-
    (1 - pnorm(abs(NDE_save[, 1] / NDE_save[, 2]))) * 2
  NIE_save[, 5] <-
    (1 - pnorm(abs(NIE_save[, 1] / NIE_save[, 2]))) * 2
  
  NIE1_save[, 6] <-
    p.adjust((1 - pnorm(abs(
      NIE1_save[, 1] / NIE1_save[, 2]
    ))) * 2, adjust_method)
  NIE2_save[, 6] <-
    p.adjust((1 - pnorm(abs(
      NIE2_save[, 1] / NIE2_save[, 2]
    ))) * 2, adjust_method)
  NDE_save[, 6] <-
    p.adjust((1 - pnorm(abs(
      NDE_save[, 1] / NDE_save[, 2]
    ))) * 2, adjust_method)
  NIE_save[, 6] <-
    p.adjust((1 - pnorm(abs(
      NIE_save[, 1] / NIE_save[, 2]
    ))) * 2, adjust_method)
  
  NIE1_save[, 7] <- NIE1_save[, 6] < fdr_rate
  NIE2_save[, 7] <- NIE2_save[, 6] < fdr_rate
  NDE_save[, 7] <- NDE_save[, 6] < fdr_rate
  NIE_save[, 7] <- NIE_save[, 6] < fdr_rate
  
  output_list <- list(
    NIE1_save = NIE1_save,
    NIE2_save = NIE2_save,
    NDE_save = NDE_save,
    NIE_save = NIE_save
  )
  
  return(output_list)
}
