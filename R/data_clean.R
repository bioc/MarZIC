data_clean <- function(MicrobData,
                       CovData,
                       lib_name,
                       y_name,
                       x_name,
                       conf_name,
                       taxa_of_interest,
                       taxDropThresh,
                       taxDropCount,
                       SDThresh,
                       SDbinThresh,
                       SDx,
                       SDy,
                       transfer_to_RA) {
  cov_names <- colnames(CovData)
  if (!(x_name %in% cov_names)) {
    stop(x_name, " is not available in dataset. Please double check.")
  }
  if (!(y_name %in% cov_names)) {
    stop(y_name, " is not available in dataset. Please double check.")
  }
  if (!(lib_name %in% cov_names)) {
    stop(lib_name, " is not available in dataset. Please double check.")
  }

  if (any(!(conf_name %in% cov_names))) {
    stop("One or more confounder is not available in dataset. Please double check.")
  }

  if (any(MicrobData > 1) && transfer_to_RA == FALSE) {
    stop(
      "One or more taxon reads are larger than 1. It's definitely not Relative abundance.
         Please double check."
    )
  }

  taxon_name <- colnames(MicrobData)

  if (any(is.na(MicrobData), is.na(CovData))) {
    message("One or more reads are missing. Subjects with missing were removed from analysis")
    complete_dat<-na.omit(cbind(MicrobData,CovData))
    MicrobData <- complete_dat[,taxon_name,drop=FALSE]
    CovData<-complete_dat[,cov_names,drop=FALSE]
  }

  if (transfer_to_RA) {
    MicrobData <- apply(MicrobData, 2, function(x) {
      x / rowSums(MicrobData)
    })
  }

  if (any(taxa_of_interest != "all")) {
    if (!all(taxa_of_interest %in% colnames(MicrobData))) {
      stop("one or more taxon in taxa_of_interest is not presented in data")
    }
    MicrobData <-
      MicrobData[, colnames(MicrobData) %in% taxa_of_interest, drop = FALSE]
    taxon_name <- colnames(MicrobData)
  }


  xi_vec <- CovData[, x_name]
  yi_vec <- CovData[, y_name]
  conf_mat <- as.matrix(CovData[, conf_name, drop = FALSE])


  SDx_real <- sd(xi_vec) / mean(xi_vec)
  SDy_real <- sd(yi_vec) / mean(yi_vec)

  nSub <- nrow(MicrobData)
  nTaxa <- ncol(MicrobData)


  if (length(taxon_name) < nTaxa) {
    warning("Taxon name is not available, or less than number of taxon.
            Taxa are renamed as rawCounts")
    colnames(MicrobData) <- paste0("rawCount", seq_len(nTaxa))
    taxon_name <- colnames(MicrobData)
  }

  if (abs(SDx_real) < SDx) {
    stop(x_name, " is nearly constant. Please double check.")
  }

  if (abs(SDy_real) < SDy) {
    stop(y_name, " is nearly constant. Please double check.")
  }


  too_much_zero_exc_name <- c()
  too_low_SD_exc_name <- c()
  # too_unbalanced_bin_exc_name<-c()
  # bin_x_ind <- table(xi_vec)==2



  for (i in seq_len(ncol(MicrobData))) {
    taxon_tobe_test <- MicrobData[, taxon_name[i]]

    taxon_sparsity <- sum(taxon_tobe_test == 0) / nSub
    if (all(taxon_sparsity > taxDropThresh, sum(taxon_tobe_test > 0) < taxDropCount)) {
      too_much_zero_exc_name <- c(too_much_zero_exc_name, taxon_name[i])
    }
    sdT <- sd(taxon_tobe_test) / mean(taxon_tobe_test)
    if (abs(sdT) < SDThresh) {
      too_low_SD_exc_name <- c(too_low_SD_exc_name, taxon_name[i])
    }
  }
  if (any(
    length(too_much_zero_exc_name) == nTaxa,
    length(too_low_SD_exc_name) == nTaxa
  )) {
    stop("All taxa are either too many zeros or nearly constant. Please double check")
  }
  if (length(too_much_zero_exc_name) > 0) {
    too_much_zero_exc_name_mes <- paste0(too_much_zero_exc_name, " ")
    warning(
      too_much_zero_exc_name_mes,
      "were excluded from analysis due to too many zeros"
    )
  }
  if (length(too_low_SD_exc_name) > 0) {
    too_low_SD_exc_name_mes <- paste0(too_low_SD_exc_name, " ")
    warning(
      too_low_SD_exc_name_mes,
      "were excluded from analysis due to nearly constant"
    )
  }

  exclude_taxon_name <-
    unique(c(too_much_zero_exc_name, too_low_SD_exc_name))
  MicrobData_clean <-
    MicrobData[, !(taxon_name %in% exclude_taxon_name), drop = FALSE]

  conf_sd <- apply(conf_mat, 2, sd) / colMeans(conf_mat)
  exclude_conf_name <- conf_name[abs(conf_sd) <= SDx]
  if (length(exclude_conf_name) > 0) {
    exclude_conf_name_mes <- paste0(exclude_conf_name, " ")
    warning(
      exclude_conf_name_mes,
      "were exluded from confounders due to nearly constant"
    )
  }

  conf_name_remain <- conf_name[abs(conf_sd) > SDx]
  return(list(
    MicrobData_clean = MicrobData_clean,
    conf_name_remain = conf_name_remain,
    CovData = CovData
  ))
}
