apply_real_data_func <-
  function(MicrobData,
           CovData,
           lib_name,
           y_name,
           x_name,
           conf_name,
           k_range,
           num_cores) {
    num_taxon <- ncol(MicrobData)
    num_sub <- nrow(MicrobData)
    taxon_ori_name <- colnames(MicrobData)

    yi_vec <- CovData[, y_name]
    xi_vec <- CovData[, x_name]
    li_vec <- CovData[, lib_name]
    conf_mat <- as.matrix(CovData[, conf_name, drop = FALSE])

    trial <- seq_len(num_taxon)

    # cl <- parallel::makeCluster(num_cores)
    #
    # parallel::clusterExport(cl = cl,
    #                         varlist = c("real_data_run_func_nomix",
    #                                     "Q_theta_cpp_nomix"),
    #                         envir = environment())

    registerDoParallel(num_cores)

    i <- numeric(0)

    list_save <- foreach(i = trial) %dopar% {
      obs_m_vec <- MicrobData[, i]
      temp_name <- colnames(MicrobData)[i]
      res_list_med <- list()
      AIC_select <- c()
      for (k in k_range) {
        if (k == 1) {
          if (sum(obs_m_vec == 0) > (0.1 * length(obs_m_vec))) {
            res_temp <-
              try(
                real_data_run_func_nomix(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, k),
                T
              )
          } else {
            res_temp <-
              try(
                real_data_run_func_nz_nomix(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, k),
                T
              )
          }
        } else {
          if (sum(obs_m_vec == 0) > (0.1 * length(obs_m_vec))) {
            res_temp <-
              try(
                real_data_run_func(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, k),
                T
              )
          } else {
            res_temp <-
              try(
                real_data_run_func_nz(yi_vec, obs_m_vec, xi_vec, li_vec, conf_mat, k),
                T
              )
          }
        }
        if (inherits(res_temp, "try-error")) {
          AIC_select[k] <- NA
        } else if (any(is.na(res_temp$mediation_effect))) {
          AIC_select[k] <- NA
        } else {
          AIC_select[k] <- res_temp$AIC_est
        }
        res_temp$taxon_name <- temp_name
        res_list_med[[k]] <- res_temp
      }
      if (all(is.na(AIC_select))) {
        res_fin_med <- NA
      } else {
        res_fin_med <- res_list_med[[which.min(AIC_select)]]
      }

      return(list(res_fin_med = res_fin_med, res_list_med = res_list_med))
    }
    # parallel::stopCluster(cl)


    # list_save<-mclapply(trial,mcapply_func,mc.cores = num_cores)

    return(
      list(
        list_save = list_save,
        nTaxa = num_taxon,
        nSub = num_sub,
        taxon_ori_name = taxon_ori_name
      )
    )
  }
