mediation_effect_cal <-
  function(para_vec,
           x_1,
           x_2,
           confound_mat,
           x4_inter = TRUE,
           x5_inter = TRUE) {
    num_conf <- ncol(confound_mat)
    k <- (length(para_vec) - 9 - 3 * num_conf) / 3
    beta_0 <- para_vec[1]
    beta_1 <- para_vec[2]
    beta_2 <- para_vec[3]
    beta_3 <- para_vec[4]
    beta_4 <- para_vec[5]
    beta_5 <- para_vec[6]
    if (!x4_inter) {
      beta_4<-0
    }
    if (!x5_inter) {
      beta_5<-0
    }
    gamma_0 <- para_vec[7]
    gamma_1 <- para_vec[8]
    phi <- para_vec[9]
    delta <- para_vec[10]
    
    alpha_0_vec <- para_vec[11:(10 + k)]
    alpha_1_vec <- para_vec[(11 + k):(10 + 2 * k)]
    
    
    if (k == 1) {
      psi_vec <- 1
      beta_conf <- para_vec[seq(13, 12 + num_conf)]
      alpha_conf <-
        para_vec[seq(13 + num_conf, 12 + 2 * num_conf)]
      gamma_conf <-
        para_vec[seq(13 + 2 * num_conf, 12 + 3 * num_conf)]
    } else {
      psi_temp <- para_vec[(11 + 2 * k):(9 + 3 * k)]
      psi_vec <- c(psi_temp, 1 - sum(psi_temp))
      beta_conf <- para_vec[seq(10 + 3 * k, 9 + 3 * k + num_conf)]
      alpha_conf <-
        para_vec[seq(10 + 3 * k + num_conf, 9 + 3 * k + 2 * num_conf)]
      gamma_conf <-
        para_vec[seq(10 + 3 * k + 2 * num_conf, 9 + 3 * k + 3 * num_conf)]
    }
    
    conf_mean <- colMeans(confound_mat)
    
    
    
    delta_1_x2 <-
      gamma_0 + gamma_1 * x_2 + sum(gamma_conf * conf_mean)
    mu_x2 <-
      alpha_0_vec + alpha_1_vec * x_2 + sum(alpha_conf * conf_mean)
    # mu_2_x2<-alpha_02+alpha_12*x_2
    
    delta_1_x1 <-
      gamma_0 + gamma_1 * x_1 + sum(gamma_conf * conf_mean)
    mu_x1 <-
      alpha_0_vec + alpha_1_vec * x_1 + sum(alpha_conf * conf_mean)
    # mu_2_x1<-alpha_02+alpha_12*x_1
    
    NIE_1 <-
      (beta_1 + beta_5 * x_2) * (expect_M_x(psi_vec, delta_1_x2, mu_x2) - expect_M_x(psi_vec, delta_1_x1, mu_x1))
    NIE_2 <-
      (beta_2 + beta_4 * x_2) * (expit(delta_1_x1) - expit(delta_1_x2))
    NDE <-
      beta_3 * (x_2 - x_1) + beta_4 * (x_2 - x_1) * (1 - expit(delta_1_x1)) + beta_5 * (x_2 - x_1) * expect_M_x(psi_vec, delta_1_x1, mu_x1)
    NIE <- NIE_1 + NIE_2
    
    return(c(NIE_1, NIE_2, NDE, NIE))
  }
