ini_bound <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  ind_M <- as.numeric(m_star_vec > 0)
  Y_mod <- glm(yi_vec ~ m_star_vec + ind_M + x_i_vec + x_i_vec * ind_M + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[1:6] <- Y_mod_ini

  ### kmeans to group
  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  km_res <- kmeans(m_nz, k)
  gp_res <- km_res$cluster
  alpha_0_vec <- numeric(k)
  alpha_1_vec <- numeric(k)
  psi_est_vec <- numeric(k)
  phi_est_vec <- numeric(k)

  for (i in 1:k) {
    m_gp <- m_nz[gp_res == i]
    x_gp <- x_i_nz[gp_res == i]

    beta_mod <- try(betareg(m_gp ~ x_gp, link = "logit"), TRUE)
    error_ind <- inherits(beta_mod, "try-error")
    if (error_ind) {
      alpha_0_vec[i] <- mean(m_gp)
      alpha_1_vec[i] <- 0
      phi_est_vec[i] <- NA
    } else {
      beta_gp <- as.numeric(beta_mod$coefficients$mean)
      alpha_0_vec[i] <- beta_gp[1]
      alpha_1_vec[i] <- beta_gp[2]
      phi_est_vec[i] <- as.numeric(beta_mod$coefficients$precision)
    }

    psi_est_vec[i] <- length(m_gp) / length(m_nz)
    # beta_mod2<-betareg(m_gp2~x_gp2,link = "logit")
    # beta_gp2<-as.numeric(beta_mod2$coefficients$mean)
  }

  alpha_0_sort <- order(alpha_0_vec, decreasing = T)
  alpha_0_sorted <- alpha_0_vec[alpha_0_sort]
  alpha_1_sorted <- alpha_1_vec[alpha_0_sort]
  psi_est_sorted <- psi_est_vec[alpha_0_sort]


  M_mod <- glm((1 - ind_M) ~ x_i_vec, family = "binomial")
  Del_est <- as.numeric(M_mod$coefficients)
  ini_par[7:8] <- Del_est
  phi_est <- mean(phi_est_vec, na.rm = T)
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11:(10 + k)] <- alpha_0_sorted
  ini_par[(11 + k):(10 + 2 * k)] <- alpha_1_sorted
  ini_par[(11 + 2 * k):(9 + 3 * k)] <- psi_est_sorted[1:(k - 1)]



  lb_est <- c(rep(-Inf, 6), rep(-10, 2), 5, 0.5, rep(-10, k), rep(-10, k), rep(0.01, k - 1))
  ub_est <- c(rep(Inf, 6), rep(10, 2), Inf, Inf, rep(10, k), rep(10, k), rep(0.99, k - 1))


  return(list(ini_par, lb_est, ub_est))
}
