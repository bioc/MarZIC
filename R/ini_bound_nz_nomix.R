ini_bound_nz_nomix <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  # ind_M<-as.numeric(m_star_vec>0)
  Y_mod <- glm(yi_vec ~ m_star_vec + x_i_vec + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[c(1, 2, 4, 6)] <- Y_mod_ini
  ini_par[c(3, 5)] <- 0

  ### kmeans to group
  m_nz <- m_star_vec[m_star_vec > 0]
  x_i_nz <- x_i_vec[m_star_vec > 0]
  beta_mod <- try(betareg(m_nz ~ x_i_nz, link = "logit"), TRUE)
  error_ind <- inherits(beta_mod, "try-error")
  if (error_ind) {
    alpha_0 <- mean(m_nz)
    alpha_1 <- 0
    phi_est <- NA
  } else {
    beta_gp <- as.numeric(beta_mod$coefficients$mean)
    alpha_0 <- beta_gp[1]
    alpha_1 <- beta_gp[2]
    phi_est <- as.numeric(beta_mod$coefficients$precision)
  }
  # M_mod<-glm((1-ind_M)~x_i_vec,family = "binomial")
  Del_est <- c(-100, 0)

  ini_par[7:8] <- Del_est
  # phi_est<-mean(phi_est_vec,na.rm = T)
  if (is.na(phi_est)) {
    phi_est <- 10
  }
  ini_par[9] <- min(phi_est, 100)
  ini_par[10] <- sqrt(summary(Y_mod)$dispersion)

  ini_par[11] <- alpha_0
  ini_par[12] <- alpha_1

  lb_est <- c(-Inf, -Inf, 0, -Inf, 0, -Inf, -100, 0, 0.1, 0.5, -10, -10)
  ub_est <- c(Inf, Inf, 0, Inf, 0, Inf, -100, 0, Inf, Inf, 10, 10)


  return(list(ini_par, lb_est, ub_est))
}
