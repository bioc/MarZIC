ini_bound_nz <- function(yi_vec, m_star_vec, x_i_vec, k) {
  ini_par <- numeric(9 + 3 * k)
  # ind_M<-as.numeric(m_star_vec>0)
  Y_mod <- glm(yi_vec ~ m_star_vec + x_i_vec + x_i_vec * m_star_vec, family = "gaussian")
  Y_mod_ini <- as.numeric(Y_mod$coefficients)
  ini_par[c(1, 2, 4, 6)] <- Y_mod_ini
  ini_par[c(3, 5)] <- 0
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





  # psi_est<-length(x_gp1)/(length(x_gp1)+length(x_gp2))

  # if(beta_gp1[1]<beta_gp2[1]) {
  #   beta_temp<-beta_gp1
  #   beta_gp1<-beta_gp2
  #   beta_gp2<-beta_temp
  #   psi_est<-1-psi_est
  # }
  # ini_par[7:8]<-beta_gp1
  # ini_par[9:10]<-beta_gp2

  # M_mod<-glm((1-ind_M)~x_i_vec,family = "binomial")
  Del_est <- c(-100, 0)
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


  lb_est <- c(-Inf, -Inf, 0, -Inf, 0, -Inf, -100, 0, 5, 0.5, rep(-10, k), rep(-10, k), rep(0.01, k - 1))
  ub_est <- c(Inf, Inf, 0, Inf, 0, Inf, -100, 0, Inf, Inf, rep(10, k), rep(10, k), rep(0.99, k - 1))


  return(list(ini_par, lb_est, ub_est))
}
