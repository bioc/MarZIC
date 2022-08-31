expect_M_x <- function(psi, delta_1, mu) {
  res <- sum(psi * (1 - expit(delta_1)) * expit(mu))
  return(res)
}
