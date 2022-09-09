library(MarZIC)
library(SummarizedExperiment)
library(dirmult)

## A make up example with 1 taxon and 100 subjects.
set.seed(1)
nSub <- 200
nTaxa <- 10
## generate covariate of interest X
X <- rbinom(nSub, 1, 0.5)
## generate mean of each taxon. All taxon are having the same mean for simplicity.
mu <- exp(-5 + X) / (1 + exp(-5 + X))
phi <- 10

## generate true RA
M_taxon<-t(sapply(mu,function(x) rdirichlet(n=1,rep(x*phi,nTaxa))))

P_zero <- exp(-3 + 0.3 * X) / (1 + exp(-3 + 0.3 * X))

non_zero_ind <- t(sapply(P_zero,function(x) 1-rbinom(nTaxa,1,rep(x,nTaxa))))

True_RA<-t(apply(M_taxon*non_zero_ind,1,function(x) x/sum(x)))

## generate outcome Y based on true RA
Y <- 1 + 100 * True_RA[,1] + 5 * (True_RA[,1] > 0) + X + rnorm(nSub)

## library size was set to 10,000 for all subjects for simplicity.
libsize <- 10000

## generate observed RA
observed_AA <- floor(M_taxon*libsize*non_zero_ind)

observed_RA <- t(apply(observed_AA,1,function(x) x/sum(x)))
colnames(observed_RA)<-paste0("rawCount",seq_len(nTaxa))
## Construct SummerizedExperiment object
CovData <- cbind(Y = Y, X = X, libsize = libsize)
test_dat <-
  SummarizedExperiment(assays = list(MicrobData = t(observed_RA)), colData = CovData)

test_IFAA <- function() {
  res <- MarZIC(
    Experiment_dat = test_dat,
    lib_name = "libsize",
    y_name = "Y",
    x_name = "X",
    num_cores = 1,
    mediator_mix_range = 1
  )

  NIE1 <- res$NIE1_save
  sig_results<-subset(NIE1,significance==TRUE)
  checkEquals(rownames(sig_results), c("rawCount1"))

}
