compute_MCest_TE <- function(dat,post_coefs_mean,post_coefs)
{
    M <- nrow(dat)
    post_coefs_mean <- matrix(post_coefs_mean,ncol=1)
    covs <- as.matrix(dat,nrow = M)
    covs_treat <- matrix(c(rep(1,M),rep(1,M),
                           covs,covs),byrow = FALSE,nrow = M)
    pE_MC <- ilogit(covs_treat%*%post_coefs_mean)
    covs_control <- matrix(c(rep(1,M),rep(0,M),
                             covs,rep(rep(0,M),ncol(covs))),byrow = FALSE,nrow = M)
    pC_MC <- ilogit(covs_control%*%post_coefs_mean)
    dat$post_mean_te <- pE_MC - pC_MC
    dat$post_prob_benefit <- postprob_pE_gt_pC(covs,post_coefs)
    est_TE_sens <- mean(dat$post_mean_te*dat$post_prob_benefit)/mean(dat$post_prob_benefit)
    est_TE_combined <- mean(dat$post_mean_te)
  return(c(est_TE_sens,est_TE_combined))
} 




compute_MCtrue_TE <- function(scenario,M=1e5)
{
  set.seed(86)
  require(mvtnorm)
  covs <- covs_gen(N=M,dat_description = scenario)
  covs2 <- response_gen(covs = covs,dat_description = scenario)
  TEs <- covs2$pE_true - covs2$pC_true
  sens_size <- mean(TEs>0)
  TE_sens <- ifelse(sens_size==0,0,mean(TEs[TEs>0]))
  TE_nonsens <- ifelse(sens_size==1,0,mean(TEs[TEs<=0]))
  TE_combined <- mean(TEs)
  return(c("true_TE_sens"=TE_sens,
           "true_TE_nonsens"=TE_nonsens,
           "true_TE_combined"=TE_combined,
           "true_sens_size"=sens_size))
}


