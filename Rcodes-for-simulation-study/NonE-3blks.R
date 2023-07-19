library(rstanarm)
rm(list = ls())
source("functions.r")
source("functions-bias.r")
source("simulation-scenarios.r")
N <- 100
iters <- 1e3
critical_prob_combined <- 0.95
for (scenario in list(scenario0,scenario1,scenario2,scenario3,scenario4))
{
  # two continuous variables from a bivariate normal distribution
  ########################  Non enrichment design ######################################
  
  cohort_var_enrichment <- matrix(NA,nrow=iters,ncol=3*length(scenario$`variable names`))
  colnames(cohort_var_enrichment) <- c(paste0("block1_mean_",scenario$`variable names`),
                                      paste0("block2_mean_",scenario$`variable names`),
                                      paste0("block3_mean_",scenario$`variable names`))
  opt_char <- matrix(NA,nrow = iters,ncol=11)
  colnames(opt_char) <- c("reject_H0","P_pE_ge_pC","P_pE_ge_pC_sens","P_pE_ge_pC_nonsens","Freq_pval",
                          "accrual_time",
                          "sensitivity","specificity",
                          "CPL","FPL","EFRR")
  est_TE <- matrix(NA,nrow=iters,ncol=2)
  colnames(est_TE) <- c("est_TE_sens","est_TE_combined")
  for(i in 1:iters)
  {
    set.seed(i*321)
    # get all patients (two blocks seamlessly without enrichment)
    block1 <- covs_gen(N=N,dat_description = scenario,pID = 1:N)
    # randomize and get the response
    block1 <- response_gen(covs = block1,dat_description = scenario,randomize = TRUE)
    
    block2 <- covs_gen(N=N,dat_description = scenario,pID = (1+N):(2*N))
    # randomize and get the response
    block2 <- response_gen(covs = block2,dat_description = scenario,randomize = TRUE)
    
    block3 <- covs_gen(N=N,dat_description = scenario,pID = (1+2*N):(3*N))
    # randomize and get the response
    block3 <- response_gen(covs = block3,dat_description = scenario,randomize = TRUE)
    
    
    block123 <- rbind(block1,block2,block3)
    block23_screenpool <- rbind(block2,block3)
    
    mod <- stan_glm(response ~ treatment+c1+c2+
                      c1:treatment+c2:treatment, 
                    data = block123,
                    family = binomial(),
                    prior_intercept = normal(0,5),
                    prior=normal(location = c(0,0,0,0,0),
                                 scale=c(5,2,2,2,2)),
                    algorithm = "optimizing",
                    init="0",draws=5000)
    post_coefs <- t(as.matrix(update(mod,iter=5000)))
    post_coefs_mean <- mod$coefficients
    
    ############################Bayesian decision rule#############################
    block123$sens <- pE_gt_pC(block123,post_coefs_mean)
    P_pE_ge_pC <- beta.ineq(block123,0)
    P_pE_ge_pC_sens <- beta.ineq(block123[block123$sens==1,],0)
    P_pE_ge_pC_nonsens <- beta.ineq(block123[block123$sens==0,],0)
    Freq_pval <- compute_freq_pval2(list(block1,block2,block3))
    reject_H0 <-  P_pE_ge_pC > critical_prob_combined
    
    ################################performance metrics ##############################
    ############################ estimated treatment effect ##########################
    est_TE[i,] <- compute_MCest_TE(subset(block123,select = scenario$`variable names`),
                                   post_coefs_mean,
                                   post_coefs)
    
    ##################################accrual time ################################
    accrual_time <- 3*N
    
    ########################  current patient loss CPL#############################
    block23_screenpool$loss <- sapply(1:nrow(block23_screenpool), function(x)
    {
      individual <- block23_screenpool[x,]
      if ((individual$pE_true >= individual$pC_true & individual$treatment==1) | (individual$pE_true < individual$pC_true & individual$treatment==0))
      {
        res <- 0
      } else
      {
        res <- abs(individual$pE_true-individual$pC_true)
      }
      return(res)
    })
    CPL <- mean(block23_screenpool$loss)
    
    #######################  sensitivity and specificity #######################
    block23 <- rbind(block2,block3)
    # number of truly sensitive patient encountered
    true_sens_num <- sum(block23_screenpool$pE_true>block23_screenpool$pC_true)
    # number of truly sensitive patient enrolled
    enroll_sens_num <- sum(block23$pE_true > block23$pC_true)
    sensitivity <- enroll_sens_num/true_sens_num
    
    # number of truly non-sensitive patient encountered
    true_nonsens_num <- sum(block23_screenpool$pE_true<=block23_screenpool$pC_true)
    # number of truly non-sensitive patient enrolled
    enroll_nonsens_num <- sum(block23$pE_true <= block23$pC_true)
    specificity <- (true_nonsens_num-enroll_nonsens_num)/true_nonsens_num
    
    ###### Expected posterior predictive accuracy for true sens patient ##########
    # EPPA <- compute_EPPA(scenario = scenario,post_coefs = post_coefs)
    
    ########################### estimate treatment effect ########################
    #TE_sens <- compute_est_TE(block12,mod)
    # TE_whole <- compute_est_TE(block12,mod)
    
    ####################### future patient metrics ##############################
    future_cohort_covs <- covs_gen(N=1e4,
                                   dat_description =scenario)
    # get the true response situation
    future_cohort <- response_gen(covs = future_cohort_covs,
                                  dat_description = scenario)
    
    future_cohort$sens_status <- pE_gt_pC(future_cohort,post_coefs_mean)
    # if declare efficacy, will treat future patient differently according to their sensitive
    # status, other wise we stick to SOC
    
    future_cohort$enrolled <- reject_H0*(future_cohort$sens_status)
    
    
    ################### expected future response rate ############################
    # patients enrolled will receive new treatment otherwise SOC
    future_cohort$response_E <- rbinom(nrow(future_cohort),1,future_cohort$pE_true)
    future_cohort$response_C <- rbinom(nrow(future_cohort),1,future_cohort$pC_true)
    future_cohort$response <- sapply(1:nrow(future_cohort), function(i)
    {
      ifelse(future_cohort[i,]$enrolled==1,future_cohort[i,]$response_E,
             future_cohort[i,]$response_C)
    })
    EFRR <- mean(future_cohort$response)
    
    #################################  FPL #####################################
    future_cohort$loss <- sapply(1:nrow(future_cohort), function(x)
    {
      individual <- future_cohort[x,]
      if ((individual$pE_true >= individual$pC_true & individual$enrolled==1) | (individual$pE_true < individual$pC_true & individual$enrolled==0))
      {
        res <- 0
      } else
      {
        res <- abs(individual$pE_true-individual$pC_true)
      }
      return(res)
    })
    FPL <- mean(future_cohort$loss)
    
    # record operating characteristics
    cohort_var_enrichment[i,] <- c(colMeans(block1[2:(1+length(scenario$`variable names`))]),
                                   colMeans(block2[2:(1+length(scenario$`variable names`))]),
                                   colMeans(block3[2:(1+length(scenario$`variable names`))]))
    
    # record operating characteristics
    opt_char[i,] <- c(reject_H0,P_pE_ge_pC,P_pE_ge_pC_sens,P_pE_ge_pC_nonsens,Freq_pval,
                      accrual_time,
                      sensitivity,specificity,
                      CPL,FPL,EFRR)
    print(i)
  }
  save(opt_char,est_TE,cohort_var_enrichment,
       file=paste0("results_3blks/NonE-",scenario$`scenario name`,".RData"))
}

