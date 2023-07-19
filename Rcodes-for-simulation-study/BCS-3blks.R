# load required packages
#library(INLA)
library(rstanarm)
rm(list = ls())
source("functions.r")
source("functions-bias.r")
source("functions-BCS.r")
source("simulation-scenarios.r")
N <- 100
critical_prob_combined <- 0.95
iters <- 1e3
alpha <- 0.2
#scenario <- scenario1
screenpool_max <- 275
for(scenario in list(scenario0,scenario1,scenario2,scenario3,scenario4,scenario5))
{
  # record trial performance
  cohort_var_enrichment <- matrix(NA,nrow=iters,ncol=3*length(scenario$`variable names`))
  colnames(cohort_var_enrichment) <- c(paste0("block1_mean_",scenario$`variable names`),
                                       paste0("block2_mean_",scenario$`variable names`),
                                       paste0("block3_mean_",scenario$`variable names`))
  opt_char <- matrix(NA,nrow = iters,ncol=12)
  colnames(opt_char) <- c("reject_H0","P_pE_ge_pC","P_pE_ge_pC_sens","P_pE_ge_pC_nonsens","Freq_pval",
                          "accrual_time","overtime",
                          "sensitivity","specificity",
                          "CPL","FPL","EFRR")
  est_TE <- matrix(NA,nrow=iters,ncol=2)
  colnames(est_TE) <- c("est_TE_sens","est_TE_nonsens")
  futility_prob_comb_IA1 <- futility_prob_D_IA1 <- futility_prob_comb_IA2 <- futility_prob_D_IA2 <- rep(NA,iters)
  ######################## start MC simulation  ##################################
  for(i in 1:iters)
  {
    set.seed(i*123)
    ##################### generate covariates for block 1#########################
    block1_covs <- covs_gen(N=N,dat_description = scenario)
    # also generate (potential) response both under control and treatment, no enrichment in block 1
    block1 <- response_gen(covs = block1_covs,dat_description = scenario,randomize = TRUE)
    
    mod1 <- stan_glm(response ~ treatment+c1+c2+
                       c1:treatment+c2:treatment, 
                     data = block1,
                     family = binomial(),
                     prior_intercept = normal(0,5),
                     prior=normal(location = c(0,0,0,0,0),
                                  scale=c(5,2,2,2,2)),
                     algorithm = "optimizing",
                     init="0",draws=5000)
    post_coefs1 <- t(as.matrix(update(mod1,iter=5000)))
    post_coefs1_mean <- mod1$coefficients
    
    # get the BCS
    BCS <- compute_BCS_RCS(dat = block1,post_coefs = post_coefs1,mod = mod1,alpha = alpha)
    # get the next block based on BCS
    # recruit D set, discard S^c set, and recruit with enrollment probability for S\D set
    futility_prob_comb_IA1[i] <- beta.ineq(block1)
    futility_prob_D_IA1[i] <- beta.ineq(block1[BCS$covs[,"In_S"],])
    ###################### get the second block of cohort#########################
    res <- get_nextblock_BCS(BCS,mod1,block1,N,scenario,screenpool_max = screenpool_max)
    block2 <- res[["next_block"]]
    block12 <- res[["accumulated_block"]]
    overtime <- res[["overtime"]]
    block2_screenpool <- res[["screenpool"]]
    
    # update of model and posterior draws of model parameters
    mod2 <- stan_glm(response ~ treatment+c1+c2+
                       c1:treatment+c2:treatment, 
                     data = block12,
                     family = binomial(),
                     prior=normal(location = c(0,0,0,0,0),
                                  scale=c(5,2,2,2,2)),
                     prior_intercept = normal(0,5),
                     algorithm = "optimizing",
                     init="0",draws=5000)
    post_coefs2 <- t(as.matrix(update(mod2,iter=5000)))
    post_coefs2_mean <- mod2$coefficients
    
    ###################### get the third block of cohort#########################
    BCS <- compute_BCS_RCS(dat = block12,post_coefs = post_coefs2,mod = mod2,alpha = alpha)
    futility_prob_comb_IA2[i] <- beta.ineq(block12)
    futility_prob_D_IA2[i] <- beta.ineq(block12[BCS$covs[,"In_S"],])
    res <- get_nextblock_BCS(BCS,mod2,block12,N,scenario,screenpool_max = screenpool_max)
    block3 <- res[["next_block"]]
    block123 <- res[["accumulated_block"]]
    overtime <- res[["overtime"]]
    block3_screenpool <- res[["screenpool"]]
    
    # update of model and posterior draws of model parameters
    mod3 <- stan_glm(response ~ treatment+c1+c2+
                       c1:treatment+c2:treatment, 
                     data = block123,
                     family = binomial(),
                     prior_intercept = normal(0,5),
                     prior=normal(location = c(0,0,0,0,0),
                                  scale=c(5,2,2,2,2)),
                     algorithm = "optimizing",
                     init="0",draws=5000)
    post_coefs3 <- t(as.matrix(update(mod3,iter=5000)))
    post_coefs3_mean <- mod3$coefficients
    
    ############################Bayesian decision rule#############################
    block123$sens <- pE_gt_pC(block123,post_coefs3_mean,0)
    P_pE_ge_pC <- beta.ineq(block123,0)
    P_pE_ge_pC_sens <- beta.ineq(block123[block123$sens==1,],0)
    P_pE_ge_pC_nonsens <- beta.ineq(block123[block123$sens==0,],0)
    Freq_pval <- compute_freq_pval2(list(block1,block2,block3))
    #  reject_H0 <-  ifelse(prop_lower==0.5,P_pE_ge_pC > critical_prob_combined,
    #                      P_pE_ge_pC_sens > critical_prob_sens)
    reject_H0 <- P_pE_ge_pC > critical_prob_combined
    
    ################################performance metrics ##############################
    ############################ estimated treatment effect ##########################
    est_TE[i,] <- compute_MCest_TE(rbind(subset(block1,select = scenario$`variable names`)
                                         ,subset(block2_screenpool,select = scenario$`variable names`),
                                         subset(block3_screenpool,select = scenario$`variable names`)),
                                   post_coefs3_mean,
                                   post_coefs3)
    
    
    ##################################accrual time ################################
    accrual_time <- max(block3_screenpool$patientID)
    
    ########################  current patient loss CPL#############################
    block23_screenpool <- rbind(block2_screenpool,block3_screenpool)
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
    if(all(is.na(block2)) & all(is.na(block3)))
    {
      sensitivity <- specificity <- 0
    } else
    {
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
    }
    
    ###### Expected posterior predictive accuracy for true sens patient ##########
    # EPPA <- compute_EPPA(scenario = scenario,post_coefs = post_coefs2)
    
    ########################### estimate treatment effect ########################
    # TE_sens <- compute_est_TE(block1,mod2)
    
    ####################### future patient metrics ##############################
    future_cohort_covs <- covs_gen(N=1e4,
                                   dat_description =scenario)
    # get the true response situation
    future_cohort <- response_gen(covs = future_cohort_covs,
                                  dat_description = scenario)
    
    future_cohort$sens_status <- pE_gt_pC(future_cohort,post_coefs3_mean)
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
    opt_char[i,] <- c(reject_H0,P_pE_ge_pC,P_pE_ge_pC_sens,P_pE_ge_pC_nonsens,Freq_pval,
                      accrual_time,overtime,
                      sensitivity,specificity,
                      CPL,FPL,EFRR)
    
    #if(i%%100==0) {print(i)}
    print(i)
  }
  df <- matrix(NA,nrow = 1e3,ncol = 4)
  colnames(df) <- c("f_prob_comb_IA1","f_prob_upper_IA1","f_prob_comb_IA2","f_prob_upper_IA2")
  df[,1] <- futility_prob_comb_IA1
  df[,2] <- futility_prob_D_IA1
  df[,3] <- futility_prob_comb_IA2
  df[,4] <- futility_prob_D_IA2
  save(opt_char,cohort_var_enrichment,est_TE,df,
       file=paste0("results_3blks/BCS-alpha",alpha,"-",scenario$`scenario name`,".RData"))
}
