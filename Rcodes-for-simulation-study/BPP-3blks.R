#------------------------------ load required packages
library(rstanarm)
rm(list = ls())
#------------------------------ load self-defined functions
source("functions.r")
source("functions-bias.r")
source("simulation-scenarios.r")
#------------------------------ general trial design parameters
N <- 100
critical_prob_combined <- 0.95
iters <- 1e3
screenpool_max <- 275
c <- 1.5
#------------------------------ record trial performance
for (scenario in list(scenario3,scenario4,scenario5))
{
  JSD <- matrix(NA,nrow=iters,ncol=2*length(scenario$`variable names`))
  colnames(JSD) <- c("block1-trt","block1-control","block12-trt","block12-control")
  cohort_var_enrichment <- matrix(NA,nrow=iters,ncol=3*length(scenario$`variable names`))
  colnames(cohort_var_enrichment) <- c(paste0("block1_mean_",scenario$`variable names`),
                                       paste0("block2_mean_",scenario$`variable names`),
                                       paste0("block3_mean_",scenario$`variable names`))
  opt_char <- matrix(NA,nrow = iters,ncol=14)
  colnames(opt_char) <- c("reject_H0","P_pE_ge_pC","P_pE_ge_pC_sens","P_pE_ge_pC_nonsens","Freq_pval",
                          "accrual_time","overtime",
                          "sensitivity","specificity",
                          "CPL","FPL","EFRR",
                          "prop_lower1","prop_lower2")
  est_TE <- matrix(NA,nrow=iters,ncol=2)
  colnames(est_TE) <- c("est_TE_sens","est_TE_nonsens")
  futility_df <- matrix(NA,nrow=iters,ncol=4)
  ######################## start MC simulation  ##################################
  for(i in 1:iters)
  {
    set.seed(i*321)
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
    
    
    # bootstrap estimate the quantiles of P(p_E > p_C) 
    block1$benefit_prob <- postprob_pE_gt_pC(covs=block1,post_coefs = post_coefs1)
    quantiles_est1 <- sapply(seq(0,1,by=0.01), function(q) quantile(block1$benefit_prob,q))
    block1$group <- block1$benefit_prob > quantiles_est1["50%"]
    JSD[i,1:2] <- compute_JSD(block1[block1$group==1,],block1[block1$group==0,])
    JSD_diff <- abs(JSD[i,1]-JSD[i,2])
    prop_lower1 <- min(0.5,exp(-JSD_diff)^c)
    futility_df[i,1:2] <- c(beta.ineq(block1),beta.ineq(block1[block1$group==1,],0))
    ###################### get the second block of cohort#########################
    res <- get_nextblock(type = "quantile",
                         post_coefs = post_coefs1,previous_blocks = block1,
                         N=N,scenario = scenario,
                         restriction = list(quantiles_est1,prop_lower1),
                         screenpool_max = screenpool_max)
    
    block2 <- res[["next_block"]]
    block12 <- res[["accumulated_block"]]
    overtime <- res[["overtime"]]
    block2_screenpool <- res[["screenpool"]]
    
    # update of model and posterior draws of model parameters
    mod2 <- stan_glm(response ~ treatment+c1+c2+
                       c1:treatment+c2:treatment, 
                     data = block12,
                     family = binomial(),
                     prior_intercept = normal(0,5),
                     prior=normal(location = c(0,0,0,0,0),
                                  scale=c(5,2,2,2,2)),
                     algorithm = "optimizing",
                     init="0",draws=5000)
    post_coefs2 <- t(as.matrix(update(mod2,iter=5000)))
    post_coefs2_mean <- mod2$coefficients
    
    # bootstrap estimate the quantiles of P(p_E > p_C) 
    
    ###################### get the third block of cohort#########################
    # bootstrap re-estimate the quantiles of P(p_E > p_C) 
    block12$benefit_prob <- postprob_pE_gt_pC(covs=block12,post_coefs = post_coefs2)
    block12_screenpool <- rbind(subset(block1,select = scenario$`variable names`)
                                ,subset(block2_screenpool,select = scenario$`variable names`))
    block12_screenpool$benefit_prob <- postprob_pE_gt_pC(covs=block12_screenpool,post_coefs = post_coefs2)
    quantiles_est2 <- sapply(seq(0,1,by=0.01), function(q) quantile(block12_screenpool$benefit_prob,q))
    block12$group <- block12$benefit_prob > quantiles_est2["50%"]
    JSD[i,3:4] <- compute_JSD(block12[block12$group==1,],block12[block12$group==0,])
    JSD_diff <- abs(JSD[i,3]-JSD[i,4])
    prop_lower2 <- min(0.5,exp(-JSD_diff)^c)
    futility_df[i,3:4] <- c(beta.ineq(block12), beta.ineq(block12[block12$group==1,],0))
    
    res <- get_nextblock(type = "quantile",
                         post_coefs = post_coefs2,previous_blocks = block12,
                         N=N,scenario = scenario,
                         restriction = list(quantiles_est2,prop_lower2),
                         screenpool_max = screenpool_max)
    
    block3 <- res[["next_block"]]
    block123 <- res[["accumulated_block"]]
    overtime <- res[["overtime"]]
    block3_screenpool <- res[["screenpool"]]
    
    #################### update model #############################################
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
                      CPL,FPL,EFRR,
                      prop_lower1,prop_lower2)
    #if(i%%100==0) {print(i)}
    print(i)
  }
  save(opt_char,JSD,cohort_var_enrichment,est_TE,futility_df,
       file=paste0("results_3blks/BPP-c",c,"-",scenario$`scenario name`,"new.RData"))
}