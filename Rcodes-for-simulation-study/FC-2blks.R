library(rstanarm)
rm(list = ls())
source("functions.r")
source("functions-bias.r")
source("simulation-scenarios.r")
# general design pamameters
#scenario <- scenario1
N <- 150
critical_prob_combined <- 0.95
iters <- 1e3
screenpool_max <- 500
scenario <- scenario5
for(eta in c(0.3,0.5,0.7))
{
  est_TE <- matrix(NA,nrow=iters,ncol=2)
  colnames(est_TE) <- c("est_TE_sens","est_TE_nonsens")
  
  cohort_var_enrichment <- matrix(NA,nrow=iters,ncol=2*length(scenario$`variable names`))
  colnames(cohort_var_enrichment) <- c(paste0("block1_mean_",scenario$`variable names`),
                                       paste0("block2_mean_",scenario$`variable names`))
  # at the end of trial, use a Bayesian rule to decide whether to declare the new intervention a success
  # if success, the future clinical practice would be treat sensitive patient with new treatment and non-sensitive with SOC
  # if not declared a success, all future patients will still use SOC (i.e, the new treatment would not be approves)
  opt_char <- matrix(NA,nrow = iters,ncol=12)
  colnames(opt_char) <- c("reject_H0","P_pE_ge_pC","P_pE_ge_pC_sens","P_pE_ge_pC_nonsens","Freq_pval",
                          "accrual_time","overtime",
                          "sensitivity","specificity",
                          "CPL","FPL","EFRR")
  for(i in 1:iters)
  {
    set.seed(i*123)
    ##################### generate covariates for block 1#########################
    block1_covs <- covs_gen(N=N,dat_description = scenario)
    # also generate (potential) response both under control and treatment, no enrichment in block 1
    block1 <- response_gen(covs = block1_covs,dat_description = scenario,randomize = TRUE)
    
    
    # model the conditional response probability via Bayesian logistic regression
    # model based on first block of data
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
    
    res <- get_nextblock(type = "fixed",
                         post_coefs = post_coefs1,previous_blocks = block1,
                         N=N,scenario = scenario,
                         restriction = list(eta),
                         screenpool_max = screenpool_max)
    block2 <- res[["next_block"]]
    block12 <- res[["accumulated_block"]]
    overtime <- res[["overtime"]]
    block2_screenpool <- res[["screenpool"]]
    
    # update of model and posterior draws of model parameters
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
    
    ############################Bayesian decision rule#############################
    block12$sens <- pE_gt_pC(block12,post_coefs2_mean,0)
    P_pE_ge_pC <- beta.ineq(block12,0)
    P_pE_ge_pC_sens <- beta.ineq(block12[block12$sens==1,],0)
    P_pE_ge_pC_nonsens <- beta.ineq(block12[block12$sens==0,],0)
    Freq_pval <- compute_freq_pval(block1,block2)
    reject_H0 <-  P_pE_ge_pC > critical_prob_combined
    
    ################################performance metrics ##############################
    ############################ estimated treatment effect ##########################
    est_TE[i,] <- compute_MCest_TE(rbind(subset(block1,select = scenario$`variable names`)
                                         ,subset(block2_screenpool,select = scenario$`variable names`)),
                                   post_coefs2_mean,
                                   post_coefs2)
    
    ##################################accrual time ################################
    accrual_time <- max(block2_screenpool$patientID)
    
    ########################  current patient loss CPL#############################
    block2_screenpool$loss <- sapply(1:nrow(block2_screenpool), function(x)
    {
      individual <- block2_screenpool[x,]
      if ((individual$pE_true >= individual$pC_true & individual$treatment==1) | (individual$pE_true < individual$pC_true & individual$treatment==0))
      {
        res <- 0
      } else
      {
        res <- abs(individual$pE_true-individual$pC_true)
      }
      return(res)
    })
    CPL <- mean(block2_screenpool$loss)
    
    #######################  sensitivity and specificity #######################
    if(all(is.na(block2)))
    {
      sensitivity <- specificity <- 0
    } else
    {
      # number of truly sensitive patient encountered
      true_sens_num <- sum(block2_screenpool$pE_true>block2_screenpool$pC_true)
      # number of truly sensitive patient enrolled
      enroll_sens_num <- sum(block2$pE_true > block2$pC_true)
      sensitivity <- enroll_sens_num/true_sens_num
      
      # number of truly non-sensitive patient encountered
      true_nonsens_num <- sum(block2_screenpool$pE_true<=block2_screenpool$pC_true)
      # number of truly non-sensitive patient enrolled
      enroll_nonsens_num <- sum(block2$pE_true <= block2$pC_true)
      specificity <- (true_nonsens_num-enroll_nonsens_num)/true_nonsens_num
    }
    
    
    ####################### future patient metrics ##############################
    future_cohort_covs <- covs_gen(N=1e4,
                                   dat_description =scenario)
    # get the true response situation
    future_cohort <- response_gen(covs = future_cohort_covs,
                                  dat_description = scenario)
    
    future_cohort$sens_status <- pE_gt_pC(future_cohort,post_coefs2_mean)
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
                                   colMeans(block2[2:(1+length(scenario$`variable names`))],
                                            na.rm = T))
    
    
    # record operating characteristics
    opt_char[i,] <- c(reject_H0,P_pE_ge_pC,P_pE_ge_pC_sens,P_pE_ge_pC_nonsens,Freq_pval,
                      accrual_time,overtime,
                      sensitivity,specificity,
                      CPL,FPL,EFRR)
    print(i)
  }
  save(opt_char,cohort_var_enrichment,est_TE,
      file=paste0("results_2blks/FC-eta0",eta*10,"-",scenario$`scenario name`,".RData"))
}