# load required packages
library(rstanarm)
rm(list = ls())
source("functions.r")
source("functions-bias.r")
source("simulation-scenarios.r")
# general design pamameters
N <- c(168,144)
critical_prob_combined <- 0.95
iters <- 1e3
scenario <- scenario2
#eta <- 0.3
# record trial performance
for(eta in c(0.3,0.5,0.7))
{
cohort_var_enrichment <- matrix(NA,nrow=iters,ncol=2*length(scenario$`variable names`))
colnames(cohort_var_enrichment) <- c(paste0("block1_mean_",scenario$`variable names`),
                                     paste0("block2_mean_",scenario$`variable names`))
opt_char <- matrix(NA,nrow = iters,ncol=12)
colnames(opt_char) <- c("reject_H0","P_pE_ge_pC","P_pE_ge_pC_sens","P_pE_ge_pC_nonsens","Freq_pval",
                        "accrual_time","overtime",
                        "sensitivity","specificity",
                        "CPL","FPL","EFRR")
est_TE <- matrix(NA,nrow=iters,ncol=2)
colnames(est_TE) <- c("est_TE_sens","est_TE_combined")
######################## start MC simulation  ##################################
for(i in 1:iters)
{
  set.seed(i*123)
  ##################### generate covariates for block 1#########################
  block1_covs <- covs_gen(N=N[1],dat_description = scenario)
  # colMeans(block1_covs)
  # also generate (potential) response both under control and treatment, no enrichment in block 1
  block1 <- response_gen(covs = block1_covs,dat_description = scenario,randomize = TRUE)
  # mean(block1$response[block1$treatment==1])
  mod1 <- stan_glm(response ~ treatment+age+sex+ECOG+PD_L1+TMB+
                     age:treatment+sex:treatment+ECOG:treatment+PD_L1:treatment+TMB:treatment, 
                   data = block1,
                   family = binomial(),
                   prior_intercept = normal(0,5),
                   prior=normal(location = c(0,0,0,0,0,0,0,0,0,0,0),
                                scale=c(5,2,2,2,2,2,2,2,2,2,2)),
                   algorithm = "optimizing",
                   init="0",draws=5000)
  post_coefs1 <- t(as.matrix(update(mod1,iter=5000)))
  post_coefs1_mean <- mod1$coefficients
  
  
  # prop_lower <- 0
  ###################### get the second block of cohort#########################
  res <- get_nextblock(type = "fixed",
                       post_coefs = post_coefs1,previous_blocks = block1,
                       N=N[2],scenario = scenario,
                       restriction = list(eta))
  
  block2 <- res[["next_block"]]
  block12 <- res[["accumulated_block"]]
  overtime <- res[["overtime"]]
  block2_screenpool <- res[["screenpool"]]
  # beta.ineq(block12)
  
  # update of model and posterior draws of model parameters
  mod2 <- stan_glm(response ~ treatment+age+sex+ECOG+PD_L1+TMB+
                     age:treatment+sex:treatment+ECOG:treatment+PD_L1:treatment+TMB:treatment, 
                   data = block12,
                   family = binomial(),
                   prior_intercept = normal(0,5),
                   prior=normal(location = c(0,0,0,0,0,0,0,0,0,0,0),
                                scale=c(5,2,2,2,2,2,2,2,2,2,2)),
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
                                 colMeans(block2[2:(1+length(scenario$`variable names`))]))
  opt_char[i,] <- c(reject_H0,P_pE_ge_pC,P_pE_ge_pC_sens,P_pE_ge_pC_nonsens,Freq_pval,
                    accrual_time,overtime,
                    sensitivity,specificity,
                    CPL,FPL,EFRR)
  #if(i%%100==0) {print(i)}
  print(i)
}
save(opt_char,cohort_var_enrichment,est_TE,
     file=paste0("results_new/FC-eta0",eta*10,"-",scenario$`scenario name`,".RData"))
}