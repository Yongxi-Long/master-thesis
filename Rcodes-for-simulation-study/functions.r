###################################
# function to generate covariates #
###################################
#-------------------------------------------------------------------------------
# N: total number of patients to generate
# TE: the main treatment effect
# var_types: a vector of covariate type (continuous or binary)
# var_names: a vector of covariate names
# var_means: a vector of covariate means (continuous) or proportions (binary)
# var_sds: a vector of standard deviations for continuous variables,sd=1 for binary
# betas: main prognostic effects of covariates
# gammas: predictive effects (interaction with treatment) of covariates
#-------------------------------------------------------------------------------
covs_gen <- function(N,dat_description,pID=1:N)
{
  require(MASS)
  var_names <- dat_description$`variable names`
  var_means <- dat_description$`variables means`
  var_sds <- dat_description$`variable standard deviations`
  rho <- dat_description$`variable correlation`
  binary <- dat_description$`binary`
  cov_num <- length(var_names)
  # whether it is a binary variable
  # covariance matrix
  Vmat <- matrix(0, cov_num, cov_num)
  for(i in 1:ncol(Vmat)){
    for(j in 1:nrow(Vmat)){
      if(i==j){Vmat[i, j] <- var_sds[i]^2}
      if(i!=j){Vmat[i, j] <- rho*var_sds[i]*var_sds[j]}
    }
  }
  # generate from MVN
  covs <- mvrnorm(n=N,mu=var_means,Sigma = Vmat)
  for(i in 1:cov_num)
  {
    if(binary[i])
    {
    covs[,i] <- 1*(covs[,i]) >= var_means[i]
    }
  }
  covs <- data.frame(cbind(pID,covs))
  colnames(covs) <- c("patientID",var_names)
  return(covs)
}
###########################################
#    function to generate response        #
###########################################
response_gen <- function(b0=-0.5,covs,dat_description,randomize=FALSE)
{
  N <- nrow(covs)
  TE <- dat_description$`treatment effect`
  var_names <- dat_description$`variable names`
  betas <- as.matrix(dat_description$`prognostic effects`,nrow=2)
  gammas <- as.matrix(dat_description$`predictive effects`,nrow=2)
  covs2 <- as.matrix(subset(covs,select=var_names),nrow=N,ncol = length(var_names))
  
  
  # potential outcomes under treatment
  p_Es <- ilogit(b0+TE+covs2%*%(betas+gammas))
  y_E <- rbinom(N,1,p_Es)
  
  # potential outcome under control
  p_Cs <- ilogit(b0+covs2%*%betas)
  y_C <- rbinom(N,1,p_Cs)
  
  if(randomize)
  {
    # randomize and get the observed outcome (1/2 to control, 1/2 to treatment)
    index <- c(rep(0,round(N/2)),rep(1,round(N/2)))
    treatment <- sample(index)
    y <- sapply(1:N, function(n)
    {
      ifelse(treatment[n],y_E[n],y_C[n])
    })
    # construct the data frame
    df <- data.frame(cbind(covs,matrix(p_Es,ncol=1),matrix(p_Cs,ncol=1),treatment,y))
    colnames(df) <- c(colnames(covs),"pE_true","pC_true","treatment","response")
    return(df)
  }
  else{
    # if randomization is not needed, only return the true response probability under treatment/control
    df <- data.frame(cbind(covs,matrix(p_Es,ncol=1),matrix(p_Cs,ncol=1)))
    colnames(df) <- c(colnames(covs),"pE_true","pC_true")
    return(df)
  }
  
}


##########################################
# function to compute inverse of a logit #
##########################################
ilogit <- function(x)
{
  if(any(x>500)) {return(1)}
  else{
    exp(x)/(1+exp(x))
  }
}

###########################################
#  function to compute P(p_E(x)>p_c(x))   #
###########################################
postprob_pE_gt_pC <- function(covs,post_coefs)
{
  
  # get the classification results-----------------------------------------------
  new_patients <- subset(covs,select = scenario$`variable names`)
  
  model_matrix_E <- as.matrix(cbind(rep(1,nrow(new_patients)),
                                    rep(1,nrow(new_patients)),
                                    new_patients,
                                    new_patients))
  post_pEs <- ilogit(model_matrix_E %*% post_coefs)
  model_matrix_C <- as.matrix(cbind(rep(1,nrow(new_patients)),
                                    rep(0,nrow(new_patients)),
                                    new_patients,
                                    matrix(rep(rep(0,nrow(new_patients)),
                                               length(scenario$`variable names`)),
                                           ncol =length(scenario$`variable names`))
  ))
  post_pCs <- ilogit(model_matrix_C %*% post_coefs)
  
  benefit_prop <- rowMeans(post_pEs > post_pCs)
  return(benefit_prop)
}

##########################################################
# function to calculate whether pe > pc                  #
##########################################################
# margin: minimum meaningful clinical treatment effect
pE_gt_pC <- function(covs,post_coefs_mean,margin=0)
{
  require(dplyr)
  post_coefs_mean <- as.matrix(post_coefs_mean,ncol=1)
  new_patients <- subset(covs,select = scenario$`variable names`)
  model_matrix_E <- as.matrix(cbind(rep(1,nrow(new_patients)),
                                    rep(1,nrow(new_patients)),
                                    new_patients,
                                    new_patients))
  pE <- ilogit(model_matrix_E %*% post_coefs_mean)
  model_matrix_C <- as.matrix(cbind(rep(1,nrow(new_patients)),
                                    rep(0,nrow(new_patients)),
                                    new_patients,
                                    matrix(rep(rep(0,nrow(new_patients)),
                                               length(scenario$`variable names`)),
                                           ncol =length(scenario$`variable names`))
                                    ))
  pC <- ilogit(model_matrix_C %*% post_coefs_mean)
  res <- 1*(pE>(pC+margin))
  return(res)
}

#####################################################
# function to simulate from estimated covariate distribution
########################################################
simulate_covs <- function(n_sample,dat)
{
  mu <- c(mean(dat$c1),mean(dat$c2))
  sd <- c(sd(dat$c1),sd(dat$c2))
  rho <- cor(dat$c1,dat$c2)
  Vmat <- matrix(NA,2,2)
  for(i in 1:ncol(Vmat)){
    for(j in 1:nrow(Vmat)){
      if(i==j){Vmat[i, j] <- sd[i]^2}
      if(i!=j){Vmat[i, j] <- rho*sd[i]*sd[j]}
    }
  }
  covs <- mvrnorm(n=n_sample,mu,Vmat)
  colnames(covs) <- c("c1","c2")
  return(covs)
}

#########################################################
# Bayesian power P(p_E>p_C), default using uniform prior
#########################################################
beta.ineq <- function(dat, delta=0)
{ 
  # calculate parameters for Pe
  a <- 1+sum(dat$response[dat$treatment==1])
  b <- 1+sum(1-dat$response[dat$treatment==1])
  c <- 1+sum(dat$response[dat$treatment==0])
  d <- 1+sum(1-dat$response[dat$treatment==0])
  if (a <=0 | b <= 0 | c <=0 | d <= 0) 
    stop("paramters has to be positive")
  if (a <0.01 | b < 0.01 | c <0.01 | d < 0.01) 
    stop("paramters are to close to 0")
  if (delta>1) 
    stop("delta>1!")
  if (delta<0) 
    stop("delta<0!")
  integrand <- function(x) { dbeta(x, a, b)*pbeta(x-delta, c, d) }
  integrate(integrand, delta, 1, rel.tol=1e-4)$value
}

#######################################################
# function to compute frequentist power               #
#######################################################
# use the inverse normal combination test statistics
# test the strong null: p_E(x) <= p_C(x) for all covariates x
# alpha: significance level
# return 0/1: reject H0 or not
compute_freq_pval <- function(block1,block2)
{
  if(all(is.na(block2)))
  {
    # vector of patient number under treatment/control for each block
    n_k <- nrow(block1)/2
    # vector of response prob in the treated population per block
    pE_k <- mean(block1$response[block1$treatment==1])
    # vector of response prob in the control population per block
    pC_k <- mean(block1$response[block1$treatment==0])
    # vector of pooled response rate per block
    pPool_k <- (pE_k+pC_k)/2
    test_stat <- (pE_k-pC_k)/(2*sqrt(pPool_k*(1-pPool_k)/n_k))
  } else {
    # vector of patient number under treatment/control for each block
    n_ks <- c(nrow(block1)/2,nrow(block2)/2)
    n <- sum(n_ks)
    # vector of response prob in the treated population per block
    pE_ks <- c(mean(block1$response[block1$treatment==1]),
               mean(block2$response[block2$treatment==1]))
    # vector of response prob in the control population per block
    pC_ks <- c(mean(block1$response[block1$treatment==0]),
               mean(block2$response[block2$treatment==0]))
    # vector of pooled response rate per block
    pPool_ks <- (pE_ks+pC_ks)/2
    inner <- (pE_ks-pC_ks)/(sqrt(2*pPool_ks*(1-pPool_ks)/n_ks)) 
    inner[is.na(inner)] <- 0
    test_stat <- 1/sqrt(n/2)*sum(sqrt(n_ks/2)*( inner))
  }
  # one-sided,reject for large value of test stat under N(0,1)
  p_val <- pnorm(q=test_stat,lower.tail = FALSE)
  return(p_val)
}

compute_freq_pval2 <- function(block_list)
{
  # remove empty blocks
  block_list <- block_list[!is.na(block_list)]
  n_ks <- sapply(block_list, function(x) nrow(x)/2)
  n <- sum(n_ks)
  # vector of response prob in the treated population per block
  pE_ks <- sapply(block_list, function(x) mean(x$response[x$treatment==1]))
  # vector of response prob in the control population per block
  pC_ks <- sapply(block_list, function(x) mean(x$response[x$treatment==0]))
  # vector of pooled response rate per block
  pPool_ks <- (pE_ks+pC_ks)/2
  inner <- (pE_ks-pC_ks)/(sqrt(2*pPool_ks*(1-pPool_ks)/n_ks)) 
  inner[is.na(inner)] <- 0
  test_stat <- 1/sqrt(n/2)*sum(sqrt(n_ks/2)*( inner))
  # one-sided,reject for large value of test stat under N(0,1)
  p_val <- pnorm(q=test_stat,lower.tail = FALSE)
  return(p_val)
}

##################################################################
# function to get next block
# function to get next block
get_nextblock <- function(type,
                          post_coefs,previous_blocks,
                          N,scenario,
                          restriction=NA,screenpool_max=500)
{
 # set.seed(i*86)
  pID_start <- max(previous_blocks$patientID)+1
  block2_covs <- covs_gen(N=4e3,dat_description = scenario,
                          pID=101:4100)
  block2_covs <- block2_covs[block2_covs$patientID>=pID_start,]
  count <- 1
  enrolled <- 0
  overtime <- 0
  block2_enrolled <- data.frame(matrix(NA,nrow=N,ncol=(1+length(scenario$`variable names`))))
  colnames(block2_enrolled) <- c("patientID",scenario$`variable names`)
  if(type=="quantile")
  {
    quantiles_est <- restriction[[1]]
    prop_lower <- restriction[[2]]
    repeat{
      patient <- block2_covs[count,]
      patient$benefit_prob <- postprob_pE_gt_pC(covs=patient,post_coefs = post_coefs)
      if(patient$benefit_prob >  quantiles_est["50%"])
      {
        enrolled <- enrolled+1
        block2_enrolled[enrolled,] <- subset(patient,select=c("patientID",scenario$`variable names`))
        count <- count+1
      } else # uncertainty
      {
        # these are the patients in the lower half quantile
        # scale the quantile of current patient's benefit prob from (0,0.5) to (0,1)
        patient$quantile <-  (which.min(abs(quantiles_est-patient$benefit_prob))-1)/100*2
        # under sample the lower quantile such that they make up prop_lower % of the future block
        # so on average their enroll prob is prop_lower/(1-prop_lower)
        # if prop_lower/(1-prop_lower <= 1/2, can enroll proportional to the quantile scaled by prop_lower/(1-prop_lower)/0.5
        p_enroll_lower_avg <- prop_lower/(1-prop_lower)
        if(p_enroll_lower_avg <= 0.5)
        {
          enroll_prob <- patient$quantile*p_enroll_lower_avg*2
        } else
        {
          # have to tun upper X% to 1 to make sure that prop_enroll on average is achieved
          cutpoint <- (2-sqrt(4-8*(1-p_enroll_lower_avg)))/2
          enroll_prob <- ifelse(patient$quantile>=cutpoint,1,patient$quantile)
        }
        enroll <- rbinom(1,1,prob = enroll_prob)
        if (enroll)
        {
          enrolled <- enrolled+1
          block2_enrolled[enrolled,] <- subset(patient,select=c("patientID",scenario$`variable names`))
          count <- count+1
        } else
        {
          count <- count + 1
        }
      }
      if(count <=(1+screenpool_max) & enrolled==N){
        break
      }
      if(count==(1+screenpool_max) & enrolled < N)
      {
        overtime <- 1
        break
        
      }
    }
  } else if(type=="quantile_eq")
  {
    quantiles_est <- restriction[[1]]
    prop_lower <- restriction[[2]]
    repeat{
      patient <- block2_covs[count,]
      patient$benefit_prob <- postprob_pE_gt_pC(covs=patient,post_coefs = post_coefs)
      if(patient$benefit_prob >  quantiles_est["50%"])
      {
        enrolled <- enrolled+1
        block2_enrolled[enrolled,] <- subset(patient,select=c("patientID",scenario$`variable names`))
        count <- count+1
      } else # uncertainty
      {
        # these are the patients in the lower half quantile
        # scale the quantile of current patient's benefit prob from (0,0.5) to (0,1)
        patient$quantile <-  (which.min(abs(quantiles_est-patient$benefit_prob))-1)/100*2
        # under sample the lower quantile such that they make up prop_lower % of the future block
        # so on average their enroll prob is prop_lower/(1-prop_lower)
        # if prop_lower/(1-prop_lower <= 1/2, can enroll proportional to the quantile scaled by prop_lower/(1-prop_lower)/0.5
        p_enroll_lower_avg <- prop_lower/(1-prop_lower)
        enroll_prob <- p_enroll_lower_avg
        enroll <- rbinom(1,1,prob = enroll_prob)
        if (enroll)
        {
          enrolled <- enrolled+1
          block2_enrolled[enrolled,] <- subset(patient,select=c("patientID",scenario$`variable names`))
          count <- count+1
        } else
        {
          count <- count + 1
        }
      }
      if(count <=(1+screenpool_max) & enrolled==N){
        break
      }
      if(count==(1+screenpool_max) & enrolled < N)
      {
        overtime <- 1
        break
        
      }
    }
  } else if (type=="fixed")
  {
    fixed_cutpoint <- restriction[[1]]
    repeat{
      patient <- block2_covs[count,]
      patient$benefit_prob <- postprob_pE_gt_pC(covs=patient,post_coefs = post_coefs)
      if(patient$benefit_prob > fixed_cutpoint)
      {
        enrolled <- enrolled+1
        block2_enrolled[enrolled,] <- subset(patient,select=c("patientID",scenario$`variable names`))
        count <- count+1
      } else # discard
      {
        count <- count + 1
      }
      if(count <=(1+screenpool_max) & enrolled==N){
        break
      }
      if(count==(1+screenpool_max) & enrolled < N)
      {
        overtime <- 1
        break
        
      }
    }
  } else
  {
    stop("unknown enrichment type, must be either 'quantile', 'quantile_eq', or 'fixed'")
  }
  # make sure we get even number of patients for 1:1 randomization
  if (enrolled%%2!=0) {
    if(enrolled==1)
    {
      block2_enrolled[1,] <- NA
    } else{
      block2_enrolled <- block2_enrolled[1:(enrolled-1),]
    }
  }
  if(all(is.na(block2_enrolled$patientID))) # if nobody is enrolled
  {
    block2_ineligible <- block2_covs[1:(count-1),]
    block2_ineligible <- response_gen(covs = block2_ineligible,dat_description = scenario,randomize = FALSE)
    block2_ineligible$treatment <- 0
    block2_ineligible$response <- rbinom(n=nrow(block2_ineligible),size = 1,prob = block2_ineligible$pC_true)
    block2_screenpool <- block2_ineligible
    block2 <- block1[0,]
    block12 <- subset(previous_blocks,select = colnames(block2_screenpool))
  } else{
    # randomize and get the outcome situation
    block2 <- response_gen(covs = block2_enrolled[!is.na(block2_enrolled$patientID),],dat_description = scenario,randomize = TRUE)
    # for the screened but ineligible patients
    indexes <- setdiff(block2_covs$patientID[1:(count-1)],block2$patientID)
    if(length(indexes)==0) # no enrichment
    {
      block2_screenpool <- block2
    }
    else
    {
      block2_ineligible <- block2_covs[indexes-pID_start+1,]
      block2_ineligible <- response_gen(covs = block2_ineligible,dat_description = scenario,randomize = FALSE)
      block2_ineligible$treatment <- 0
      block2_ineligible$response <- rbinom(n=nrow(block2_ineligible),size = 1,prob = block2_ineligible$pC_true)
      block2_screenpool <- rbind(block2,block2_ineligible)
    }
    block12 <- rbind(subset(previous_blocks,select = colnames(block2)),block2)
  }
  res <- list("next_block"=block2,"accumulated_block"=block12,"screenpool"=block2_screenpool,
              "overtime"=overtime)
  return(res)
}

beta.ineq2 <- function(block_list,delta=0)
{
  block1 <- block_list[[1]]
  block2 <- block_list[[2]]
  if (all(is.na(block2)))
  {
    a <- 1+sum(block1$response[block1$treatment==1])
    b <- 1+sum(1-block1$response[block1$treatment==1])
    c <- 1+sum(block1$response[block1$treatment==0])
    d <- 1+sum(1-block1$response[block1$treatment==0])
  } else{
  #discount_para <- mean(block1$response)/(mean(block1$response)+mean(block2$response))
 # discount_para <- min(1,(mean(block1$response)/mean(block2$response))^2)
  discount_para <- min(1,(mean(block1$response)/mean(block2$response))^2)
  # calculate parameters for Pe
  a <- 1+sum(block1$response[block1$treatment==1])*discount_para + sum(block2$response[block2$treatment==1])
  b <- 1+sum(1-block1$response[block1$treatment==1])*discount_para + sum(1-block2$response[block2$treatment==1])
  c <- 1+sum(block1$response[block1$treatment==0])*discount_para + sum(block2$response[block2$treatment==0])
  d <- 1+sum(1-block1$response[block1$treatment==0])*discount_para + sum(1-block2$response[block2$treatment==0])
  }
  if (a <=0 | b <= 0 | c <=0 | d <= 0) 
    stop("paramters has to be positive")
  if (a <0.01 | b < 0.01 | c <0.01 | d < 0.01) 
    stop("paramters are to close to 0")
  if (delta>1) 
    stop("delta>1!")
  if (delta<0) 
    stop("delta<0!")
  integrand <- function(x) { dbeta(x, a, b)*pbeta(x-delta, c, d) }
  integrate(integrand, delta, 1, rel.tol=1e-4)$value
}

################################################
compute_JSD <- function(dat1,dat2,margin=0.05)
{
  # compute JSD for pE in the upper vs. lower quantile group
  # distributional parameters pE in dat 1
  r_EU <- 1 + sum(dat1$response[dat1$treatment==1])
  s_EU <- 1 + sum(1-dat1$response[dat1$treatment==1])
  # distributional parameters pE in dat 2
  r_EL <- 1 + sum(dat2$response[dat2$treatment==1])
  s_EL <- 1 + sum(1-dat2$response[dat2$treatment==1])
  #  distributional parameters pE in combined data
  dat <- rbind(dat1,dat2)
  r_E <- 1 + sum(dat$response[dat$treatment==1])
  s_E <- 1 + sum(1-dat$response[dat$treatment==1])
  # KL divergence D(p_EU|p_E)
  int <- try(integrate(f=function(x)
  {
    dbeta(x,r_EU,s_EU)*log(dbeta(x,r_EU,s_EU)/dbeta(x,r_E,s_E))
  },
  lower = 0,upper = 1,
  subdivisions = 2000,rel.tol=.Machine$double.eps^.2),silent = TRUE)
  
  if(inherits(int ,'try-error')){
    D_KL_p_EU_p_E <- integrate(f=function(x)
    {
      dbeta(x,r_EU,s_EU)*log(dbeta(x,r_EU,s_EU)/dbeta(x,r_E,s_E))
    },
    lower = margin,upper = 1-margin,
    subdivisions = 2000,rel.tol=.Machine$double.eps^.1)$value
  } else {
    D_KL_p_EU_p_E <- int$value
  }
  # KL divergence D(p_EL|p_E)
  int <- try(integrate(f=function(x)
  {
    dbeta(x,r_EL,s_EL)*log(dbeta(x,r_EL,s_EL)/dbeta(x,r_E,s_E))
  },
  lower = 0,upper = 1,
  subdivisions = 2000,rel.tol=.Machine$double.eps^.2),silent = TRUE)
  if(inherits(int ,'try-error')){
    D_KL_p_EL_p_E <- integrate(f=function(x)
    {
      dbeta(x,r_EL,s_EL)*log(dbeta(x,r_EL,s_EL)/dbeta(x,r_E,s_E))
    },
    lower = margin,upper = 1-margin,
    subdivisions = 2000,rel.tol=.Machine$double.eps^.1)$value
  } else {
    D_KL_p_EL_p_E <- int$value
  }
  
  JSD_p_EU_p_EL <- 0.5*D_KL_p_EU_p_E + 0.5*D_KL_p_EL_p_E
  
  
  # compute JSD for pC in the upper vs. lower quantile group
  # distributional parameters pC in dat 1
  r_CU <- 1 + sum(dat1$response[dat1$treatment==0])
  s_CU <- 1 + sum(1-dat1$response[dat1$treatment==0])
  # distributional parameters pC in dat 2
  r_CL <- 1 + sum(dat2$response[dat2$treatment==0])
  s_CL <- 1 + sum(1-dat2$response[dat2$treatment==0])
  #  distributional parameters pC in combined data
  dat <- rbind(dat1,dat2)
  r_C <- 1 + sum(dat$response[dat$treatment==0])
  s_C <- 1 + sum(1-dat$response[dat$treatment==0])
  # KL divergence D(p_CU|p_C)
  int <- try(integrate(f=function(x)
  {
    dbeta(x,r_CU,s_CU)*log(dbeta(x,r_CU,s_CU)/dbeta(x,r_C,s_C))
  },
  lower = 0,upper = 1,
  subdivisions = 2000,rel.tol=.Machine$double.eps^.2),silent = TRUE)
  if(inherits(int ,'try-error')){
    D_KL_p_CU_p_C <- integrate(f=function(x)
    {
      dbeta(x,r_CU,s_CU)*log(dbeta(x,r_CU,s_CU)/dbeta(x,r_C,s_C))
    },
    lower = margin,upper = 1-margin,
    rel.tol=0.01)$value
  } else {
    D_KL_p_CU_p_C <- int$value
  }
  # KL divergence D(p_CL|p_C)
  int <- try(integrate(f=function(x)
  {
    dbeta(x,r_CL,s_CL)*log(dbeta(x,r_CL,s_CL)/dbeta(x,r_C,s_C))
  },
  lower = 0,upper = 1,
  subdivisions = 2000,rel.tol=.Machine$double.eps^.2),silent = TRUE)
  if(inherits(int ,'try-error')){
    D_KL_p_CL_p_C <- integrate(f=function(x)
    {
      dbeta(x,r_CL,s_CL)*log(dbeta(x,r_CL,s_CL)/dbeta(x,r_C,s_C))
    },
    lower = margin,upper = 1-margin,
    subdivisions = 2000,rel.tol=.Machine$double.eps^.1)$value
  } else {
    D_KL_p_CL_p_C <- int$value
  }
  JSD_p_CU_p_CL <- 0.5*D_KL_p_CU_p_C + 0.5*D_KL_p_CL_p_C
  return(c("JSD-treatment"=JSD_p_EU_p_EL,"JSD-control"=JSD_p_CU_p_CL))
}
