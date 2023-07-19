library(MASS)
library(Matrix)
library(matrixcalc)
library(matrixStats)
library(mvtnorm)
##################################################
# compute Bayesian credible subgroup using       #
# restricted covariate space method              #
##################################################
# q: dimension of predictive covariates (including main treatment effect)
compute_BCS_RCS <- function(dat,post_coefs,mod,alpha)
{
  post_coefs_sub <- post_coefs[c("treatment","treatment:c1","treatment:c2"),]
  post_mean <- matrix(rep(mod$coefficients[c("treatment","treatment:c1","treatment:c2")],ncol(post_coefs)),
                      ncol=ncol(post_coefs))
  varcov_m <- vcov(mod)[c("treatment","treatment:c1","treatment:c2"),c("treatment","treatment:c1","treatment:c2")]
  
  covs <- matrix(NA,byrow = FALSE,nrow=nrow(dat),ncol=9)
  colnames(covs) <- c("1","c1","c2","TE_se","TE_mean","lower","upper","In_S","In_D")
  covs[,c("c1","c2")] <- as.matrix(dat[,c("c1","c2")])
  covs[,"1"] <- rep(1,nrow(covs))
  covs[,"TE_mean"] <- covs[,1:3] %*% post_mean[,1]
  covs[,"TE_se"] <-  sqrt(rowSums((covs[,1:3] %*% varcov_m) * covs[,1:3]))
  
  MC_dist <- t(abs(covs[,1:3] %*% (post_coefs_sub - post_mean)) / covs[,"TE_se"])
  sups <- rowMaxs(MC_dist)
  MC_alpha <- quantile(sups,1- alpha)
  covs[,"lower"] <- as.vector(covs[,"TE_mean"] - MC_alpha * covs[,"TE_se"])
  covs[,"upper"] <- as.vector(covs[,"TE_mean"] + MC_alpha * covs[,"TE_se"])
  # The S set
  covs[,"In_S"] <- 1*(covs[,"upper"] >= 0)
  # the D set
  covs[,"In_D"] <- 1*(covs[,"lower"] > 0)
  # colMeans(covs)
  return(list("MC_alpha"=MC_alpha,"covs"=covs))
}

#### classify a patient into D set or S set
classify_D_S <- function(BCS,mod,cov)
{
  require(dplyr)
  post_mean <- matrix(mod$coefficients[c("treatment","treatment:c1","treatment:c2")],
                      ncol=1)
  varcov_m <- vcov(mod)[c("treatment","treatment:c1","treatment:c2"),c("treatment","treatment:c1","treatment:c2")]
  x <- matrix(c(1,cov[,"c1"],cov[,"c2"]),nrow=1)
  se <- sqrt(rowSums((x %*% varcov_m) * x))
  upper <- x %*% post_mean + BCS$MC_alpha*se
  lower <- x %*% post_mean - BCS$MC_alpha*se
  In_S <- upper >= 0
  In_D <- lower > 0
  enroll_prob <- case_when (
    In_D ~ 1,
    !In_S ~ 0,
    In_S & !In_D ~ upper/(upper-lower)
  )
  res <- c(lower,upper,In_S,In_D,enroll_prob)
  names(res) <- c("lower","upper","In_S","In_D","enroll_prob")
  return(res)
}

get_nextblock_BCS <- function(BCS,mod,previous_blocks,
                              N,scenario,
                              screenpool_max=500)
{
  # set.seed(i*86)
  pID_start <- max(previous_blocks$patientID)+1
  block2_covs <- covs_gen(N=1e3,dat_description = scenario,
                          pID=101:1100)
  block2_covs <- block2_covs[block2_covs$patientID>=pID_start,]
  count <- 1
  enrolled <- 0
  overtime <- 0
  block2_enrolled <- data.frame(matrix(NA,nrow=N,ncol=(1+length(scenario$`variable names`))))
  colnames(block2_enrolled) <- c("patientID",scenario$`variable names`)
  
  repeat{
    patient <- block2_covs[count,]
    enroll <- rbinom(n=1,size=1,
                     prob=classify_D_S(BCS = BCS,mod = mod,cov = patient)["enroll_prob"])
    
    if(enroll)
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
