# description of data to be simulated
scenario0 <- list(
  "scenario name" = "null_scenario",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario0_00 <- list(
  "scenario name" = "scenario0_00",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario0_01 <- list(
  "scenario name" = "scenario0_01",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0.25,0.25),
  "predictive effects"=c(0,0)
)

scenario0_02 <- list(
  "scenario name" = "scenario0_02",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0.5,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0.5,0.5),
  "predictive effects"=c(0,0)
)

scenario0_03 <- list(
  "scenario name" = "scenario0_03",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0.5,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0.75,0.75),
  "predictive effects"=c(0,0)
)

scenario0_04 <- list(
  "scenario name" = "scenario0_04",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0.5,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(1,1),
  "predictive effects"=c(0,0)
)

# description of data to be simulated
scenario1 <- list(
  "scenario name" = "scenario1",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.5,0.5)
)
scenario1_00 <- list(
  "scenario name" = "scenario1_00",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)
scenario1_01 <- list(
  "scenario name" = "scenario1_01",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.1,0.1)
)
scenario1_02 <- list(
  "scenario name" = "scenario1_02",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.2,0.2)
)

scenario1_03 <- list(
  "scenario name" = "scenario1_03",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.3,0.3)
)

scenario1_04 <- list(
  "scenario name" = "scenario1_04",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.4,0.4)
)
scenario1_05 <- list(
  "scenario name" = "scenario1_05",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.5,0.5)
)
scenario1_06 <- list(
  "scenario name" = "scenario1_06",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.6,0.6)
)

# global TE, no subgroup
scenario2 <- list(
  "scenario name" = "scenario2",
  "treatment effect"= 0.5,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario2_00 <- list(
  "scenario name" = "scenario2_00",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario2_01 <- list(
  "scenario name" = "scenario2_01",
  "treatment effect"= 0.25,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario2_02 <- list(
  "scenario name" = "scenario2_02",
  "treatment effect"= 0.5,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario2_03 <- list(
  "scenario name" = "scenario2_03",
  "treatment effect"= 0.75,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario2_04 <- list(
  "scenario name" = "scenario2_04",
  "treatment effect"= 1,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0,0)
)

scenario3 <- list(
  "scenario name" = "scenario3",
  "treatment effect"= 0.2,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.5,0.5)
)

scenario4 <- list(
  "scenario name" = "scenario4",
  "treatment effect"= 0.5,
  "variable names"=c("c1","c2"),
  "binary" = c(FALSE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(0,0),
  "predictive effects"=c(0.1,0.1)
)

scenario5 <- list(
  "scenario name" = "scenario5",
  "treatment effect"= 0,
  "variable names"=c("c1","c2"),
  "binary" = c(TRUE,FALSE),
  "variables means"=c(0,0.5),
  "variable standard deviations"=c(1,1.5),
  "variable correlation" = 0.2,
  "prognostic effects"=c(-1.2,0),
  "predictive effects"=c(1.2,0)
)

