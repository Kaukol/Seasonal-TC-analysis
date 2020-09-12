###Set seed------------------------------------
setwd("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-data")

set.seed(1)

source("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-analysis/gibbssearch.R")

###Search program---------------------------------------

fastsearch <- function(StartPoint){
  
  #Model setting
  StartModel <- StartPoint$StartModelIndex
  InverseT <- StartPoint$InverseT
  Seq <- StartPoint$Seq
  family <- StartPoint$family
  info <- StartPoint$info
  z <- StartPoint$database
  
  return(GibbsSampler_Search_MH(z[,1],z[,-1],StartModel,Seq, InverseT, info = info, family = family))
}

fastsearch_SA <- function(StartPoint){
  
  #Model setting
  StartModel <- StartPoint$StartModelIndex
  InverseT <- StartPoint$InverseT
  Seq <- StartPoint$Seq
  family <- StartPoint$family
  info <- StartPoint$info
  z <- StartPoint$database
  
  return(GibbsSampler_Search_SA(z[,1],z[,-1],StartModel,Seq, InverseT, info = info, family = family))
}

ModelSelected <- function(modelmatrix){return(colSums(BurnSeq(modelmatrix,1000))/1000)}

###Data-------------------------------------------------------------------

###Tropicial cyclone counts
tc.counts <- read.csv("countsrevised.csv")
###covariates data in Aug Sep Oct
x1 <- read.csv("data_aug_1970.csv")
x2 <- read.csv("data_sep_1970.csv")
x3 <- read.csv("data_oct_1970.csv")
x.r <- cbind(x1[,3:14], x2[,3:14], x3[,3:14])

tc <- tc.counts$AR.W

n <- dim(x.r)[1]
p <- dim(x.r)[2]

###Random start model -----------------------------------------------------
Seq_length <- 5000
#Cut_length <- 1000
Loop <- 200

###Parallel
ncores <- 32
cl <- makeCluster(ncores, type = "FORK")

###Linear Model -------------------------------------------------------------------------------------------

InverseT <- 1
####Linear---------------------------------------------------
###AR------------------------------------------------------------------
y <- tc
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(database = z, StartModelIndex = rbinom(p,1,0.5), InverseT = InverseT , Seq = Seq_length, info = "AICc", family = "gaussian")))
  i <- i+1
}

###Searching
SelectedLModel <- parLapply(cl,StartPoint,fastsearch)

SelectedLModel_SA <- parLapply(cl,StartPoint,fastsearch_SA)

####Poisson regression-------------------------------------------
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(database = z, StartModelIndex = rbinom(p,1,0.5), InverseT = InverseT , Seq = Seq_length, info = "AICc", family = "poisson")))
  i <- i+1
}

###Searching
SelectedPModel <- parLapply(cl,StartPoint,fastsearch)

SelectedPModel_SA <- parLapply(cl,StartPoint,fastsearch_SA)

####Linear---------------------------------------------------
InverseT <- 3
###AR------------------------------------------------------------------
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(database = z, StartModelIndex = rbinom(p,1,0.5), InverseT = InverseT , Seq = Seq_length, info = "AICc", family = "gaussian")))
  i <- i+1
}

###Searching
SelectedLModel.3 <- parLapply(cl,StartPoint,fastsearch)

####Poisson regression-------------------------------------------
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(database = z, StartModelIndex = rbinom(p,1,0.5), InverseT = InverseT , Seq = Seq_length, info = "AICc", family = "poisson")))
  i <- i+1
}

###Searching
SelectedPModel.3 <- parLapply(cl,StartPoint,fastsearch)

##END
stopCluster(cl)

LModelCompare <- t(sapply(SelectedLModel, ModelSelected))
summary(LModelCompare)
PModelCompare <- t(sapply(SelectedPModel, ModelSelected))
summary(PModelCompare)
LModelCompare.3 <- t(sapply(SelectedLModel.3, ModelSelected))
summary(LModelCompare.3)
PModelCompare.3 <- t(sapply(SelectedPModel.3, ModelSelected))
summary(PModelCompare.3)
LModelCompare_SA <- t(sapply(SelectedLModel_SA, ModelSelected))
summary(LModelCompare_SA)
PModelCompare_SA <- t(sapply(SelectedPModel_SA, ModelSelected))
summary(PModelCompare_SA)

save.image("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/RData-results/tcarw.RData")
