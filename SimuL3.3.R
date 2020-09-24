###Set seed------------------------------------
setwd("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-data")

set.seed(1)

source("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-analysis/gibbssearch.R")
source("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-analysis/SimuX.R")

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

###number of observations
n <- 100
###number of parameter
p <- 30
###Gibbs length
MCMC_Length <- 20000
###Random start points
Loop <- 500

###Variable data
x <-  SimuX(n,30,0.9)
###coefficients
beta <- c(seq(2,30,length.out = 15)/15, rep(0,15))

###Parallel
ncores <- 32
cl <- makeCluster(ncores, type = "FORK")

#########Simulation Linear -----------------

##Noise
sigma_0 <- 1
##Response
y <- x %*% beta + sigma_0 * rnorm(n)

z <- as.data.frame(cbind(y,x))
colnames(z)[1] <- "y"

FullModel<- glm(y~., data = z, family = "gaussian")

###Setting up start point
InverseT <- 3

i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(database = z, StartModelIndex = rbinom(p,1,0.5), InverseT = InverseT , Seq = MCMC_Length, info = "BIC", family = "gaussian")))
  i <- i+1
}

SelectedModel <- parLapply(cl,StartPoint,fastsearch)

###Burn Length
Cut_Length <- 10000

###Functions
###Burn sequence function
BurnSeq_Cut <- function(ModelMatrix){return(BurnSeq(ModelMatrix, Cut_Length))}
###BIC value for an index
BICForModel <- function(Model){
  ModelTF <- Model == 1
  return(InfoCritForModel(ModelTF, z, info = "BIC", family = "gaussian"))
}
###Burn coefs
CoefFrequency <- function(BurnCoefs) {return(colSums(BurnCoefs)/nrow(BurnCoefs))}
###Most frequency variable model
MostFreqModel_0.9 <- function(CoefFreq){return(Model <- glm(y~., data = z[,CoefFreq > 0.9], family = "gaussian"))}
MostFreqModel_0.8 <- function(CoefFreq){return(Model <- glm(y~., data = z[,CoefFreq > 0.8], family = "gaussian"))}
MostFreqModel_0.5 <- function(CoefFreq){return(Model <- glm(y~., data = z[,CoefFreq > 0.5], family = "gaussian"))}
###BIC values function
CalBIC <- function(BurnCoefs){return(apply(as.matrix(BurnCoefs),1, BICForModel))}
###Best Model
CalBest <- function(BurnCoefs){
  ModelBIC <- CalBIC(BurnCoefs)
  return(BestModel <- glm(y~., data = z[,BurnCoefs[which.min(ModelBIC),] == 1], family = "gaussian"))
}
###Model BIC frequency
CalBICFreq <- function(ModelBIC){return(table(round(ModelBIC,2))/length(ModelBIC))}

###Each sequence
###Burn sequence in lists
BurnCoefList <- lapply(SelectedModel, BurnSeq_Cut)
###Coef freq list
CoefFreqList <- lapply(BurnCoefList, CoefFrequency)
###Most freq model list
ModelList_0.9 <- lapply(CoefFreqList, MostFreqModel_0.9)
ModelList_0.8 <- lapply(CoefFreqList, MostFreqModel_0.8)
ModelList_0.5 <- lapply(CoefFreqList, MostFreqModel_0.5)
###BIC Model
ModelListBIC <- lapply(BurnCoefList, CalBIC)
###Model frequency list
ModelFrequencyList <- lapply(ModelListBIC, CalBICFreq)
###Best Model list
BestModelList <- lapply(BurnCoefList, CalBest)
###BIC for each point
ModelBICpoints <- lapply(BestModelList, BIC)

###Whole sequence
###Burn sequence in whole
BurnCoefs <- do.call(rbind, BurnCoefList)
###Coef frequency
CoefFreq <- CoefFrequency(BurnCoefs)
###Model by frequency
Model_0.9 <- MostFreqModel_0.9(CoefFreq)
Model_0.8 <- MostFreqModel_0.8(CoefFreq)
Model_0.5 <- MostFreqModel_0.5(CoefFreq)
###BIC values for burn index
ModelBIC <- CalBIC(BurnCoefs)
###Model frequency
ModelFrequency <- CalBICFreq(ModelBIC)
###Best Model
BestModel <- glm(y~., data = z[,BurnCoefs[which.min(ModelBIC),] == 1], family = "gaussian")

stopCluster(cl)

###save data
save.image("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/RData-results/SimuL3.3.RData")

###End

