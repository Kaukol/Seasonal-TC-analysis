#Library
#library(dr)
library(glmnet)
#library(SISIR)
#library(plsRglm)
library(AICcmodavg)
#library(ggplot2)
library(doParallel)
library(boot)

###Set seed------------------------------------
set.seed(1)

####
#sample0 <- function(x, ...) x[sample.int(length(x), ...)]
cv.glm1 <- function(data, glmfit, cost=function(y,yhat) mean((y-yhat)^2), K=n){
  n <- nrow(data)
  #    glm.f <- formula(glmfit)
  glm.y <- glmfit$y
  cost.0 <- cost(glm.y, fitted(glmfit))
  seq_len <- 1:n
  CV <- 0
  Hindcast <- vector()
  Call <- glmfit$call
  for(i in 1:n) {
    j.out <- seq_len == i
    j.in <- seq_len != i
    ## we want data from here but formula from the parent.
    Call$data <- data[j.in, , drop=FALSE]
    d.glm <- eval.parent(Call)
    p.alpha <- 1/n
    Hindcast[i] <- predict(d.glm, data[j.out, , drop=FALSE], type = "response")
    cost.i <- cost(glm.y[j.out], Hindcast[i])
    CV <- CV + p.alpha * cost.i
    cost.0 <- cost.0 - p.alpha *
      cost(glm.y, predict(d.glm, data, type = "response"))
  }
  list(K = K,
       delta = as.numeric(c(CV, CV + cost.0)),  # drop any names
       Hindcast = Hindcast)
}

### Recent year forecast
recentHindcast <- function(data, glmfit, cost=function(y,yhat) mean((y-yhat)^2), K = 5){
  call <- match.call()
  n <- nrow(data)
  if ((K > n) || (K <= 1))
    stop("Fail to do hindcast")
  k <- round(K)
  glm.f <- formula(glmfit)
  glm.y <- glmfit$y
  CV <- 0
  Call <- glmfit$call
  Forecast <- vector()
  for(i in 1:k) {
    j.out <- n-k+i
    j.in <- 1:(n-k+i-1)
    ## we want data from here but formula from the parent.
    Call$data <- data[j.in, , drop=FALSE]
    d.glm <- eval.parent(Call)
    p.alpha <- 1/(n-k+1)
    Forecast[i] <- predict(d.glm, data[j.out, , drop=FALSE], type = "response")
    cost.i <- cost(glm.y[j.out], Forecast[i])
    CV <- CV + p.alpha * cost.i
  }
  list(K = K,
       delta = as.numeric(CV),  # drop any names
       Forecast = Forecast)
}

###Search program---------------------------------------
####AICc
fastsearch.aicc <- function(StartPoint){
  
  #Model setting
  StartModel <- StartPoint$StartModelIndex
  Multiplier <- StartPoint$Multiplier
  Seq_Length <- StartPoint$Seq
  Cut_Length <- StartPoint$Cut
  threshold <- StartPoint$threshold
  distribution <- StartPoint$distribution
  z <- StartPoint$database
  
  #sigmoid function with multiplied by n
  Sigmoid <- function(x,n){ return(1/(1+exp(-n*x)))}
  
  #conditional probablity based on AICc
  ConditionalProbabilityAICc <- function(glm1, glm2, n){ return(Sigmoid(AICc(glm2)-AICc(glm1),n))}
  
  #Model selection 
  ModelChoice1 <- function(ModelSeletionMatrix, threshold, Cut_Length){
    CutModelMatrix <- ModelSeletionMatrix[(nrow(ModelSeletionMatrix) - Cut_Length + 1):nrow(ModelSeletionMatrix),]
    return(colSums(CutModelMatrix)/Cut_Length > threshold)
  }
  
  #Gibbs sampler
  GibbsSampler.fs <- function(StartModel, Seq_Length, Multiplier, ConditionalProbability, data, distribution){
    
    #Start
    database <- data
    n <- Multiplier
    StartIndex <- StartModel
    StartModel <- c(1, StartIndex)
    ModelSeletionMatrix <- StartModel
    
    #Loop
    j <- 0
    while(j < Seq_Length){
      for(i in 1:ncol(database) - 1){
        TempIndex <- StartIndex
        TempDual <- TempIndex
        TempDual[i] <- 1 - TempIndex[i]
        
        TempCheck <- c(1,TempIndex) == 1
        TempCheckDual <- c(1, TempDual) == 1
        
        if(sum(!(TempIndex == 0)) == 0){
          y <- database[,TempCheck]
          Model_One <- glm(y ~ 1,family = distribution)
        }else{
          Model_One <- glm(y ~ .,family = distribution, data = database[,TempCheck])
        }
        
        if(sum(!(TempDual == 0)) == 0){
          y <- database[,TempCheckDual]
          Model_Two <- glm(y ~ 1,family = distribution)
        }else{
          Model_Two <- glm(y ~ .,family = distribution, data = database[,TempCheckDual])
        }
        
        #Random jump
        BinomProb <- ConditionalProbability(Model_One, Model_Two, n)
        
        TossUp <- rbinom(1,1,BinomProb)
        if(TossUp == 1){
          TempIndex[i] <- TempIndex[i]
        }
        else{
          TempIndex[i] <- 1 - TempIndex[i]
        }
        StartIndex <- TempIndex
      }
      ModelSeletionMatrix <- rbind(ModelSeletionMatrix, c(1,StartIndex))
      j = j + 1
    }
    return(ModelSeletionMatrix)
  }
  
  #Model Selection
  ModelSeletionMatrix <- GibbsSampler.fs(StartModel, Seq_Length, Multiplier, ConditionalProbabilityAICc, z, distribution)
  SelectedModel <- ModelChoice1(ModelSeletionMatrix, threshold, Cut_Length) 
  return(SelectedModel)
  
}
AICcselect <- function(SelectedModel){ return(AICc(glm(y~., data = z[,SelectedModel], family = distribution)))}
###Parallel
ncores <- 20

###Data-------------------------------------------------------------------
setwd("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/Seasonal-TC-data")
###Tropicial cyclone counts
tc.counts <- read.csv("countsrevised.csv")
###covariates data in Aug Sep Oct
x1 <- read.csv("data_aug_1970.csv")
x2 <- read.csv("data_sep_1970.csv")
x3 <- read.csv("data_oct_1970.csv")
x.r <- cbind(x1[,3:14], x2[,3:14], x3[,3:14])

tc.ar <- tc.counts$AR.1
tc.are <- tc.counts$AR.E1
tc.arw <- tc.counts$AR.W

tc.sp <- tc.counts$SP
tc.spe <- tc.counts$SP.E1
tc.spw <- tc.counts$SP.W1

Year <-tc.counts$Year
n <- dim(x.r)[1]
p <- dim(x.r)[2]

###Random start model -----------------------------------------------------
Multiplier <- 5
Seq_length <- 300
Cut_length <- 100
Loop <- 100
threshold_num <- 0.8

###Linear Model -------------------------------------------------------------------------------------------
distribution <- gaussian

###AICc -----------------------------------------------------------------------------

###AR------------------------------------------------------------------
y <- tc.ar
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedLModel.aicc.ar <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.ar <- SelectedLModel.aicc.ar[,colSums(SelectedLModel.aicc.ar)!=1]

LModel.aicc.ar.table <- apply(SelectedLModel.aicc.ar,2,AICcselect)
table(LModel.aicc.ar.table)

LModel.aicc.ar <- glm(y~., data = z[,SelectedLModel.aicc.ar[,which.min(LModel.aicc.ar.table)]], family = distribution)
summary(LModel.aicc.ar)
stopCluster(cl)

###ARW
y <- tc.arw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel.aicc.arw <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.arw <- SelectedLModel.aicc.arw[,colSums(SelectedLModel.aicc.arw)!=1]

LModel.aicc.arw.table <- apply(SelectedLModel.aicc.arw,2,AICcselect)
table(LModel.aicc.arw.table)

LModel.aicc.arw <- glm(y~., data = z[,SelectedLModel.aicc.arw[,which.min(LModel.aicc.arw.table)]], family = distribution)
summary(LModel.aicc.arw)
stopCluster(cl)

###ARE
y <- tc.are
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel.aicc.are <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.are <- SelectedLModel.aicc.are[,colSums(SelectedLModel.aicc.are)!=1]

LModel.aicc.are.table <- apply(SelectedLModel.aicc.are,2,AICcselect)
table(LModel.aicc.are.table)

LModel.aicc.are <- glm(y~., data = z[,SelectedLModel.aicc.are[,which.min(LModel.aicc.are.table)]], family = distribution)
summary(LModel.aicc.are)
stopCluster(cl)

###SP------------------------------------------------------------------------
y <- tc.sp
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedLModel.aicc.sp <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.sp <- SelectedLModel.aicc.sp[,colSums(SelectedLModel.aicc.sp)!=1]

LModel.aicc.sp.table <- apply(SelectedLModel.aicc.sp,2,AICcselect)
table(LModel.aicc.sp.table)

LModel.aicc.sp <- glm(y~., data = z[,SelectedLModel.aicc.sp[,which.min(LModel.aicc.sp.table)]], family = distribution)
summary(LModel.aicc.sp)
stopCluster(cl)

###SPW
y <- tc.spw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedLModel.aicc.spw <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.spw <- SelectedLModel.aicc.spw[,colSums(SelectedLModel.aicc.spw)!=1]

LModel.aicc.spw.table <- apply(SelectedLModel.aicc.spw,2,AICcselect)
table(LModel.aicc.spw.table)

LModel.aicc.spw <- glm(y~., data = z[,SelectedLModel.aicc.spw[,which.min(LModel.aicc.spw.table)]], family = distribution)
summary(LModel.aicc.spw)
stopCluster(cl)

###SPE
y <- tc.spe
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedLModel.aicc.spe <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.spe <- SelectedLModel.aicc.spe[,colSums(SelectedLModel.aicc.spe)!=1]

LModel.aicc.spe.table <- apply(SelectedLModel.aicc.spe,2,AICcselect)
table(LModel.aicc.spe.table)

LModel.aicc.spe <- glm(y~., data = z[,SelectedLModel.aicc.spe[,which.min(LModel.aicc.spe.table)]], family = distribution)
summary(LModel.aicc.spe)
stopCluster(cl)



####Poisson Model --------------------------------------------------------------------------------------
distribution <- poisson

###AICc -----------------------------------------------------------------------------

###AR------------------------------------------------------------------
y <- tc.ar
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedPModel.aicc.ar <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.ar <- SelectedPModel.aicc.ar[,colSums(SelectedPModel.aicc.ar)!=1]

PModel.aicc.ar.table <- apply(SelectedPModel.aicc.ar,2,AICcselect)
table(PModel.aicc.ar.table)

PModel.aicc.ar <- glm(y~., data = z[,SelectedPModel.aicc.ar[,which.min(PModel.aicc.ar.table)]], family = distribution)
summary(PModel.aicc.ar)
stopCluster(cl)

###ARW
y <- tc.arw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedPModel.aicc.arw <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.arw <- SelectedPModel.aicc.arw[,colSums(SelectedPModel.aicc.arw)!=1]

PModel.aicc.arw.table <- apply(SelectedPModel.aicc.arw,2,AICcselect)
table(PModel.aicc.arw.table)

PModel.aicc.arw <- glm(y~., data = z[,SelectedPModel.aicc.arw[,which.min(PModel.aicc.arw.table)]], family = distribution)
summary(PModel.aicc.arw)
stopCluster(cl)

###ARE
y <- tc.are
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedPModel.aicc.are <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.are <- SelectedPModel.aicc.are[,colSums(SelectedPModel.aicc.are)!=1]

PModel.aicc.are.table <- apply(SelectedPModel.aicc.are,2,AICcselect)
table(PModel.aicc.are.table)

PModel.aicc.are <- glm(y~., data = z[,SelectedPModel.aicc.are[,which.min(PModel.aicc.are.table)]], family = distribution)
summary(PModel.aicc.are)
stopCluster(cl)

###SP------------------------------------------------------------------------
y <- tc.sp
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedPModel.aicc.sp <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.sp <- SelectedPModel.aicc.sp[,colSums(SelectedPModel.aicc.sp)!=1]

PModel.aicc.sp.table <- apply(SelectedPModel.aicc.sp,2,AICcselect)
table(PModel.aicc.sp.table)

PModel.aicc.sp <- glm(y~., data = z[,SelectedPModel.aicc.sp[,which.min(PModel.aicc.sp.table)]], family = distribution)
summary(PModel.aicc.sp)
stopCluster(cl)

###SPW
y <- tc.spw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedPModel.aicc.spw <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.spw <- SelectedPModel.aicc.spw[,colSums(SelectedPModel.aicc.spw)!=1]

PModel.aicc.spw.table <- apply(SelectedPModel.aicc.spw,2,AICcselect)
table(PModel.aicc.spw.table)

PModel.aicc.spw <- glm(y~., data = z[,SelectedPModel.aicc.spw[,which.min(PModel.aicc.spw.table)]], family = distribution)
summary(PModel.aicc.spw)
stopCluster(cl)

###SPE
y <- tc.spe
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))

SelectedPModel.aicc.spe <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.spe <- SelectedPModel.aicc.spe[,colSums(SelectedPModel.aicc.spe)!=1]

PModel.aicc.spe.table <- apply(SelectedPModel.aicc.spe,2,AICcselect)
table(PModel.aicc.spe.table)

PModel.aicc.spe <- glm(y~., data = z[,SelectedPModel.aicc.spe[,which.min(PModel.aicc.spe.table)]], family = distribution)
summary(PModel.aicc.spe)
stopCluster(cl)



###Summary-----------------------------------------------------------------------------------


### LOCCV(Hindcasts)
y <- tc.ar
z <- cbind(y, x.r)
NULModel.ar <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.ar <- cv.glm1(z, NULModel.ar)
LOOCV.mae.NULModel.ar <- cv.glm1(z, NULModel.ar, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.ar <- recentHindcast(z, NULModel.ar)

z1 <- z[, SelectedLModel.aicc.ar[, which.min(LModel.aicc.ar.table)]]
LModel.aicc.ar.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.ar <- cv.glm1(z1, LModel.aicc.ar.simple)
LOOCV.mae.LModel.aicc.ar <- cv.glm1(z1, LModel.aicc.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.ar <- recentHindcast(z1, LModel.aicc.ar.simple)

z1 <- z[, SelectedPModel.aicc.ar[, which.min(PModel.aicc.ar.table)]]
PModel.aicc.ar.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.ar <- cv.glm1(z1, PModel.aicc.ar.simple)
LOOCV.mae.PModel.aicc.ar <- cv.glm1(z1, PModel.aicc.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.ar <- recentHindcast(z1, PModel.aicc.ar.simple)


###ARW------------------------------------------------------------------

### LOCCV(Hindcasts)
y <- tc.arw
z <- cbind(y, x.r)
NULModel.arw <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.arw <- cv.glm1(z, NULModel.arw)
LOOCV.mae.NULModel.arw <- cv.glm1(z, NULModel.arw, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.arw <- recentHindcast(z, NULModel.arw)

z1 <- z[, SelectedLModel.aicc.arw[, which.min(LModel.aicc.arw.table)]]
LModel.aicc.arw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.arw <- cv.glm1(z1, LModel.aicc.arw.simple)
LOOCV.mae.LModel.aicc.arw <- cv.glm1(z1, LModel.aicc.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.arw <- recentHindcast(z1, LModel.aicc.arw.simple)

z1 <- z[, SelectedPModel.aicc.arw[, which.min(PModel.aicc.arw.table)]]
PModel.aicc.arw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.arw <- cv.glm1(z1, PModel.aicc.arw.simple)
LOOCV.mae.PModel.aicc.arw <- cv.glm1(z1, PModel.aicc.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.arw <- recentHindcast(z1, PModel.aicc.arw.simple)

###ARE ----------------------------------------------------------------------------

### LOCCV(Hindcasts)
y <- tc.are
z <- cbind(y, x.r)
NULModel.are <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.are <- cv.glm1(z, NULModel.are)
LOOCV.mae.NULModel.are <- cv.glm1(z, NULModel.are, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.are <- recentHindcast(z, NULModel.are)

z1 <- z[, SelectedLModel.aicc.are[, which.min(LModel.aicc.are.table)]]
LModel.aicc.are.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.are <- cv.glm1(z1, LModel.aicc.are.simple)
LOOCV.mae.LModel.aicc.are <- cv.glm1(z1, LModel.aicc.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.are <- recentHindcast(z1, LModel.aicc.are.simple)

z1 <- z[, SelectedPModel.aicc.are[, which.min(PModel.aicc.are.table)]]
PModel.aicc.are.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.are <- cv.glm1(z1, PModel.aicc.are.simple)
LOOCV.mae.PModel.aicc.are <- cv.glm1(z1, PModel.aicc.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.are <- recentHindcast(z1, PModel.aicc.are.simple)

###SP ------------------------------------------------

### LOCCV(Hindcasts)
y <- tc.sp
z <- cbind(y, x.r)
NULModel.sp <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.sp <- cv.glm1(z, NULModel.sp)
LOOCV.mae.NULModel.sp <- cv.glm1(z, NULModel.sp, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.sp <- recentHindcast(z, NULModel.sp)

z1 <- z[, SelectedLModel.aicc.sp[, which.min(LModel.aicc.sp.table)]]
LModel.aicc.sp.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.sp <- cv.glm1(z1, LModel.aicc.sp.simple)
LOOCV.mae.LModel.aicc.sp <- cv.glm1(z1, LModel.aicc.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.sp <- recentHindcast(z1, LModel.aicc.sp.simple)

z1 <- z[, SelectedPModel.aicc.sp[, which.min(PModel.aicc.sp.table)]]
PModel.aicc.sp.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.sp <- cv.glm1(z1, PModel.aicc.sp.simple)
LOOCV.mae.PModel.aicc.sp <- cv.glm1(z1, PModel.aicc.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.sp <- recentHindcast(z1, PModel.aicc.sp.simple)

###SPW -------------------------------------------------------------------

### LOCCV(Hindcasts)
y <- tc.spw
z <- cbind(y, x.r)
NULModel.spw <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.spw <- cv.glm1(z, NULModel.spw)
LOOCV.mae.NULModel.spw <- cv.glm1(z, NULModel.spw, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.spw <- recentHindcast(z, NULModel.spw)

z1 <- z[, SelectedLModel.aicc.spw[, which.min(LModel.aicc.spw.table)]]
LModel.aicc.spw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.spw <- cv.glm1(z1, LModel.aicc.spw.simple)
LOOCV.mae.LModel.aicc.spw <- cv.glm1(z1, LModel.aicc.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spw <- recentHindcast(z1, LModel.aicc.spw.simple)

z1 <- z[, SelectedPModel.aicc.spw[, which.min(PModel.aicc.spw.table)]]
PModel.aicc.spw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.spw <- cv.glm1(z1, PModel.aicc.spw.simple)
LOOCV.mae.PModel.aicc.spw <- cv.glm1(z1, PModel.aicc.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spw <- recentHindcast(z1, PModel.aicc.spw.simple)

###SPE------------------------------------------------------------------------------

### LOCCV(Hindcasts)
y <- tc.spe
z <- cbind(y, x.r)
NULModel.spe <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.spe <- cv.glm1(z, NULModel.spe)
LOOCV.mae.NULModel.spe <- cv.glm1(z, NULModel.spe, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.spe <- recentHindcast(z, NULModel.spe)

z1 <- z[, SelectedLModel.aicc.spe[, which.min(LModel.aicc.spe.table)]]
LModel.aicc.spe.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.spe <- cv.glm1(z1, LModel.aicc.spe.simple)
LOOCV.mae.LModel.aicc.spe <- cv.glm1(z1, LModel.aicc.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spe <- recentHindcast(z1, LModel.aicc.spe.simple)

z1 <- z[, SelectedPModel.aicc.spe[, which.min(PModel.aicc.spe.table)]]
PModel.aicc.spe.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.PModel.aicc.spe <- cv.glm1(z1, PModel.aicc.spe.simple)
LOOCV.mae.PModel.aicc.spe <- cv.glm1(z1, PModel.aicc.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spe <- recentHindcast(z1, PModel.aicc.spe.simple)

save.image("/mount/autofs/home_ad1/student.unimelb.edu.au/lizhongc/myGit/RData-results/test.RData")

