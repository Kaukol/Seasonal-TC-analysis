#Library
library(dr)
library(glmnet)
library(SISIR)
library(plsRglm)
library(AICcmodavg)
library(ggplot2)
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
    p.alpha <- 1/k
    Forecast[i] <- predict(d.glm, data[j.out, , drop=FALSE], type = "response")
    cost.i <- cost(glm.y[j.out], Forecast[i])
    CV <- CV + p.alpha * cost.i
  }
  list(K = K,
       delta = as.numeric(CV),  # drop any names
       Forecast = Forecast)
}

###Search program---------------------------------------
####AIC
fastsearch.aic <- function(StartPoint){
  
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
  #conditional probablity based on AIC
  ConditionalProbabilityAIC <- function(glm1, glm2, n){ return(Sigmoid(AIC(glm2)-AIC(glm1),n))}
  #conditional probablity based on AIC
  ConditionalProbabilityBIC <- function(glm1, glm2, n){ return(Sigmoid(BIC(glm2)-AIC(glm1),n))}
  #conditional probablity based on AIC
  ConditionalProbabilityAICc <- function(glm1, glm2, n){ return(Sigmoid(AICc(glm2)-AIC(glm1),n))}
  
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
  ModelSeletionMatrix <- GibbsSampler.fs(StartModel, Seq_Length, Multiplier, ConditionalProbabilityAIC, z, distribution)
  SelectedModel <- ModelChoice1(ModelSeletionMatrix, threshold, Cut_Length) 
  return(SelectedModel)
  
}
AICselect <- function(SelectedModel){ return(AIC(glm(y~., data = z[,SelectedModel], family = distribution)))}
####AIC
fastsearch.bic <- function(StartPoint){
  
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

  #conditional probablity based on AIC
  ConditionalProbabilityBIC <- function(glm1, glm2, n){ return(Sigmoid(BIC(glm2)-BIC(glm1),n))}
  
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
  ModelSeletionMatrix <- GibbsSampler.fs(StartModel, Seq_Length, Multiplier, ConditionalProbabilityBIC, z, distribution)
  SelectedModel <- ModelChoice1(ModelSeletionMatrix, threshold, Cut_Length) 
  return(SelectedModel)
  
}
BICselect <- function(SelectedModel){ return(BIC(glm(y~., data = z[,SelectedModel], family = distribution)))}
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
ncores <- detectCores() - 1

###Data-------------------------------------------------------------------
setwd("C:/Users/nealf/OneDrive/My R Work and Data/TCNewData")
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
Loop <- 200
threshold_num <- 0.8


###Linear Model -------------------------------------------------------------------------------------------
distribution <- gaussian

###AIC -------------------------------------------------------------------------------------

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
SelectedLModel.aic.ar <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.ar <- SelectedLModel.aic.ar[,colSums(SelectedLModel.aic.ar)!=1]

LModel.aic.ar.table <- apply(SelectedLModel.aic.ar,2,AICselect)
table(LModel.aic.ar.table)

LModel.aic.ar <- glm(y~., data = z[,SelectedLModel.aic.ar[,which.min(LModel.aic.ar.table)]], family = distribution)
summary(LModel.aic.ar)
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
SelectedLModel.aic.arw <- parSapply(cl,StartPoint,fastsearch.aic)

##Exclude null model
SelectedLModel.aic.arw <- SelectedLModel.aic.arw[,colSums(SelectedLModel.aic.arw)!=1]

LModel.aic.arw.table <- apply(SelectedLModel.aic.arw,2,AICselect)
table(LModel.aic.arw.table)

LModel.aic.arw <- glm(y~., data = z[,SelectedLModel.aic.arw[,which.min(LModel.aic.arw.table)]], family = distribution)
summary(LModel.aic.arw)
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
SelectedLModel.aic.are <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.are <- SelectedLModel.aic.are[,colSums(SelectedLModel.aic.are)!=1]

LModel.aic.are.table <- apply(SelectedLModel.aic.are,2,AICselect)
table(LModel.aic.are.table)

LModel.aic.are <- glm(y~., data = z[,SelectedLModel.aic.are[,which.min(LModel.aic.are.table)]], family = distribution)
summary(LModel.aic.are)
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
SelectedLModel.aic.sp <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.sp <- SelectedLModel.aic.sp[,colSums(SelectedLModel.aic.sp)!=1]

LModel.aic.sp.table <- apply(SelectedLModel.aic.sp,2,AICselect)
table(LModel.aic.sp.table)

LModel.aic.sp <- glm(y~., data = z[,SelectedLModel.aic.sp[,which.min(LModel.aic.sp.table)]], family = distribution)
summary(LModel.aic.sp)
stopCluster(cl)

###SPW
y <- tc.spw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:(Loop+500)){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.aic.spw <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.spw <- SelectedLModel.aic.spw[,colSums(SelectedLModel.aic.spw)!=1]

LModel.aic.spw.table <- apply(SelectedLModel.aic.spw,2,AICselect)
table(LModel.aic.spw.table)

LModel.aic.spw <- glm(y~., data = z[,SelectedLModel.aic.spw[,which.min(LModel.aic.spw.table)]], family = distribution)
summary(LModel.aic.spw)
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
SelectedLModel.aic.spe <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.spe <- SelectedLModel.aic.spe[,colSums(SelectedLModel.aic.spe)!=1]

LModel.aic.spe.table <- apply(SelectedLModel.aic.spe,2,AICselect)
table(LModel.aic.spe.table)

LModel.aic.spe <- glm(y~., data = z[,SelectedLModel.aic.spe[,which.min(LModel.aic.spe.table)]], family = distribution)
summary(LModel.aic.spe)
stopCluster(cl)



###BIC ---------------------------------------------------------------------------------

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
SelectedLModel.bic.ar <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.ar <- SelectedLModel.bic.ar[,colSums(SelectedLModel.bic.ar)!=1]

LModel.bic.ar.table <- apply(SelectedLModel.bic.ar,2,BICselect)
table(LModel.bic.ar.table)

LModel.bic.ar <- glm(y~., data = z[,SelectedLModel.bic.ar[,which.min(LModel.bic.ar.table)]], family = distribution)
summary(LModel.bic.ar)
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
SelectedLModel.bic.arw <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.arw <- SelectedLModel.bic.arw[,colSums(SelectedLModel.bic.arw)!=1]

LModel.bic.arw.table <- apply(SelectedLModel.bic.arw,2,BICselect)
table(LModel.bic.arw.table)

LModel.bic.arw <- glm(y~., data = z[,SelectedLModel.bic.arw[,which.min(LModel.bic.arw.table)]], family = distribution)
summary(LModel.bic.arw)
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
SelectedLModel.bic.are <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.are <- SelectedLModel.bic.are[,colSums(SelectedLModel.bic.are)!=1]

LModel.bic.are.table <- apply(SelectedLModel.bic.are,2,BICselect)
table(LModel.bic.are.table)

LModel.bic.are <- glm(y~., data = z[,SelectedLModel.bic.are[,which.min(LModel.bic.are.table)]], family = distribution)
summary(LModel.bic.are)
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
SelectedLModel.bic.sp <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.sp <- SelectedLModel.bic.sp[,colSums(SelectedLModel.bic.sp)!=1]

LModel.bic.sp.table <- apply(SelectedLModel.bic.sp,2,BICselect)
table(LModel.bic.sp.table)

LModel.bic.sp <- glm(y~., data = z[,SelectedLModel.bic.sp[,which.min(LModel.bic.sp.table)]], family = distribution)
summary(LModel.bic.sp)
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
SelectedLModel.bic.spw <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.spw <- SelectedLModel.bic.spw[,colSums(SelectedLModel.bic.spw)!=1]

LModel.bic.spw.table <- apply(SelectedLModel.bic.spw,2,BICselect)
table(LModel.bic.spw.table)

LModel.bic.spw <- glm(y~., data = z[,SelectedLModel.bic.spw[,which.min(LModel.bic.spw.table)]], family = distribution)
summary(LModel.bic.spw)
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
SelectedLModel.bic.spe <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.spe <- SelectedLModel.bic.spe[,colSums(SelectedLModel.bic.spe)!=1]

LModel.bic.spe.table <- apply(SelectedLModel.bic.spe,2,BICselect)
table(LModel.bic.spe.table)

LModel.bic.spe <- glm(y~., data = z[,SelectedLModel.bic.spe[,which.min(LModel.bic.spe.table)]], family = distribution)
summary(LModel.bic.spe)
stopCluster(cl)



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

###AIC -------------------------------------------------------------------------------------

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
SelectedPModel.aic.ar <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.ar <- SelectedPModel.aic.ar[,colSums(SelectedPModel.aic.ar)!=1]

PModel.aic.ar.table <- apply(SelectedPModel.aic.ar,2,AICselect)
table(PModel.aic.ar.table)

PModel.aic.ar <- glm(y~., data = z[,SelectedPModel.aic.ar[,which.min(PModel.aic.ar.table)]], family = distribution)
summary(PModel.aic.ar)
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
SelectedPModel.aic.arw <- parSapply(cl,StartPoint,fastsearch.aic)

##Exclude null model
SelectedPModel.aic.arw <- SelectedPModel.aic.arw[,colSums(SelectedPModel.aic.arw)!=1]

PModel.aic.arw.table <- apply(SelectedPModel.aic.arw,2,AICselect)
table(PModel.aic.arw.table)

PModel.aic.arw <- glm(y~., data = z[,SelectedPModel.aic.arw[,which.min(PModel.aic.arw.table)]], family = distribution)
summary(PModel.aic.arw)
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
SelectedPModel.aic.are <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.are <- SelectedPModel.aic.are[,colSums(SelectedPModel.aic.are)!=1]

PModel.aic.are.table <- apply(SelectedPModel.aic.are,2,AICselect)
table(PModel.aic.are.table)

PModel.aic.are <- glm(y~., data = z[,SelectedPModel.aic.are[,which.min(PModel.aic.are.table)]], family = distribution)
summary(PModel.aic.are)
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
SelectedPModel.aic.sp <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.sp <- SelectedPModel.aic.sp[,colSums(SelectedPModel.aic.sp)!=1]

PModel.aic.sp.table <- apply(SelectedPModel.aic.sp,2,AICselect)
table(PModel.aic.sp.table)

PModel.aic.sp <- glm(y~., data = z[,SelectedPModel.aic.sp[,which.min(PModel.aic.sp.table)]], family = distribution)
summary(PModel.aic.sp)
stopCluster(cl)

###SPW
y <- tc.spw
z <- as.data.frame(cbind(y,x.r))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:(Loop+500)){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.aic.spw <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.spw <- SelectedPModel.aic.spw[,colSums(SelectedPModel.aic.spw)!=1]

PModel.aic.spw.table <- apply(SelectedPModel.aic.spw,2,AICselect)
table(PModel.aic.spw.table)

PModel.aic.spw <- glm(y~., data = z[,SelectedPModel.aic.spw[,which.min(PModel.aic.spw.table)]], family = distribution)
summary(PModel.aic.spw)
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
SelectedPModel.aic.spe <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.spe <- SelectedPModel.aic.spe[,colSums(SelectedPModel.aic.spe)!=1]

PModel.aic.spe.table <- apply(SelectedPModel.aic.spe,2,AICselect)
table(PModel.aic.spe.table)

PModel.aic.spe <- glm(y~., data = z[,SelectedPModel.aic.spe[,which.min(PModel.aic.spe.table)]], family = distribution)
summary(PModel.aic.spe)
stopCluster(cl)



###BIC ----------------------------------------------------------------------------------------------

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
SelectedPModel.bic.ar <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.ar <- SelectedPModel.bic.ar[,colSums(SelectedPModel.bic.ar)!=1]

PModel.bic.ar.table <- apply(SelectedPModel.bic.ar,2,BICselect)
table(PModel.bic.ar.table)

PModel.bic.ar <- glm(y~., data = z[,SelectedPModel.bic.ar[,which.min(PModel.bic.ar.table)]], family = distribution)
summary(PModel.bic.ar)
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
SelectedPModel.bic.arw <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.arw <- SelectedPModel.bic.arw[,colSums(SelectedPModel.bic.arw)!=1]

PModel.bic.arw.table <- apply(SelectedPModel.bic.arw,2,BICselect)
table(PModel.bic.arw.table)

PModel.bic.arw <- glm(y~., data = z[,SelectedPModel.bic.arw[,which.min(PModel.bic.arw.table)]], family = distribution)
summary(PModel.bic.arw)
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
SelectedPModel.bic.are <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.are <- SelectedPModel.bic.are[,colSums(SelectedPModel.bic.are)!=1]

PModel.bic.are.table <- apply(SelectedPModel.bic.are,2,BICselect)
table(PModel.bic.are.table)

PModel.bic.are <- glm(y~., data = z[,SelectedPModel.bic.are[,which.min(PModel.bic.are.table)]], family = distribution)
summary(PModel.bic.are)
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
SelectedPModel.bic.sp <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.sp <- SelectedPModel.bic.sp[,colSums(SelectedPModel.bic.sp)!=1]

PModel.bic.sp.table <- apply(SelectedPModel.bic.sp,2,BICselect)
table(PModel.bic.sp.table)

PModel.bic.sp <- glm(y~., data = z[,SelectedPModel.bic.sp[,which.min(PModel.bic.sp.table)]], family = distribution)
summary(PModel.bic.sp)
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
SelectedPModel.bic.spw <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.spw <- SelectedPModel.bic.spw[,colSums(SelectedPModel.bic.spw)!=1]

PModel.bic.spw.table <- apply(SelectedPModel.bic.spw,2,BICselect)
table(PModel.bic.spw.table)

PModel.bic.spw <- glm(y~., data = z[,SelectedPModel.bic.spw[,which.min(PModel.bic.spw.table)]], family = distribution)
summary(PModel.bic.spw)
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
SelectedPModel.bic.spe <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.spe <- SelectedPModel.bic.spe[,colSums(SelectedPModel.bic.spe)!=1]

PModel.bic.spe.table <- apply(SelectedPModel.bic.spe,2,BICselect)
table(PModel.bic.spe.table)

PModel.bic.spe <- glm(y~., data = z[,SelectedPModel.bic.spe[,which.min(PModel.bic.spe.table)]], family = distribution)
summary(PModel.bic.spe)
stopCluster(cl)




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

###AR--------------------------------------------

table(LModel.aic.ar.table)
summary(LModel.aic.ar)
table(LModel.bic.ar.table)
summary(LModel.bic.ar)
table(LModel.aicc.ar.table)
summary(LModel.aicc.ar)

table(PModel.aic.ar.table)
summary(PModel.aic.ar)
table(PModel.aic.ar.table)
summary(PModel.bic.ar)
table(PModel.aic.ar.table)
summary(PModel.aicc.ar)

LInfo.ar <- rbind(c(AIC(LModel.aic.ar),AIC(LModel.bic.ar),AIC(LModel.aicc.ar),AIC(PModel.aic.ar),AIC(PModel.bic.ar),AIC(PModel.aicc.ar)), 
                  c(BIC(LModel.aic.ar),BIC(LModel.bic.ar),BIC(LModel.aicc.ar),BIC(PModel.aic.ar),BIC(PModel.bic.ar),BIC(PModel.aicc.ar)),
                  c(AICc(LModel.aic.ar),AICc(LModel.bic.ar),AICc(LModel.aicc.ar),AICc(PModel.aic.ar),AICc(PModel.bic.ar),AICc(PModel.aicc.ar)))

rownames(LInfo.ar) <- c("AIC","BIC","AICc")
colnames(LInfo.ar) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.ar

###Fitting plot
plot(Year, tc.ar, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.ar), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.ar), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.ar), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.ar, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.ar), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.ar), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.ar), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.ar
z <- cbind(y, x.r)
NULModel.ar <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.ar <- cv.glm1(z, NULModel.ar)
LOOCV.mae.NULModel.ar <- cv.glm1(z, NULModel.ar, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.ar <- recentHindcast(z, NULModel.ar)

z1 <- z[, SelectedLModel.aic.ar[, which.min(LModel.aic.ar.table)]]
LModel.aic.ar.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.ar <- cv.glm1(z1, LModel.aic.ar.simple)
LOOCV.mae.LModel.aic.ar <- cv.glm1(z1, LModel.aic.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.ar <- recentHindcast(z1, LModel.aic.ar.simple)

z1 <- z[, SelectedLModel.bic.ar[, which.min(LModel.bic.ar.table)]]
LModel.bic.ar.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.ar <- cv.glm1(z1, LModel.bic.ar.simple)
LOOCV.mae.LModel.bic.ar <- cv.glm1(z1, LModel.bic.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.ar <- recentHindcast(z1, LModel.bic.ar.simple)

z1 <- z[, SelectedLModel.aicc.ar[, which.min(LModel.aicc.ar.table)]]
LModel.aicc.ar.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.ar <- cv.glm1(z1, LModel.aicc.ar.simple)
LOOCV.mae.LModel.aicc.ar <- cv.glm1(z1, LModel.aicc.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.ar <- recentHindcast(z1, LModel.aicc.ar.simple)

z1 <- z[, SelectedPModel.aic.ar[, which.min(PModel.aic.ar.table)]]
PModel.aic.ar.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.ar <- cv.glm1(z1, PModel.aic.ar.simple)
LOOCV.mae.PModel.aic.ar <- cv.glm1(z1, PModel.aic.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.ar <- recentHindcast(z1, PModel.aic.ar.simple)

z1 <- z[, SelectedPModel.bic.ar[, which.min(PModel.bic.ar.table)]]
PModel.bic.ar.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.ar <- cv.glm1(z1, PModel.bic.ar.simple)
LOOCV.mae.PModel.bic.ar <- cv.glm1(z1, PModel.bic.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.ar <- recentHindcast(z1, PModel.bic.ar.simple)

z1 <- z[, SelectedPModel.aicc.ar[, which.min(PModel.aicc.ar.table)]]
PModel.aicc.ar.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.ar <- cv.glm1(z1, PModel.aicc.ar.simple)
LOOCV.mae.PModel.aicc.ar <- cv.glm1(z1, PModel.aicc.ar.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.ar <- recentHindcast(z1, PModel.aicc.ar.simple)

Hindcast.mse.ar <- cbind(LOOCV.mse.NULModel.ar$delta, LOOCV.mse.LModel.aic.ar$delta,LOOCV.mse.LModel.bic.ar$delta,LOOCV.mse.LModel.aicc.ar$delta,
                     LOOCV.mse.PModel.aic.ar$delta,LOOCV.mse.PModel.bic.ar$delta,LOOCV.mse.PModel.aicc.ar$delta)
rownames(Hindcast.mse.ar) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.ar) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.ar

Hindcast.mae.ar <- cbind(LOOCV.mae.NULModel.ar$delta, LOOCV.mae.LModel.aic.ar$delta,LOOCV.mae.LModel.bic.ar$delta,LOOCV.mae.LModel.aicc.ar$delta,
                     LOOCV.mae.PModel.aic.ar$delta,LOOCV.mae.PModel.bic.ar$delta,LOOCV.mae.PModel.aicc.ar$delta)
rownames(Hindcast.mae.ar) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.ar) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.ar

Forecast.ar <- cbind(Forecast.NULModel.ar$delta, Forecast.LModel.aic.ar$delta,Forecast.LModel.bic.ar$delta,Forecast.LModel.aicc.ar$delta,
                     Forecast.PModel.aic.ar$delta,Forecast.PModel.bic.ar$delta,Forecast.PModel.aicc.ar$delta)
rownames(Forecast.ar) <- c("MSE")
colnames(Forecast.ar) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.ar

###Forecast plot
plot(Year, tc.ar, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.ar$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.ar$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.ar$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.ar, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.ar$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.ar$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.ar$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.ar[44:48], ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.ar$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.ar$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.ar$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.ar[44:48], ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.ar$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.ar$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.ar$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###ARW------------------------------------------------------------------

table(LModel.aic.arw.table)
summary(LModel.aic.arw)
table(LModel.bic.arw.table)
summary(LModel.bic.arw)
table(LModel.aicc.arw.table)
summary(LModel.aicc.arw)

table(PModel.aic.arw.table)
summary(PModel.aic.arw)
table(PModel.bic.arw.table)
summary(PModel.bic.arw)
table(PModel.aicc.arw.table)
summary(PModel.aicc.arw)

LInfo.arw <- rbind(c(AIC(LModel.aic.arw),AIC(LModel.bic.arw),AIC(LModel.aicc.arw),AIC(PModel.aic.arw),AIC(PModel.bic.arw),AIC(PModel.aicc.arw)), 
                   c(BIC(LModel.aic.arw),BIC(LModel.bic.arw),BIC(LModel.aicc.arw),BIC(PModel.aic.arw),BIC(PModel.bic.arw),BIC(PModel.aicc.arw)),
                   c(AICc(LModel.aic.arw),AICc(LModel.bic.arw),AICc(LModel.aicc.arw),AICc(PModel.aic.arw),AICc(PModel.bic.arw),AICc(PModel.aicc.arw)))

rownames(LInfo.arw) <- c("AIC","BIC","AICc")
colnames(LInfo.arw) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.arw)

plot(Year, tc.arw, ylab = "TC counts", main = "ARW LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.arw), col = c("red"), lwd = 2)
 lines(Year, fitted.values(LModel.bic.arw), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.arw), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.arw, ylab = "TC counts", main = "ARW PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.arw), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.arw), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.arw), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.arw
z <- cbind(y, x.r)
NULModel.arw <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.arw <- cv.glm1(z, NULModel.arw)
LOOCV.mae.NULModel.arw <- cv.glm1(z, NULModel.arw, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.arw <- recentHindcast(z, NULModel.arw)

z1 <- z[, SelectedLModel.aic.arw[, which.min(LModel.aic.arw.table)]]
LModel.aic.arw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.arw <- cv.glm1(z1, LModel.aic.arw.simple)
LOOCV.mae.LModel.aic.arw <- cv.glm1(z1, LModel.aic.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.arw <- recentHindcast(z1, LModel.aic.arw.simple)

z1 <- z[, SelectedLModel.bic.arw[, which.min(LModel.bic.arw.table)]]
LModel.bic.arw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.arw <- cv.glm1(z1, LModel.bic.arw.simple)
LOOCV.mae.LModel.bic.arw <- cv.glm1(z1, LModel.bic.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.arw <- recentHindcast(z1, LModel.bic.arw.simple)

z1 <- z[, SelectedLModel.aicc.arw[, which.min(LModel.aicc.arw.table)]]
LModel.aicc.arw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.arw <- cv.glm1(z1, LModel.aicc.arw.simple)
LOOCV.mae.LModel.aicc.arw <- cv.glm1(z1, LModel.aicc.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.arw <- recentHindcast(z1, LModel.aicc.arw.simple)

z1 <- z[, SelectedPModel.aic.arw[, which.min(PModel.aic.arw.table)]]
PModel.aic.arw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.arw <- cv.glm1(z1, PModel.aic.arw.simple)
LOOCV.mae.PModel.aic.arw <- cv.glm1(z1, PModel.aic.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.arw <- recentHindcast(z1, PModel.aic.arw.simple)

z1 <- z[, SelectedPModel.bic.arw[, which.min(PModel.bic.arw.table)]]
PModel.bic.arw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.arw <- cv.glm1(z1, PModel.bic.arw.simple)
LOOCV.mae.PModel.bic.arw <- cv.glm1(z1, PModel.bic.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.arw <- recentHindcast(z1, PModel.bic.arw.simple)

z1 <- z[, SelectedPModel.aicc.arw[, which.min(PModel.aicc.arw.table)]]
PModel.aicc.arw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.arw <- cv.glm1(z1, PModel.aicc.arw.simple)
LOOCV.mae.PModel.aicc.arw <- cv.glm1(z1, PModel.aicc.arw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.arw <- recentHindcast(z1, PModel.aicc.arw.simple)

Hindcast.mse.arw <- cbind(LOOCV.mse.NULModel.arw$delta, LOOCV.mse.LModel.aic.arw$delta,LOOCV.mse.LModel.bic.arw$delta,LOOCV.mse.LModel.aicc.arw$delta,
                          LOOCV.mse.PModel.aic.arw$delta,LOOCV.mse.PModel.bic.arw$delta,LOOCV.mse.PModel.aicc.arw$delta)
rownames(Hindcast.mse.arw) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.arw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.arw

Hindcast.mae.arw <- cbind(LOOCV.mae.NULModel.arw$delta, LOOCV.mae.LModel.aic.arw$delta,LOOCV.mae.LModel.bic.arw$delta,LOOCV.mae.LModel.aicc.arw$delta,
                          LOOCV.mae.PModel.aic.arw$delta,LOOCV.mae.PModel.bic.arw$delta,LOOCV.mae.PModel.aicc.arw$delta)
rownames(Hindcast.mae.arw) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.arw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.arw

Forecast.arw <- cbind(Forecast.NULModel.arw$delta, Forecast.LModel.aic.arw$delta,Forecast.LModel.bic.arw$delta,Forecast.LModel.aicc.arw$delta,
                      Forecast.PModel.aic.arw$delta,Forecast.PModel.bic.arw$delta,Forecast.PModel.aicc.arw$delta)
rownames(Forecast.arw) <- c("MSE")
colnames(Forecast.arw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.arw

###Forecast plot
plot(Year, tc.arw, ylab = "TC counts", main = "ARW LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.arw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.arw$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.arw$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.arw, ylab = "TC counts", main = "ARW PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.arw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.arw$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.arw$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.arw[44:48], ylab = "TC counts", main = "ARW LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.arw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.arw$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.arw$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.arw[44:48], ylab = "TC counts", main = "ARW PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.arw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.arw$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.arw$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)
###ARE ----------------------------------------------------------------------------

table(LModel.aic.are.table)
summary(LModel.aic.are)
table(LModel.bic.are.table)
summary(LModel.bic.are)
table(LModel.aicc.are.table)
summary(LModel.aicc.are)

table(PModel.aic.are.table)
summary(PModel.aic.are)
table(PModel.bic.are.table)
summary(PModel.bic.are)
table(PModel.aicc.are.table)
summary(PModel.aicc.are)

LInfo.are <- rbind(c(AIC(LModel.aic.are),AIC(LModel.bic.are),AIC(LModel.aicc.are),AIC(PModel.aic.are),AIC(PModel.bic.are),AIC(PModel.aicc.are)), 
                   c(BIC(LModel.aic.are),BIC(LModel.bic.are),BIC(LModel.aicc.are),BIC(PModel.aic.are),BIC(PModel.bic.are),BIC(PModel.aicc.are)),
                   c(AICc(LModel.aic.are),AICc(LModel.bic.are),AICc(LModel.aicc.are),AICc(PModel.aic.are),AICc(PModel.bic.are),AICc(PModel.aicc.are)))

rownames(LInfo.are) <- c("AIC","BIC","AICc")
colnames(LInfo.are) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.are)

plot(Year, tc.are, ylab = "TC counts", main = "ARE LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.are), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.are), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.are), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.are, ylab = "TC counts", main = "ARE PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.are), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.are), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.are), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.are
z <- cbind(y, x.r)
NULModel.are <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.are <- cv.glm1(z, NULModel.are)
LOOCV.mae.NULModel.are <- cv.glm1(z, NULModel.are, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.are <- recentHindcast(z, NULModel.are)

z1 <- z[, SelectedLModel.aic.are[, which.min(LModel.aic.are.table)]]
LModel.aic.are.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.are <- cv.glm1(z1, LModel.aic.are.simple)
LOOCV.mae.LModel.aic.are <- cv.glm1(z1, LModel.aic.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.are <- recentHindcast(z1, LModel.aic.are.simple)

z1 <- z[, SelectedLModel.bic.are[, which.min(LModel.bic.are.table)]]
LModel.bic.are.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.are <- cv.glm1(z1, LModel.bic.are.simple)
LOOCV.mae.LModel.bic.are <- cv.glm1(z1, LModel.bic.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.are <- recentHindcast(z1, LModel.bic.are.simple)

z1 <- z[, SelectedLModel.aicc.are[, which.min(LModel.aicc.are.table)]]
LModel.aicc.are.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.are <- cv.glm1(z1, LModel.aicc.are.simple)
LOOCV.mae.LModel.aicc.are <- cv.glm1(z1, LModel.aicc.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.are <- recentHindcast(z1, LModel.aicc.are.simple)

z1 <- z[, SelectedPModel.aic.are[, which.min(PModel.aic.are.table)]]
PModel.aic.are.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.are <- cv.glm1(z1, PModel.aic.are.simple)
LOOCV.mae.PModel.aic.are <- cv.glm1(z1, PModel.aic.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.are <- recentHindcast(z1, PModel.aic.are.simple)

z1 <- z[, SelectedPModel.bic.are[, which.min(PModel.bic.are.table)]]
PModel.bic.are.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.are <- cv.glm1(z1, PModel.bic.are.simple)
LOOCV.mae.PModel.bic.are <- cv.glm1(z1, PModel.bic.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.are <- recentHindcast(z1, PModel.bic.are.simple)

z1 <- z[, SelectedPModel.aicc.are[, which.min(PModel.aicc.are.table)]]
PModel.aicc.are.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.are <- cv.glm1(z1, PModel.aicc.are.simple)
LOOCV.mae.PModel.aicc.are <- cv.glm1(z1, PModel.aicc.are.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.are <- recentHindcast(z1, PModel.aicc.are.simple)

Hindcast.mse.are <- cbind(LOOCV.mse.NULModel.are$delta, LOOCV.mse.LModel.aic.are$delta,LOOCV.mse.LModel.bic.are$delta,LOOCV.mse.LModel.aicc.are$delta,
                          LOOCV.mse.PModel.aic.are$delta,LOOCV.mse.PModel.bic.are$delta,LOOCV.mse.PModel.aicc.are$delta)
rownames(Hindcast.mse.are) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.are) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.are

Hindcast.mae.are <- cbind(LOOCV.mae.NULModel.are$delta, LOOCV.mae.LModel.aic.are$delta,LOOCV.mae.LModel.bic.are$delta,LOOCV.mae.LModel.aicc.are$delta,
                          LOOCV.mae.PModel.aic.are$delta,LOOCV.mae.PModel.bic.are$delta,LOOCV.mae.PModel.aicc.are$delta)
rownames(Hindcast.mae.are) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.are) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.are

Forecast.are <- cbind(Forecast.NULModel.are$delta, Forecast.LModel.aic.are$delta,Forecast.LModel.bic.are$delta,Forecast.LModel.aicc.are$delta,
                      Forecast.PModel.aic.are$delta,Forecast.PModel.bic.are$delta,Forecast.PModel.aicc.are$delta)
rownames(Forecast.are) <- c("MSE")
colnames(Forecast.are) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.are

###Forecast plot
plot(Year, tc.are, ylab = "TC counts", main = "ARE LR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.are$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.are$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.are$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.are, ylab = "TC counts", main = "ARE PR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.are$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.are$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.are$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.are[44:48], ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.are$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.are$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.are$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.are[44:48], ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.are$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.are$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.are$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###SP ------------------------------------------------

summary(LModel.aic.sp)
summary(LModel.bic.sp)
summary(LModel.aicc.sp)

summary(PModel.aic.sp)
summary(PModel.bic.sp)
summary(PModel.aicc.sp)

LInfo.sp <- rbind(c(AIC(LModel.aic.sp),AIC(LModel.bic.sp),AIC(LModel.aicc.sp),AIC(PModel.aic.sp),AIC(PModel.bic.sp),AIC(PModel.aicc.sp)), 
                  c(BIC(LModel.aic.sp),BIC(LModel.bic.sp),BIC(LModel.aicc.sp),BIC(PModel.aic.sp),BIC(PModel.bic.sp),BIC(PModel.aicc.sp)),
                  c(AICc(LModel.aic.sp),AICc(LModel.bic.sp),AICc(LModel.aicc.sp),AICc(PModel.aic.sp),AICc(PModel.bic.sp),AICc(PModel.aicc.sp)))

rownames(LInfo.sp) <- c("AIC","BIC","AICc")
colnames(LInfo.sp) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.sp)

plot(Year, tc.sp, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.sp), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.sp), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.sp), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.sp, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.sp), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.sp), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.sp), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.sp
z <- cbind(y, x.r)
NULModel.sp <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.sp <- cv.glm1(z, NULModel.sp)
LOOCV.mae.NULModel.sp <- cv.glm1(z, NULModel.sp, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.sp <- recentHindcast(z, NULModel.sp)

z1 <- z[, SelectedLModel.aic.sp[, which.min(LModel.aic.sp.table)]]
LModel.aic.sp.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.sp <- cv.glm1(z1, LModel.aic.sp.simple)
LOOCV.mae.LModel.aic.sp <- cv.glm1(z1, LModel.aic.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.sp <- recentHindcast(z1, LModel.aic.sp.simple)

z1 <- z[, SelectedLModel.bic.sp[, which.min(LModel.bic.sp.table)]]
LModel.bic.sp.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.sp <- cv.glm1(z1, LModel.bic.sp.simple)
LOOCV.mae.LModel.bic.sp <- cv.glm1(z1, LModel.bic.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.sp <- recentHindcast(z1, LModel.bic.sp.simple)

z1 <- z[, SelectedLModel.aicc.sp[, which.min(LModel.aicc.sp.table)]]
LModel.aicc.sp.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.sp <- cv.glm1(z1, LModel.aicc.sp.simple)
LOOCV.mae.LModel.aicc.sp <- cv.glm1(z1, LModel.aicc.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.sp <- recentHindcast(z1, LModel.aicc.sp.simple)

z1 <- z[, SelectedPModel.aic.sp[, which.min(PModel.aic.sp.table)]]
PModel.aic.sp.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.sp <- cv.glm1(z1, PModel.aic.sp.simple)
LOOCV.mae.PModel.aic.sp <- cv.glm1(z1, PModel.aic.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.sp <- recentHindcast(z1, PModel.aic.sp.simple)

z1 <- z[, SelectedPModel.bic.sp[, which.min(PModel.bic.sp.table)]]
PModel.bic.sp.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.sp <- cv.glm1(z1, PModel.bic.sp.simple)
LOOCV.mae.PModel.bic.sp <- cv.glm1(z1, PModel.bic.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.sp <- recentHindcast(z1, PModel.bic.sp.simple)

z1 <- z[, SelectedPModel.aicc.sp[, which.min(PModel.aicc.sp.table)]]
PModel.aicc.sp.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.sp <- cv.glm1(z1, PModel.aicc.sp.simple)
LOOCV.mae.PModel.aicc.sp <- cv.glm1(z1, PModel.aicc.sp.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.sp <- recentHindcast(z1, PModel.aicc.sp.simple)

Hindcast.mse.sp <- cbind(LOOCV.mse.NULModel.sp$delta, LOOCV.mse.LModel.aic.sp$delta,LOOCV.mse.LModel.bic.sp$delta,LOOCV.mse.LModel.aicc.sp$delta,
                         LOOCV.mse.PModel.aic.sp$delta,LOOCV.mse.PModel.bic.sp$delta,LOOCV.mse.PModel.aicc.sp$delta)
rownames(Hindcast.mse.sp) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.sp) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.sp

Hindcast.mae.sp <- cbind(LOOCV.mae.NULModel.sp$delta, LOOCV.mae.LModel.aic.sp$delta,LOOCV.mae.LModel.bic.sp$delta,LOOCV.mae.LModel.aicc.sp$delta,
                         LOOCV.mae.PModel.aic.sp$delta,LOOCV.mae.PModel.bic.sp$delta,LOOCV.mae.PModel.aicc.sp$delta)
rownames(Hindcast.mae.sp) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.sp) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.sp

Forecast.sp <- cbind(Forecast.NULModel.sp$delta, Forecast.LModel.aic.sp$delta,Forecast.LModel.bic.sp$delta,Forecast.LModel.aicc.sp$delta,
                     Forecast.PModel.aic.sp$delta,Forecast.PModel.bic.sp$delta,Forecast.PModel.aicc.sp$delta)
rownames(Forecast.sp) <- c("MSE")
colnames(Forecast.sp) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.sp

###Forecast plot
plot(Year, tc.sp, ylab = "TC counts", main = "SP LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.sp$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.sp$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.sp$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.sp, ylab = "TC counts", main = "SP PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.sp$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.sp$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.sp$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.sp[44:48], ylab = "TC counts", main = "SP LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.sp$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.sp$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.sp$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.sp[44:48], ylab = "TC counts", main = "SP PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.sp$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.sp$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.sp$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###SPW -------------------------------------------------------------------

table(LModel.aic.spw.table)
summary(LModel.aic.spw)
table(LModel.bic.spw.table)
summary(LModel.bic.spw)
table(LModel.aicc.spw.table)
summary(LModel.aicc.spw)

table(PModel.aic.spw.table)
summary(PModel.aic.spw)
table(PModel.bic.spw.table)
summary(PModel.bic.spw)
table(PModel.aicc.spw.table)
summary(PModel.aicc.spw)


LInfo.spw <- rbind(c(AIC(LModel.aic.spw),AIC(LModel.bic.spw),AIC(LModel.aicc.spw),AIC(PModel.aic.spw),AIC(PModel.bic.spw),AIC(PModel.aicc.spw)), 
                   c(BIC(LModel.aic.spw),BIC(LModel.bic.spw),BIC(LModel.aicc.spw),BIC(PModel.aic.spw),BIC(PModel.bic.spw),BIC(PModel.aicc.spw)),
                   c(AICc(LModel.aic.spw),AICc(LModel.bic.spw),AICc(LModel.aicc.spw),AICc(PModel.aic.spw),AICc(PModel.bic.spw),AICc(PModel.aicc.spw)))

rownames(LInfo.spw) <- c("AIC","BIC","AICc")
colnames(LInfo.spw) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.spw)

plot(Year, tc.spw, ylab = "TC counts", main = "ARW LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.spw), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.spw), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.spw), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.spw, ylab = "TC counts", main = "ARW PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.spw), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.spw), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.spw), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.spw
z <- cbind(y, x.r)
NULModel.spw <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.spw <- cv.glm1(z, NULModel.spw)
LOOCV.mae.NULModel.spw <- cv.glm1(z, NULModel.spw, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.spw <- recentHindcast(z, NULModel.spw)

z1 <- z[, SelectedLModel.aic.spw[, which.min(LModel.aic.spw.table)]]
LModel.aic.spw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.spw <- cv.glm1(z1, LModel.aic.spw.simple)
LOOCV.mae.LModel.aic.spw <- cv.glm1(z1, LModel.aic.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.spw <- recentHindcast(z1, LModel.aic.spw.simple)

z1 <- z[, SelectedLModel.bic.spw[, which.min(LModel.bic.spw.table)]]
LModel.bic.spw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.spw <- cv.glm1(z1, LModel.bic.spw.simple)
LOOCV.mae.LModel.bic.spw <- cv.glm1(z1, LModel.bic.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.spw <- recentHindcast(z1, LModel.bic.spw.simple)

z1 <- z[, SelectedLModel.aicc.spw[, which.min(LModel.aicc.spw.table)]]
LModel.aicc.spw.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.spw <- cv.glm1(z1, LModel.aicc.spw.simple)
LOOCV.mae.LModel.aicc.spw <- cv.glm1(z1, LModel.aicc.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spw <- recentHindcast(z1, LModel.aicc.spw.simple)

z1 <- z[, SelectedPModel.aic.spw[, which.min(PModel.aic.spw.table)]]
PModel.aic.spw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.spw <- cv.glm1(z1, PModel.aic.spw.simple)
LOOCV.mae.PModel.aic.spw <- cv.glm1(z1, PModel.aic.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.spw <- recentHindcast(z1, PModel.aic.spw.simple)

z1 <- z[, SelectedPModel.bic.spw[, which.min(PModel.bic.spw.table)]]
PModel.bic.spw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.spw <- cv.glm1(z1, PModel.bic.spw.simple)
LOOCV.mae.PModel.bic.spw <- cv.glm1(z1, PModel.bic.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.spw <- recentHindcast(z1, PModel.bic.spw.simple)

z1 <- z[, SelectedPModel.aicc.spw[, which.min(PModel.aicc.spw.table)]]
PModel.aicc.spw.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.spw <- cv.glm1(z1, PModel.aicc.spw.simple)
LOOCV.mae.PModel.aicc.spw <- cv.glm1(z1, PModel.aicc.spw.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spw <- recentHindcast(z1, PModel.aicc.spw.simple)

Hindcast.mse.spw <- cbind(LOOCV.mse.NULModel.spw$delta, LOOCV.mse.LModel.aic.spw$delta,LOOCV.mse.LModel.bic.spw$delta,LOOCV.mse.LModel.aicc.spw$delta,
                          LOOCV.mse.PModel.aic.spw$delta,LOOCV.mse.PModel.bic.spw$delta,LOOCV.mse.PModel.aicc.spw$delta)
rownames(Hindcast.mse.spw) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.spw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.spw

Hindcast.mae.spw <- cbind(LOOCV.mae.NULModel.spw$delta, LOOCV.mae.LModel.aic.spw$delta,LOOCV.mae.LModel.bic.spw$delta,LOOCV.mae.LModel.aicc.spw$delta,
                          LOOCV.mae.PModel.aic.spw$delta,LOOCV.mae.PModel.bic.spw$delta,LOOCV.mae.PModel.aicc.spw$delta)
rownames(Hindcast.mae.spw) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.spw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.spw

Forecast.spw <- cbind(Forecast.NULModel.spw$delta, Forecast.LModel.aic.spw$delta,Forecast.LModel.bic.spw$delta,Forecast.LModel.aicc.spw$delta,
                      Forecast.PModel.aic.spw$delta,Forecast.PModel.bic.spw$delta,Forecast.PModel.aicc.spw$delta)
rownames(Forecast.spw) <- c("MSE")
colnames(Forecast.spw) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.spw

###Forecast plot
plot(Year, tc.spw, ylab = "TC counts", main = "SPW LR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.spw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.spw$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.spw$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.spw, ylab = "TC counts", main = "SPW PR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.spw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.spw$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.spw$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.spw[44:48], ylab = "TC counts", main = "SPW LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.spw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.spw$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.spw$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.spw[44:48], ylab = "TC counts", main = "SPW PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.spw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.spw$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.spw$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)


###SPE------------------------------------------------------------------------------

table(LModel.aic.spe.table)
summary(LModel.aic.spe)
table(LModel.bic.spe.table)
summary(LModel.bic.spe)
table(LModel.aicc.spe.table)
summary(LModel.aicc.spe)

table(PModel.aic.spe.table)
summary(PModel.aic.spe)
table(PModel.bic.spe.table)
summary(PModel.bic.spe)
table(PModel.aicc.spe.table)
summary(PModel.aicc.spe)


LInfo.spe <- rbind(c(AIC(LModel.aic.spe),AIC(LModel.bic.spe),AIC(LModel.aicc.spe),AIC(PModel.aic.spe),AIC(PModel.bic.spe),AIC(PModel.aicc.spe)), 
                   c(BIC(LModel.aic.spe),BIC(LModel.bic.spe),BIC(LModel.aicc.spe),BIC(PModel.aic.spe),BIC(PModel.bic.spe),BIC(PModel.aicc.spe)),
                   c(AICc(LModel.aic.spe),AICc(LModel.bic.spe),AICc(LModel.aicc.spe),AICc(PModel.aic.spe),AICc(PModel.bic.spe),AICc(PModel.aicc.spe)))

rownames(LInfo.spe) <- c("AIC","BIC","AICc")
colnames(LInfo.spe) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.spe

plot(Year, tc.spe, ylab = "TC counts", main = "ARE LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.spe), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.spe), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.spe), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.spe, ylab = "TC counts", main = "ARE PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.spe), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.spe), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.spe), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.spe
z <- cbind(y, x.r)
NULModel.spe <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.spe <- cv.glm1(z, NULModel.spe)
LOOCV.mae.NULModel.spe <- cv.glm1(z, NULModel.spe, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.spe <- recentHindcast(z, NULModel.spe)

z1 <- z[, SelectedLModel.aic.spe[, which.min(LModel.aic.spe.table)]]
LModel.aic.spe.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.spe <- cv.glm1(z1, LModel.aic.spe.simple)
LOOCV.mae.LModel.aic.spe <- cv.glm1(z1, LModel.aic.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.spe <- recentHindcast(z1, LModel.aic.spe.simple)

z1 <- z[, SelectedLModel.bic.spe[, which.min(LModel.bic.spe.table)]]
LModel.bic.spe.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.spe <- cv.glm1(z1, LModel.bic.spe.simple)
LOOCV.mae.LModel.bic.spe <- cv.glm1(z1, LModel.bic.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.spe <- recentHindcast(z1, LModel.bic.spe.simple)

z1 <- z[, SelectedLModel.aicc.spe[, which.min(LModel.aicc.spe.table)]]
LModel.aicc.spe.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.spe <- cv.glm1(z1, LModel.aicc.spe.simple)
LOOCV.mae.LModel.aicc.spe <- cv.glm1(z1, LModel.aicc.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spe <- recentHindcast(z1, LModel.aicc.spe.simple)

z1 <- z[, SelectedPModel.aic.spe[, which.min(PModel.aic.spe.table)]]
PModel.aic.spe.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.spe <- cv.glm1(z1, PModel.aic.spe.simple)
LOOCV.mae.PModel.aic.spe <- cv.glm1(z1, PModel.aic.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.spe <- recentHindcast(z1, PModel.aic.spe.simple)

z1 <- z[, SelectedPModel.bic.spe[, which.min(PModel.bic.spe.table)]]
PModel.bic.spe.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.spe <- cv.glm1(z1, PModel.bic.spe.simple)
LOOCV.mae.PModel.bic.spe <- cv.glm1(z1, PModel.bic.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.spe <- recentHindcast(z1, PModel.bic.spe.simple)

z1 <- z[, SelectedPModel.aicc.spe[, which.min(PModel.aicc.spe.table)]]
PModel.aicc.spe.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.spe <- cv.glm1(z1, PModel.aicc.spe.simple)
LOOCV.mae.PModel.aicc.spe <- cv.glm1(z1, PModel.aicc.spe.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spe <- recentHindcast(z1, PModel.aicc.spe.simple)

Hindcast.mse.spe <- cbind(LOOCV.mse.NULModel.spe$delta, LOOCV.mse.LModel.aic.spe$delta,LOOCV.mse.LModel.bic.spe$delta,LOOCV.mse.LModel.aicc.spe$delta,
                          LOOCV.mse.PModel.aic.spe$delta,LOOCV.mse.PModel.bic.spe$delta,LOOCV.mse.PModel.aicc.spe$delta)
rownames(Hindcast.mse.spe) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.spe) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.spe

Hindcast.mae.spe <- cbind(LOOCV.mae.NULModel.spe$delta, LOOCV.mae.LModel.aic.spe$delta,LOOCV.mae.LModel.bic.spe$delta,LOOCV.mae.LModel.aicc.spe$delta,
                          LOOCV.mae.PModel.aic.spe$delta,LOOCV.mae.PModel.bic.spe$delta,LOOCV.mae.PModel.aicc.spe$delta)
rownames(Hindcast.mae.spe) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.spe) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.spe

Forecast.spe <- cbind(Forecast.NULModel.spe$delta, Forecast.LModel.aic.spe$delta,Forecast.LModel.bic.spe$delta,Forecast.LModel.aicc.spe$delta,
                      Forecast.PModel.aic.spe$delta,Forecast.PModel.bic.spe$delta,Forecast.PModel.aicc.spe$delta)
rownames(Forecast.spe) <- c("MSE")
colnames(Forecast.spe) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.spe

###Forecast plot
plot(Year, tc.spe, ylab = "TC counts", main = "SPE LR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.spe$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.spe$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.spe$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.spe, ylab = "TC counts", main = "SPE PR", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.spe$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.spe$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.spe$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.spe[44:48], ylab = "TC counts", main = "SPE LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.LModel.aic.spe$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.LModel.bic.spe$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.LModel.aicc.spe$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[44:48], tc.spe[44:48], ylab = "TC counts", main = "SPE PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[44:48], Forecast.PModel.aic.spe$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.bic.spe$Forecast, col = c("blue"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.spe$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)




###Best Model------------------------------------------

##AR

Hindcast.mse.ar
Hindcast.mae.ar
Forecast.ar

####Linear with AICc is best

table(LModel.aicc.ar.table)
summary(LModel.aicc.ar.simple)

table(PModel.aicc.ar.table)
summary(PModel.aicc.ar.simple)

#This model has nine variables and it selects same the same variables with AIC
#The poisson model contains less variables, only three for aicc

##ARW
Hindcast.mse.arw
Hindcast.mae.arw
Forecast.arw

####Linear with AICc is best

table(LModel.aicc.arw.table)
summary(LModel.aicc.arw.simple)

table(PModel.aicc.arw.table)
summary(PModel.aicc.arw.simple)

#The best model contains 7 variables
##The poisson model contains five variables

##ARE

Hindcast.mse.are
Hindcast.mae.are
Forecast.are

####Linear with AICc is best

table(LModel.aicc.are.table)
summary(LModel.aicc.are.simple)

table(PModel.aicc.are.table)
summary(PModel.aicc.are.simple)

##SP

Hindcast.mse.sp
Hindcast.mae.sp
Forecast.sp

####Linear with AICc is best

table(LModel.aicc.sp.table)
summary(LModel.aicc.sp.simple)

table(PModel.aicc.sp.table)
summary(PModel.aicc.sp.simple)

#This model has nine variables and it selects same the same variables with AIC
#The poisson model contains less variables, only three for aicc

##ARW
Hindcast.mse.spw
Hindcast.mae.spw
Forecast.spw

####Linear with AICc is best

table(LModel.aicc.spw.table)
summary(LModel.aicc.spw.simple)

table(PModel.aicc.spw.table)
summary(PModel.aicc.spw.simple)

#The best model contains 7 variables
##The poisson model contains five variables

##ARE

Hindcast.mse.spe
Hindcast.mae.spe
Forecast.spe

####Linear with AICc is best

table(LModel.aicc.spe.table)
summary(LModel.aicc.spe.simple)

table(PModel.aicc.spe.table)
summary(PModel.aicc.spe.simple)
#The best model contains six variables
##Poisson contains three variables

###Model average-----------------------------------------

##AR-------------------------------- 
distribution <- gaussian
#AICc Linear
y <- tc.ar
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.ar[, which.min(LModel.aicc.ar.table)]
SelectedLModel.aicc.ar.second <- SelectedLModel.aicc.ar[,LModel.aicc.ar.table != LModel.aicc.ar.table[which.min(LModel.aicc.ar.table)] ]
LModel.aicc.ar.second.table <- apply(SelectedLModel.aicc.ar.second,2,AICcselect)
M2 <- SelectedLModel.aicc.ar.second[, which.min(LModel.aicc.ar.second.table)]
SelectedLModel.aicc.ar.third <- SelectedLModel.aicc.ar.second[,LModel.aicc.ar.second.table != LModel.aicc.ar.second.table[which.min(LModel.aicc.ar.second.table)] ]
LModel.aicc.ar.third.table <- apply(SelectedLModel.aicc.ar.third,2,AICcselect)
M3 <- SelectedLModel.aicc.ar.third[, which.min(LModel.aicc.ar.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.ar.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.ar.1 <- cv.glm1(z.1, LModel.aicc.ar.simple.1)
LOOCV.mae.LModel.aicc.ar.1 <- cv.glm1(z.1, LModel.aicc.ar.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.ar.1 <- recentHindcast(z.1, LModel.aicc.ar.simple.1)

LModel.aicc.ar.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.ar.2 <- cv.glm1(z.2, LModel.aicc.ar.simple.2)
LOOCV.mae.LModel.aicc.ar.2 <- cv.glm1(z.2, LModel.aicc.ar.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.ar.2 <- recentHindcast(z.2, LModel.aicc.ar.simple.2)

LModel.aicc.ar.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.ar.3 <- cv.glm1(z.3, LModel.aicc.ar.simple.3)
LOOCV.mae.LModel.aicc.ar.3 <- cv.glm1(z.3, LModel.aicc.ar.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.ar.3 <- recentHindcast(z.3, LModel.aicc.ar.simple.3)

LOOCV.LModel.aicc.ar.avg <- (LOOCV.mse.LModel.aicc.ar.1$Hindcast + LOOCV.mse.LModel.aicc.ar.2$Hindcast + LOOCV.mse.LModel.aicc.ar.3$Hindcast)/3
LOOCV.mse.LModel.aicc.ar.avg <- mean((tc.ar - LOOCV.LModel.aicc.ar.avg )^2)
Forecast.LModel.aicc.ar.avg <- (Forecast.LModel.aicc.ar.1$Forecast+Forecast.LModel.aicc.ar.2$Forecast+Forecast.LModel.aicc.ar.3$Forecast)/3
Forecast.mse.LModel.aicc.ar.avg <- mean((tc.ar[44:48]- Forecast.LModel.aicc.ar.avg)^2)

cbind(LOOCV.LModel.aicc.ar$delta[1],LOOCV.mse.LModel.aicc.ar.1$delta[1],LOOCV.mse.LModel.aicc.ar.2$delta[1],LOOCV.mse.LModel.aicc.ar.3$delta[1], LOOCV.mse.LModel.aicc.ar.avg )
cbind(Forecast.LModel.aicc.ar$delta,Forecast.LModel.aicc.ar.1$delta,Forecast.LModel.aicc.ar.2$delta,Forecast.LModel.aicc.ar.3$delta,Forecast.mse.LModel.aicc.ar.avg)

LModel.aicc.ar.simple.1
LModel.aicc.ar.simple.2
LModel.aicc.ar.simple.3

#The model averging does not improve the performance: 1, the AIC for three models are almost same, 2, they have similar number of variables.

##ARW-----------------------------------------
distribution <- gaussian
#AICc Linear
y <- tc.arw
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.arw[, which.min(LModel.aicc.arw.table)]
SelectedLModel.aicc.arw.second <- SelectedLModel.aicc.arw[,LModel.aicc.arw.table != LModel.aicc.arw.table[which.min(LModel.aicc.arw.table)] ]
LModel.aicc.arw.second.table <- apply(SelectedLModel.aicc.arw.second,2,AICcselect)
M2 <- SelectedLModel.aicc.arw.second[, which.min(LModel.aicc.arw.second.table)]
SelectedLModel.aicc.arw.third <- SelectedLModel.aicc.arw.second[,LModel.aicc.arw.second.table != LModel.aicc.arw.second.table[which.min(LModel.aicc.arw.second.table)] ]
LModel.aicc.arw.third.table <- apply(SelectedLModel.aicc.arw.third,2,AICcselect)
M3 <- SelectedLModel.aicc.arw.third[, which.min(LModel.aicc.arw.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.arw.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.arw.1 <- cv.glm1(z.1, LModel.aicc.arw.simple.1)
LOOCV.mae.LModel.aicc.arw.1 <- cv.glm1(z.1, LModel.aicc.arw.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.arw.1 <- recentHindcast(z.1, LModel.aicc.arw.simple.1)

LModel.aicc.arw.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.arw.2 <- cv.glm1(z.2, LModel.aicc.arw.simple.2)
LOOCV.mae.LModel.aicc.arw.2 <- cv.glm1(z.2, LModel.aicc.arw.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.arw.2 <- recentHindcast(z.2, LModel.aicc.arw.simple.2)

LModel.aicc.arw.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.arw.3 <- cv.glm1(z.3, LModel.aicc.arw.simple.3)
LOOCV.mae.LModel.aicc.arw.3 <- cv.glm1(z.3, LModel.aicc.arw.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.arw.3 <- recentHindcast(z.3, LModel.aicc.arw.simple.3)

LOOCV.LModel.aicc.arw.avg <- (LOOCV.mse.LModel.aicc.arw.1$Hindcast + LOOCV.mse.LModel.aicc.arw.2$Hindcast + LOOCV.mse.LModel.aicc.arw.3$Hindcast)/3
LOOCV.mse.LModel.aicc.arw.avg <- mean((tc.arw - LOOCV.LModel.aicc.arw.avg )^2)
Forecast.LModel.aicc.arw.avg <- (Forecast.LModel.aicc.arw.1$Forecast+Forecast.LModel.aicc.arw.2$Forecast+Forecast.LModel.aicc.arw.3$Forecast)/3
Forecast.mse.LModel.aicc.arw.avg <- mean((tc.arw[44:48]- Forecast.LModel.aicc.arw.avg)^2)

cbind(LOOCV.LModel.aicc.arw$delta[1],LOOCV.mse.LModel.aicc.arw.1$delta[1],LOOCV.mse.LModel.aicc.arw.2$delta[1],LOOCV.mse.LModel.aicc.arw.3$delta[1], LOOCV.mse.LModel.aicc.arw.avg )
cbind(Forecast.LModel.aicc.arw$delta,Forecast.LModel.aicc.arw.1$delta,Forecast.LModel.aicc.arw.2$delta,Forecast.LModel.aicc.arw.3$delta,Forecast.mse.LModel.aicc.arw.avg)

LModel.aicc.arw.simple.1
LModel.aicc.arw.simple.2
LModel.aicc.arw.simple.3

##For ARW, the model averaging improve slightly

##ARE---------------------------------------
distribution <- gaussian
#AICc Linear
y <- tc.are
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.are[, which.min(LModel.aicc.are.table)]
SelectedLModel.aicc.are.second <- SelectedLModel.aicc.are[,LModel.aicc.are.table != LModel.aicc.are.table[which.min(LModel.aicc.are.table)] ]
LModel.aicc.are.second.table <- apply(SelectedLModel.aicc.are.second,2,AICcselect)
M2 <- SelectedLModel.aicc.are.second[, which.min(LModel.aicc.are.second.table)]
SelectedLModel.aicc.are.third <- SelectedLModel.aicc.are.second[,LModel.aicc.are.second.table != LModel.aicc.are.second.table[which.min(LModel.aicc.are.second.table)] ]
LModel.aicc.are.third.table <- apply(SelectedLModel.aicc.are.third,2,AICcselect)
M3 <- SelectedLModel.aicc.are.third[, which.min(LModel.aicc.are.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.are.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.are.1 <- cv.glm1(z.1, LModel.aicc.are.simple.1)
LOOCV.mae.LModel.aicc.are.1 <- cv.glm1(z.1, LModel.aicc.are.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.are.1 <- recentHindcast(z.1, LModel.aicc.are.simple.1)

LModel.aicc.are.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.are.2 <- cv.glm1(z.2, LModel.aicc.are.simple.2)
LOOCV.mae.LModel.aicc.are.2 <- cv.glm1(z.2, LModel.aicc.are.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.are.2 <- recentHindcast(z.2, LModel.aicc.are.simple.2)

LModel.aicc.are.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.are.3 <- cv.glm1(z.3, LModel.aicc.are.simple.3)
LOOCV.mae.LModel.aicc.are.3 <- cv.glm1(z.3, LModel.aicc.are.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.are.3 <- recentHindcast(z.3, LModel.aicc.are.simple.3)

LOOCV.LModel.aicc.are.avg <- (LOOCV.mse.LModel.aicc.are.1$Hindcast + LOOCV.mse.LModel.aicc.are.2$Hindcast + LOOCV.mse.LModel.aicc.are.3$Hindcast)/3
LOOCV.mse.LModel.aicc.are.avg <- mean((tc.are - LOOCV.LModel.aicc.are.avg )^2)
Forecast.LModel.aicc.are.avg <- (Forecast.LModel.aicc.are.1$Forecast+Forecast.LModel.aicc.are.2$Forecast+Forecast.LModel.aicc.are.3$Forecast)/3
Forecast.mse.LModel.aicc.are.avg <- mean((tc.are[44:48]- Forecast.LModel.aicc.are.avg)^2)

cbind(LOOCV.LModel.aicc.are$delta[1],LOOCV.mse.LModel.aicc.are.1$delta[1],LOOCV.mse.LModel.aicc.are.2$delta[1],LOOCV.mse.LModel.aicc.are.3$delta[1], LOOCV.mse.LModel.aicc.are.avg )
cbind(Forecast.LModel.aicc.are$delta,Forecast.LModel.aicc.are.1$delta,Forecast.LModel.aicc.are.2$delta,Forecast.LModel.aicc.are.3$delta,Forecast.mse.LModel.aicc.are.avg)

LModel.aicc.are.simple.1
LModel.aicc.are.simple.2
LModel.aicc.are.simple.3

#The best models contains six variables and the second best contains only three



##Poisson

##AR-------------------------------- 
distribution <- poisson
#AICc Linear
y <- tc.ar
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.ar[, which.min(PModel.aicc.ar.table)]
SelectedPModel.aicc.ar.second <- SelectedPModel.aicc.ar[,PModel.aicc.ar.table != PModel.aicc.ar.table[which.min(PModel.aicc.ar.table)] ]
PModel.aicc.ar.second.table <- apply(SelectedPModel.aicc.ar.second,2,AICcselect)
M2 <- SelectedPModel.aicc.ar.second[, which.min(PModel.aicc.ar.second.table)]
SelectedPModel.aicc.ar.third <- SelectedPModel.aicc.ar.second[,PModel.aicc.ar.second.table != PModel.aicc.ar.second.table[which.min(PModel.aicc.ar.second.table)] ]
PModel.aicc.ar.third.table <- apply(SelectedPModel.aicc.ar.third,2,AICcselect)
M3 <- SelectedPModel.aicc.ar.third[, which.min(PModel.aicc.ar.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.ar.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.ar.1 <- cv.glm1(z.1, PModel.aicc.ar.simple.1)
LOOCV.mae.PModel.aicc.ar.1 <- cv.glm1(z.1, PModel.aicc.ar.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.ar.1 <- recentHindcast(z.1, PModel.aicc.ar.simple.1)

PModel.aicc.ar.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.ar.2 <- cv.glm1(z.2, PModel.aicc.ar.simple.2)
LOOCV.mae.PModel.aicc.ar.2 <- cv.glm1(z.2, PModel.aicc.ar.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.ar.2 <- recentHindcast(z.2, PModel.aicc.ar.simple.2)

PModel.aicc.ar.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.ar.3 <- cv.glm1(z.3, PModel.aicc.ar.simple.3)
LOOCV.mae.PModel.aicc.ar.3 <- cv.glm1(z.3, PModel.aicc.ar.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.ar.3 <- recentHindcast(z.3, PModel.aicc.ar.simple.3)

LOOCV.PModel.aicc.ar.avg <- (LOOCV.mse.PModel.aicc.ar.1$Hindcast + LOOCV.mse.PModel.aicc.ar.2$Hindcast + LOOCV.mse.PModel.aicc.ar.3$Hindcast)/3
LOOCV.mse.PModel.aicc.ar.avg <- mean((tc.ar - LOOCV.PModel.aicc.ar.avg )^2)
Forecast.PModel.aicc.ar.avg <- (Forecast.PModel.aicc.ar.1$Forecast+Forecast.PModel.aicc.ar.2$Forecast+Forecast.PModel.aicc.ar.3$Forecast)/3
Forecast.mse.PModel.aicc.ar.avg <- mean((tc.ar[44:48]- Forecast.PModel.aicc.ar.avg)^2)

cbind(LOOCV.PModel.aicc.ar$delta[1],LOOCV.mse.PModel.aicc.ar.1$delta[1],LOOCV.mse.PModel.aicc.ar.2$delta[1],LOOCV.mse.PModel.aicc.ar.3$delta[1], LOOCV.mse.PModel.aicc.ar.avg )
cbind(Forecast.PModel.aicc.ar$delta,Forecast.PModel.aicc.ar.1$delta,Forecast.PModel.aicc.ar.2$delta,Forecast.PModel.aicc.ar.3$delta,Forecast.mse.PModel.aicc.ar.avg)

PModel.aicc.ar.simple.1
PModel.aicc.ar.simple.2
PModel.aicc.ar.simple.3

#The model averging does not improve the performance: 1, the AIC for three models are almost same, 2, they have similar number of variables.

##ARW-----------------------------------------
distribution <- poisson
#AICc Linear
y <- tc.arw
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.arw[, which.min(PModel.aicc.arw.table)]
SelectedPModel.aicc.arw.second <- SelectedPModel.aicc.arw[,PModel.aicc.arw.table != PModel.aicc.arw.table[which.min(PModel.aicc.arw.table)] ]
PModel.aicc.arw.second.table <- apply(SelectedPModel.aicc.arw.second,2,AICcselect)
M2 <- SelectedPModel.aicc.arw.second[, which.min(PModel.aicc.arw.second.table)]
SelectedPModel.aicc.arw.third <- SelectedPModel.aicc.arw.second[,PModel.aicc.arw.second.table != PModel.aicc.arw.second.table[which.min(PModel.aicc.arw.second.table)] ]
PModel.aicc.arw.third.table <- apply(SelectedPModel.aicc.arw.third,2,AICcselect)
M3 <- SelectedPModel.aicc.arw.third[, which.min(PModel.aicc.arw.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.arw.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.arw.1 <- cv.glm1(z.1, PModel.aicc.arw.simple.1)
LOOCV.mae.PModel.aicc.arw.1 <- cv.glm1(z.1, PModel.aicc.arw.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.arw.1 <- recentHindcast(z.1, PModel.aicc.arw.simple.1)

PModel.aicc.arw.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.arw.2 <- cv.glm1(z.2, PModel.aicc.arw.simple.2)
LOOCV.mae.PModel.aicc.arw.2 <- cv.glm1(z.2, PModel.aicc.arw.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.arw.2 <- recentHindcast(z.2, PModel.aicc.arw.simple.2)

PModel.aicc.arw.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.arw.3 <- cv.glm1(z.3, PModel.aicc.arw.simple.3)
LOOCV.mae.PModel.aicc.arw.3 <- cv.glm1(z.3, PModel.aicc.arw.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.arw.3 <- recentHindcast(z.3, PModel.aicc.arw.simple.3)

LOOCV.PModel.aicc.arw.avg <- (LOOCV.mse.PModel.aicc.arw.1$Hindcast + LOOCV.mse.PModel.aicc.arw.2$Hindcast + LOOCV.mse.PModel.aicc.arw.3$Hindcast)/3
LOOCV.mse.PModel.aicc.arw.avg <- mean((tc.arw - LOOCV.PModel.aicc.arw.avg )^2)
Forecast.PModel.aicc.arw.avg <- (Forecast.PModel.aicc.arw.1$Forecast+Forecast.PModel.aicc.arw.2$Forecast+Forecast.PModel.aicc.arw.3$Forecast)/3
Forecast.mse.PModel.aicc.arw.avg <- mean((tc.arw[44:48]- Forecast.PModel.aicc.arw.avg)^2)

cbind(LOOCV.PModel.aicc.arw$delta[1],LOOCV.mse.PModel.aicc.arw.1$delta[1],LOOCV.mse.PModel.aicc.arw.2$delta[1],LOOCV.mse.PModel.aicc.arw.3$delta[1], LOOCV.mse.PModel.aicc.arw.avg )
cbind(Forecast.PModel.aicc.arw$delta,Forecast.PModel.aicc.arw.1$delta,Forecast.PModel.aicc.arw.2$delta,Forecast.PModel.aicc.arw.3$delta,Forecast.mse.PModel.aicc.arw.avg)

PModel.aicc.arw.simple.1
PModel.aicc.arw.simple.2
PModel.aicc.arw.simple.3

##For ARW, the model averaging improve slightly

##ARE---------------------------------------
distribution <- poisson
#AICc Linear
y <- tc.are
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.are[, which.min(PModel.aicc.are.table)]
SelectedPModel.aicc.are.second <- SelectedPModel.aicc.are[,PModel.aicc.are.table != PModel.aicc.are.table[which.min(PModel.aicc.are.table)] ]
PModel.aicc.are.second.table <- apply(SelectedPModel.aicc.are.second,2,AICcselect)
M2 <- SelectedPModel.aicc.are.second[, which.min(PModel.aicc.are.second.table)]
SelectedPModel.aicc.are.third <- SelectedPModel.aicc.are.second[,PModel.aicc.are.second.table != PModel.aicc.are.second.table[which.min(PModel.aicc.are.second.table)] ]
PModel.aicc.are.third.table <- apply(SelectedPModel.aicc.are.third,2,AICcselect)
M3 <- SelectedPModel.aicc.are.third[, which.min(PModel.aicc.are.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.are.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.are.1 <- cv.glm1(z.1, PModel.aicc.are.simple.1)
LOOCV.mae.PModel.aicc.are.1 <- cv.glm1(z.1, PModel.aicc.are.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.are.1 <- recentHindcast(z.1, PModel.aicc.are.simple.1)

PModel.aicc.are.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.are.2 <- cv.glm1(z.2, PModel.aicc.are.simple.2)
LOOCV.mae.PModel.aicc.are.2 <- cv.glm1(z.2, PModel.aicc.are.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.are.2 <- recentHindcast(z.2, PModel.aicc.are.simple.2)

PModel.aicc.are.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.are.3 <- cv.glm1(z.3, PModel.aicc.are.simple.3)
LOOCV.mae.PModel.aicc.are.3 <- cv.glm1(z.3, PModel.aicc.are.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.are.3 <- recentHindcast(z.3, PModel.aicc.are.simple.3)

LOOCV.PModel.aicc.are.avg <- (LOOCV.mse.PModel.aicc.are.1$Hindcast + LOOCV.mse.PModel.aicc.are.2$Hindcast + LOOCV.mse.PModel.aicc.are.3$Hindcast)/3
LOOCV.mse.PModel.aicc.are.avg <- mean((tc.are - LOOCV.PModel.aicc.are.avg )^2)
Forecast.PModel.aicc.are.avg <- (Forecast.PModel.aicc.are.1$Forecast+Forecast.PModel.aicc.are.2$Forecast+Forecast.PModel.aicc.are.3$Forecast)/3
Forecast.mse.PModel.aicc.are.avg <- mean((tc.are[44:48]- Forecast.PModel.aicc.are.avg)^2)

cbind(LOOCV.PModel.aicc.are$delta[1],LOOCV.mse.PModel.aicc.are.1$delta[1],LOOCV.mse.PModel.aicc.are.2$delta[1],LOOCV.mse.PModel.aicc.are.3$delta[1], LOOCV.mse.PModel.aicc.are.avg )
cbind(Forecast.PModel.aicc.are$delta,Forecast.PModel.aicc.are.1$delta,Forecast.PModel.aicc.are.2$delta,Forecast.PModel.aicc.are.3$delta,Forecast.mse.PModel.aicc.are.avg)

PModel.aicc.are.simple.1
PModel.aicc.are.simple.2
PModel.aicc.are.simple.3

#The best models contains six variables and the second best contains only three

####Linear vs Poisson----------------------------------------------

cbind(LOOCV.LModel.aicc.ar$delta[1],LOOCV.mse.LModel.aicc.ar.1$delta[1],LOOCV.mse.LModel.aicc.ar.2$delta[1],LOOCV.mse.LModel.aicc.ar.3$delta[1], LOOCV.mse.LModel.aicc.ar.avg )
cbind(Forecast.LModel.aicc.ar$delta,Forecast.LModel.aicc.ar.1$delta,Forecast.LModel.aicc.ar.2$delta,Forecast.LModel.aicc.ar.3$delta,Forecast.mse.LModel.aicc.ar.avg)

cbind(LOOCV.PModel.aicc.ar$delta[1],LOOCV.mse.PModel.aicc.ar.1$delta[1],LOOCV.mse.PModel.aicc.ar.2$delta[1],LOOCV.mse.PModel.aicc.ar.3$delta[1], LOOCV.mse.PModel.aicc.ar.avg )
cbind(Forecast.PModel.aicc.ar$delta,Forecast.PModel.aicc.ar.1$delta,Forecast.PModel.aicc.ar.2$delta,Forecast.PModel.aicc.ar.3$delta,Forecast.mse.PModel.aicc.ar.avg)

cbind(LOOCV.LModel.aicc.arw$delta[1],LOOCV.mse.LModel.aicc.arw.1$delta[1],LOOCV.mse.LModel.aicc.arw.2$delta[1],LOOCV.mse.LModel.aicc.arw.3$delta[1], LOOCV.mse.LModel.aicc.arw.avg )
cbind(Forecast.LModel.aicc.arw$delta,Forecast.LModel.aicc.arw.1$delta,Forecast.LModel.aicc.arw.2$delta,Forecast.LModel.aicc.arw.3$delta,Forecast.mse.LModel.aicc.arw.avg)

cbind(LOOCV.PModel.aicc.arw$delta[1],LOOCV.mse.PModel.aicc.arw.1$delta[1],LOOCV.mse.PModel.aicc.arw.2$delta[1],LOOCV.mse.PModel.aicc.arw.3$delta[1], LOOCV.mse.PModel.aicc.arw.avg )
cbind(Forecast.PModel.aicc.arw$delta,Forecast.PModel.aicc.arw.1$delta,Forecast.PModel.aicc.arw.2$delta,Forecast.PModel.aicc.arw.3$delta,Forecast.mse.PModel.aicc.arw.avg)

cbind(LOOCV.LModel.aicc.are$delta[1],LOOCV.mse.LModel.aicc.are.1$delta[1],LOOCV.mse.LModel.aicc.are.2$delta[1],LOOCV.mse.LModel.aicc.are.3$delta[1], LOOCV.mse.LModel.aicc.are.avg )
cbind(Forecast.LModel.aicc.are$delta,Forecast.LModel.aicc.are.1$delta,Forecast.LModel.aicc.are.2$delta,Forecast.LModel.aicc.are.3$delta,Forecast.mse.LModel.aicc.are.avg)

cbind(LOOCV.PModel.aicc.are$delta[1],LOOCV.mse.PModel.aicc.are.1$delta[1],LOOCV.mse.PModel.aicc.are.2$delta[1],LOOCV.mse.PModel.aicc.are.3$delta[1], LOOCV.mse.PModel.aicc.are.avg )
cbind(Forecast.PModel.aicc.are$delta,Forecast.PModel.aicc.are.1$delta,Forecast.PModel.aicc.are.2$delta,Forecast.PModel.aicc.are.3$delta,Forecast.mse.PModel.aicc.are.avg)




##SP------------------------------ 
distribution <- gaussian
#AICc Linear
y <- tc.sp
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.sp[, which.min(LModel.aicc.sp.table)]
SelectedLModel.aicc.sp.second <- SelectedLModel.aicc.sp[,LModel.aicc.sp.table != LModel.aicc.sp.table[which.min(LModel.aicc.sp.table)] ]
LModel.aicc.sp.second.table <- apply(SelectedLModel.aicc.sp.second,2,AICcselect)
M2 <- SelectedLModel.aicc.sp.second[, which.min(LModel.aicc.sp.second.table)]
SelectedLModel.aicc.sp.third <- SelectedLModel.aicc.sp.second[,LModel.aicc.sp.second.table != LModel.aicc.sp.second.table[which.min(LModel.aicc.sp.second.table)] ]
LModel.aicc.sp.third.table <- apply(SelectedLModel.aicc.sp.third,2,AICcselect)
M3 <- SelectedLModel.aicc.sp.third[, which.min(LModel.aicc.sp.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.sp.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.sp.1 <- cv.glm1(z.1, LModel.aicc.sp.simple.1)
LOOCV.mae.LModel.aicc.sp.1 <- cv.glm1(z.1, LModel.aicc.sp.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.sp.1 <- recentHindcast(z.1, LModel.aicc.sp.simple.1)

LModel.aicc.sp.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.sp.2 <- cv.glm1(z.2, LModel.aicc.sp.simple.2)
LOOCV.mae.LModel.aicc.sp.2 <- cv.glm1(z.2, LModel.aicc.sp.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.sp.2 <- recentHindcast(z.2, LModel.aicc.sp.simple.2)

LModel.aicc.sp.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.sp.3 <- cv.glm1(z.3, LModel.aicc.sp.simple.3)
LOOCV.mae.LModel.aicc.sp.3 <- cv.glm1(z.3, LModel.aicc.sp.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.sp.3 <- recentHindcast(z.3, LModel.aicc.sp.simple.3)

LOOCV.LModel.aicc.sp.avg <- (LOOCV.mse.LModel.aicc.sp.1$Hindcast + LOOCV.mse.LModel.aicc.sp.2$Hindcast + LOOCV.mse.LModel.aicc.sp.3$Hindcast)/3
LOOCV.mse.LModel.aicc.sp.avg <- mean((tc.sp - LOOCV.LModel.aicc.sp.avg )^2)
Forecast.LModel.aicc.sp.avg <- (Forecast.LModel.aicc.sp.1$Forecast+Forecast.LModel.aicc.sp.2$Forecast+Forecast.LModel.aicc.sp.3$Forecast)/3
Forecast.mse.LModel.aicc.sp.avg <- mean((tc.sp[44:48]- Forecast.LModel.aicc.sp.avg)^2)

cbind(LOOCV.LModel.aicc.sp$delta[1],LOOCV.mse.LModel.aicc.sp.1$delta[1],LOOCV.mse.LModel.aicc.sp.2$delta[1],LOOCV.mse.LModel.aicc.sp.3$delta[1], LOOCV.mse.LModel.aicc.sp.avg )
cbind(Forecast.LModel.aicc.sp$delta,Forecast.LModel.aicc.sp.1$delta,Forecast.LModel.aicc.sp.2$delta,Forecast.LModel.aicc.sp.3$delta,Forecast.mse.LModel.aicc.sp.avg)

LModel.aicc.sp.simple.1
LModel.aicc.sp.simple.2
LModel.aicc.sp.simple.3

#The model averging does not improve the performance: 1, the AIC for three models are almost same, 2, they have similar number of variables.

##SPW---------------------------------------
distribution <- gaussian
#AICc Linear
y <- tc.spw
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.spw[, which.min(LModel.aicc.spw.table)]
SelectedLModel.aicc.spw.second <- SelectedLModel.aicc.spw[,LModel.aicc.spw.table != LModel.aicc.spw.table[which.min(LModel.aicc.spw.table)] ]
LModel.aicc.spw.second.table <- apply(SelectedLModel.aicc.spw.second,2,AICcselect)
M2 <- SelectedLModel.aicc.spw.second[, which.min(LModel.aicc.spw.second.table)]
SelectedLModel.aicc.spw.third <- SelectedLModel.aicc.spw.second[,LModel.aicc.spw.second.table != LModel.aicc.spw.second.table[which.min(LModel.aicc.spw.second.table)] ]
LModel.aicc.spw.third.table <- apply(SelectedLModel.aicc.spw.third,2,AICcselect)
M3 <- SelectedLModel.aicc.spw.third[, which.min(LModel.aicc.spw.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.spw.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.spw.1 <- cv.glm1(z.1, LModel.aicc.spw.simple.1)
LOOCV.mae.LModel.aicc.spw.1 <- cv.glm1(z.1, LModel.aicc.spw.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spw.1 <- recentHindcast(z.1, LModel.aicc.spw.simple.1)

LModel.aicc.spw.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.spw.2 <- cv.glm1(z.2, LModel.aicc.spw.simple.2)
LOOCV.mae.LModel.aicc.spw.2 <- cv.glm1(z.2, LModel.aicc.spw.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spw.2 <- recentHindcast(z.2, LModel.aicc.spw.simple.2)

LModel.aicc.spw.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.spw.3 <- cv.glm1(z.3, LModel.aicc.spw.simple.3)
LOOCV.mae.LModel.aicc.spw.3 <- cv.glm1(z.3, LModel.aicc.spw.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spw.3 <- recentHindcast(z.3, LModel.aicc.spw.simple.3)

LOOCV.LModel.aicc.spw.avg <- (LOOCV.mse.LModel.aicc.spw.1$Hindcast + LOOCV.mse.LModel.aicc.spw.2$Hindcast + LOOCV.mse.LModel.aicc.spw.3$Hindcast)/3
LOOCV.mse.LModel.aicc.spw.avg <- mean((tc.spw - LOOCV.LModel.aicc.spw.avg )^2)
Forecast.LModel.aicc.spw.avg <- (Forecast.LModel.aicc.spw.1$Forecast+Forecast.LModel.aicc.spw.2$Forecast+Forecast.LModel.aicc.spw.3$Forecast)/3
Forecast.mse.LModel.aicc.spw.avg <- mean((tc.spw[44:48]- Forecast.LModel.aicc.spw.avg)^2)

cbind(LOOCV.LModel.aicc.spw$delta[1],LOOCV.mse.LModel.aicc.spw.1$delta[1],LOOCV.mse.LModel.aicc.spw.2$delta[1],LOOCV.mse.LModel.aicc.spw.3$delta[1], LOOCV.mse.LModel.aicc.spw.avg )
cbind(Forecast.LModel.aicc.spw$delta,Forecast.LModel.aicc.spw.1$delta,Forecast.LModel.aicc.spw.2$delta,Forecast.LModel.aicc.spw.3$delta,Forecast.mse.LModel.aicc.spw.avg)

LModel.aicc.spw.simple.1
LModel.aicc.spw.simple.2
LModel.aicc.spw.simple.3

##For ARW, the model averaging improve slightly

##SPE------------------------------------
distribution <- gaussian
#AICc Linear
y <- tc.spe
z <- cbind(y, x.r)

M1 <- SelectedLModel.aicc.spe[, which.min(LModel.aicc.spe.table)]
SelectedLModel.aicc.spe.second <- SelectedLModel.aicc.spe[,LModel.aicc.spe.table != LModel.aicc.spe.table[which.min(LModel.aicc.spe.table)] ]
LModel.aicc.spe.second.table <- apply(SelectedLModel.aicc.spe.second,2,AICcselect)
M2 <- SelectedLModel.aicc.spe.second[, which.min(LModel.aicc.spe.second.table)]
SelectedLModel.aicc.spe.third <- SelectedLModel.aicc.spe.second[,LModel.aicc.spe.second.table != LModel.aicc.spe.second.table[which.min(LModel.aicc.spe.second.table)] ]
LModel.aicc.spe.third.table <- apply(SelectedLModel.aicc.spe.third,2,AICcselect)
M3 <- SelectedLModel.aicc.spe.third[, which.min(LModel.aicc.spe.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

LModel.aicc.spe.simple.1 <- glm(formula = y ~ ., family = gaussian, data = z.1)
LOOCV.mse.LModel.aicc.spe.1 <- cv.glm1(z.1, LModel.aicc.spe.simple.1)
LOOCV.mae.LModel.aicc.spe.1 <- cv.glm1(z.1, LModel.aicc.spe.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spe.1 <- recentHindcast(z.1, LModel.aicc.spe.simple.1)

LModel.aicc.spe.simple.2 <- glm(formula = y ~ ., family = gaussian, data = z.2)
LOOCV.mse.LModel.aicc.spe.2 <- cv.glm1(z.2, LModel.aicc.spe.simple.2)
LOOCV.mae.LModel.aicc.spe.2 <- cv.glm1(z.2, LModel.aicc.spe.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spe.2 <- recentHindcast(z.2, LModel.aicc.spe.simple.2)

LModel.aicc.spe.simple.3 <- glm(formula = y ~ ., family = gaussian, data = z.3)
LOOCV.mse.LModel.aicc.spe.3 <- cv.glm1(z.3, LModel.aicc.spe.simple.3)
LOOCV.mae.LModel.aicc.spe.3 <- cv.glm1(z.3, LModel.aicc.spe.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.spe.3 <- recentHindcast(z.3, LModel.aicc.spe.simple.3)

LOOCV.LModel.aicc.spe.avg <- (LOOCV.mse.LModel.aicc.spe.1$Hindcast + LOOCV.mse.LModel.aicc.spe.2$Hindcast + LOOCV.mse.LModel.aicc.spe.3$Hindcast)/3
LOOCV.mse.LModel.aicc.spe.avg <- mean((tc.spe - LOOCV.LModel.aicc.spe.avg )^2)
Forecast.LModel.aicc.spe.avg <- (Forecast.LModel.aicc.spe.1$Forecast+Forecast.LModel.aicc.spe.2$Forecast+Forecast.LModel.aicc.spe.3$Forecast)/3
Forecast.mse.LModel.aicc.spe.avg <- mean((tc.spe[44:48]- Forecast.LModel.aicc.spe.avg)^2)

cbind(LOOCV.LModel.aicc.spe$delta[1],LOOCV.mse.LModel.aicc.spe.1$delta[1],LOOCV.mse.LModel.aicc.spe.2$delta[1],LOOCV.mse.LModel.aicc.spe.3$delta[1], LOOCV.mse.LModel.aicc.spe.avg )
cbind(Forecast.LModel.aicc.spe$delta,Forecast.LModel.aicc.spe.1$delta,Forecast.LModel.aicc.spe.2$delta,Forecast.LModel.aicc.spe.3$delta,Forecast.mse.LModel.aicc.spe.avg)

LModel.aicc.spe.simple.1
LModel.aicc.spe.simple.2
LModel.aicc.spe.simple.3

#The best models contains six variables and the second best contains only three




##SP------------------------------ 
distribution <- poisson
#AICc Linear
y <- tc.sp
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.sp[, which.min(PModel.aicc.sp.table)]
SelectedPModel.aicc.sp.second <- SelectedPModel.aicc.sp[,PModel.aicc.sp.table != PModel.aicc.sp.table[which.min(PModel.aicc.sp.table)] ]
PModel.aicc.sp.second.table <- apply(SelectedPModel.aicc.sp.second,2,AICcselect)
M2 <- SelectedPModel.aicc.sp.second[, which.min(PModel.aicc.sp.second.table)]
SelectedPModel.aicc.sp.third <- SelectedPModel.aicc.sp.second[,PModel.aicc.sp.second.table != PModel.aicc.sp.second.table[which.min(PModel.aicc.sp.second.table)] ]
PModel.aicc.sp.third.table <- apply(SelectedPModel.aicc.sp.third,2,AICcselect)
M3 <- SelectedPModel.aicc.sp.third[, which.min(PModel.aicc.sp.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.sp.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.sp.1 <- cv.glm1(z.1, PModel.aicc.sp.simple.1)
LOOCV.mae.PModel.aicc.sp.1 <- cv.glm1(z.1, PModel.aicc.sp.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.sp.1 <- recentHindcast(z.1, PModel.aicc.sp.simple.1)

PModel.aicc.sp.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.sp.2 <- cv.glm1(z.2, PModel.aicc.sp.simple.2)
LOOCV.mae.PModel.aicc.sp.2 <- cv.glm1(z.2, PModel.aicc.sp.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.sp.2 <- recentHindcast(z.2, PModel.aicc.sp.simple.2)

PModel.aicc.sp.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.sp.3 <- cv.glm1(z.3, PModel.aicc.sp.simple.3)
LOOCV.mae.PModel.aicc.sp.3 <- cv.glm1(z.3, PModel.aicc.sp.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.sp.3 <- recentHindcast(z.3, PModel.aicc.sp.simple.3)

LOOCV.PModel.aicc.sp.avg <- (LOOCV.mse.PModel.aicc.sp.1$Hindcast + LOOCV.mse.PModel.aicc.sp.2$Hindcast + LOOCV.mse.PModel.aicc.sp.3$Hindcast)/3
LOOCV.mse.PModel.aicc.sp.avg <- mean((tc.sp - LOOCV.PModel.aicc.sp.avg )^2)
Forecast.PModel.aicc.sp.avg <- (Forecast.PModel.aicc.sp.1$Forecast+Forecast.PModel.aicc.sp.2$Forecast+Forecast.PModel.aicc.sp.3$Forecast)/3
Forecast.mse.PModel.aicc.sp.avg <- mean((tc.sp[44:48]- Forecast.PModel.aicc.sp.avg)^2)

cbind(LOOCV.PModel.aicc.sp$delta[1],LOOCV.mse.PModel.aicc.sp.1$delta[1],LOOCV.mse.PModel.aicc.sp.2$delta[1],LOOCV.mse.PModel.aicc.sp.3$delta[1], LOOCV.mse.PModel.aicc.sp.avg )
cbind(Forecast.PModel.aicc.sp$delta,Forecast.PModel.aicc.sp.1$delta,Forecast.PModel.aicc.sp.2$delta,Forecast.PModel.aicc.sp.3$delta,Forecast.mse.PModel.aicc.sp.avg)

PModel.aicc.sp.simple.1
PModel.aicc.sp.simple.2
PModel.aicc.sp.simple.3

#The model averging does not improve the performance: 1, the AIC for three models are almost same, 2, they have similar number of variables.

##SPW--------------------------------------
distribution <- poisson
#AICc Linear
y <- tc.spw
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.spw[, which.min(PModel.aicc.spw.table)]
SelectedPModel.aicc.spw.second <- SelectedPModel.aicc.spw[,PModel.aicc.spw.table != PModel.aicc.spw.table[which.min(PModel.aicc.spw.table)] ]
PModel.aicc.spw.second.table <- apply(SelectedPModel.aicc.spw.second,2,AICcselect)
M2 <- SelectedPModel.aicc.spw.second[, which.min(PModel.aicc.spw.second.table)]
SelectedPModel.aicc.spw.third <- SelectedPModel.aicc.spw.second[,PModel.aicc.spw.second.table != PModel.aicc.spw.second.table[which.min(PModel.aicc.spw.second.table)] ]
PModel.aicc.spw.third.table <- apply(SelectedPModel.aicc.spw.third,2,AICcselect)
M3 <- SelectedPModel.aicc.spw.third[, which.min(PModel.aicc.spw.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.spw.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.spw.1 <- cv.glm1(z.1, PModel.aicc.spw.simple.1)
LOOCV.mae.PModel.aicc.spw.1 <- cv.glm1(z.1, PModel.aicc.spw.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spw.1 <- recentHindcast(z.1, PModel.aicc.spw.simple.1)

PModel.aicc.spw.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.spw.2 <- cv.glm1(z.2, PModel.aicc.spw.simple.2)
LOOCV.mae.PModel.aicc.spw.2 <- cv.glm1(z.2, PModel.aicc.spw.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spw.2 <- recentHindcast(z.2, PModel.aicc.spw.simple.2)

PModel.aicc.spw.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.spw.3 <- cv.glm1(z.3, PModel.aicc.spw.simple.3)
LOOCV.mae.PModel.aicc.spw.3 <- cv.glm1(z.3, PModel.aicc.spw.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spw.3 <- recentHindcast(z.3, PModel.aicc.spw.simple.3)

LOOCV.PModel.aicc.spw.avg <- (LOOCV.mse.PModel.aicc.spw.1$Hindcast + LOOCV.mse.PModel.aicc.spw.2$Hindcast + LOOCV.mse.PModel.aicc.spw.3$Hindcast)/3
LOOCV.mse.PModel.aicc.spw.avg <- mean((tc.spw - LOOCV.PModel.aicc.spw.avg )^2)
Forecast.PModel.aicc.spw.avg <- (Forecast.PModel.aicc.spw.1$Forecast+Forecast.PModel.aicc.spw.2$Forecast+Forecast.PModel.aicc.spw.3$Forecast)/3
Forecast.mse.PModel.aicc.spw.avg <- mean((tc.spw[44:48]- Forecast.PModel.aicc.spw.avg)^2)

cbind(LOOCV.PModel.aicc.spw$delta[1],LOOCV.mse.PModel.aicc.spw.1$delta[1],LOOCV.mse.PModel.aicc.spw.2$delta[1],LOOCV.mse.PModel.aicc.spw.3$delta[1], LOOCV.mse.PModel.aicc.spw.avg )
cbind(Forecast.PModel.aicc.spw$delta,Forecast.PModel.aicc.spw.1$delta,Forecast.PModel.aicc.spw.2$delta,Forecast.PModel.aicc.spw.3$delta,Forecast.mse.PModel.aicc.spw.avg)

PModel.aicc.spw.simple.1
PModel.aicc.spw.simple.2
PModel.aicc.spw.simple.3

##For ARW, the model averaging improve slightly

##SPE------------------------------------
distribution <- poisson
#AICc Linear
y <- tc.spe
z <- cbind(y, x.r)

M1 <- SelectedPModel.aicc.spe[, which.min(PModel.aicc.spe.table)]
SelectedPModel.aicc.spe.second <- SelectedPModel.aicc.spe[,PModel.aicc.spe.table != PModel.aicc.spe.table[which.min(PModel.aicc.spe.table)] ]
PModel.aicc.spe.second.table <- apply(SelectedPModel.aicc.spe.second,2,AICcselect)
M2 <- SelectedPModel.aicc.spe.second[, which.min(PModel.aicc.spe.second.table)]
SelectedPModel.aicc.spe.third <- SelectedPModel.aicc.spe.second[,PModel.aicc.spe.second.table != PModel.aicc.spe.second.table[which.min(PModel.aicc.spe.second.table)] ]
PModel.aicc.spe.third.table <- apply(SelectedPModel.aicc.spe.third,2,AICcselect)
M3 <- SelectedPModel.aicc.spe.third[, which.min(PModel.aicc.spe.third.table)]

z.1 <- z[,M1]
z.2 <- z[,M2]
z.3 <- z[,M3]

PModel.aicc.spe.simple.1 <- glm(formula = y ~ ., family = poisson, data = z.1)
LOOCV.mse.PModel.aicc.spe.1 <- cv.glm1(z.1, PModel.aicc.spe.simple.1)
LOOCV.mae.PModel.aicc.spe.1 <- cv.glm1(z.1, PModel.aicc.spe.simple.1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spe.1 <- recentHindcast(z.1, PModel.aicc.spe.simple.1)

PModel.aicc.spe.simple.2 <- glm(formula = y ~ ., family = poisson, data = z.2)
LOOCV.mse.PModel.aicc.spe.2 <- cv.glm1(z.2, PModel.aicc.spe.simple.2)
LOOCV.mae.PModel.aicc.spe.2 <- cv.glm1(z.2, PModel.aicc.spe.simple.2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spe.2 <- recentHindcast(z.2, PModel.aicc.spe.simple.2)

PModel.aicc.spe.simple.3 <- glm(formula = y ~ ., family = poisson, data = z.3)
LOOCV.mse.PModel.aicc.spe.3 <- cv.glm1(z.3, PModel.aicc.spe.simple.3)
LOOCV.mae.PModel.aicc.spe.3 <- cv.glm1(z.3, PModel.aicc.spe.simple.3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.spe.3 <- recentHindcast(z.3, PModel.aicc.spe.simple.3)

LOOCV.PModel.aicc.spe.avg <- (LOOCV.mse.PModel.aicc.spe.1$Hindcast + LOOCV.mse.PModel.aicc.spe.2$Hindcast + LOOCV.mse.PModel.aicc.spe.3$Hindcast)/3
LOOCV.mse.PModel.aicc.spe.avg <- mean((tc.spe - LOOCV.PModel.aicc.spe.avg )^2)
Forecast.PModel.aicc.spe.avg <- (Forecast.PModel.aicc.spe.1$Forecast+Forecast.PModel.aicc.spe.2$Forecast+Forecast.PModel.aicc.spe.3$Forecast)/3
Forecast.mse.PModel.aicc.spe.avg <- mean((tc.spe[44:48]- Forecast.PModel.aicc.spe.avg)^2)

cbind(LOOCV.PModel.aicc.spe$delta[1],LOOCV.mse.PModel.aicc.spe.1$delta[1],LOOCV.mse.PModel.aicc.spe.2$delta[1],LOOCV.mse.PModel.aicc.spe.3$delta[1], LOOCV.mse.PModel.aicc.spe.avg )
cbind(Forecast.PModel.aicc.spe$delta,Forecast.PModel.aicc.spe.1$delta,Forecast.PModel.aicc.spe.2$delta,Forecast.PModel.aicc.spe.3$delta,Forecast.mse.PModel.aicc.spe.avg)

PModel.aicc.spe.simple.1
PModel.aicc.spe.simple.2
PModel.aicc.spe.simple.3

#The best models contains six variables and the second best contains only three

####Linear vs Poisson----------------------------------------------

cbind(LOOCV.LModel.aicc.sp$delta[1],LOOCV.mse.LModel.aicc.sp.1$delta[1],LOOCV.mse.LModel.aicc.sp.2$delta[1],LOOCV.mse.LModel.aicc.sp.3$delta[1], LOOCV.mse.LModel.aicc.sp.avg )
cbind(Forecast.LModel.aicc.sp$delta,Forecast.LModel.aicc.sp.1$delta,Forecast.LModel.aicc.sp.2$delta,Forecast.LModel.aicc.sp.3$delta,Forecast.mse.LModel.aicc.sp.avg)

cbind(LOOCV.PModel.aicc.sp$delta[1],LOOCV.mse.PModel.aicc.sp.1$delta[1],LOOCV.mse.PModel.aicc.sp.2$delta[1],LOOCV.mse.PModel.aicc.sp.3$delta[1], LOOCV.mse.PModel.aicc.sp.avg )
cbind(Forecast.PModel.aicc.sp$delta,Forecast.PModel.aicc.sp.1$delta,Forecast.PModel.aicc.sp.2$delta,Forecast.PModel.aicc.sp.3$delta,Forecast.mse.PModel.aicc.sp.avg)

cbind(LOOCV.LModel.aicc.spw$delta[1],LOOCV.mse.LModel.aicc.spw.1$delta[1],LOOCV.mse.LModel.aicc.spw.2$delta[1],LOOCV.mse.LModel.aicc.spw.3$delta[1], LOOCV.mse.LModel.aicc.spw.avg )
cbind(Forecast.LModel.aicc.spw$delta,Forecast.LModel.aicc.spw.1$delta,Forecast.LModel.aicc.spw.2$delta,Forecast.LModel.aicc.spw.3$delta,Forecast.mse.LModel.aicc.spw.avg)

cbind(LOOCV.PModel.aicc.spw$delta[1],LOOCV.mse.PModel.aicc.spw.1$delta[1],LOOCV.mse.PModel.aicc.spw.2$delta[1],LOOCV.mse.PModel.aicc.spw.3$delta[1], LOOCV.mse.PModel.aicc.spw.avg )
cbind(Forecast.PModel.aicc.spw$delta,Forecast.PModel.aicc.spw.1$delta,Forecast.PModel.aicc.spw.2$delta,Forecast.PModel.aicc.spw.3$delta,Forecast.mse.PModel.aicc.spw.avg)

cbind(LOOCV.LModel.aicc.spe$delta[1],LOOCV.mse.LModel.aicc.spe.1$delta[1],LOOCV.mse.LModel.aicc.spe.2$delta[1],LOOCV.mse.LModel.aicc.spe.3$delta[1], LOOCV.mse.LModel.aicc.spe.avg )
cbind(Forecast.LModel.aicc.spe$delta,Forecast.LModel.aicc.spe.1$delta,Forecast.LModel.aicc.spe.2$delta,Forecast.LModel.aicc.spe.3$delta,Forecast.mse.LModel.aicc.spe.avg)

cbind(LOOCV.PModel.aicc.spe$delta[1],LOOCV.mse.PModel.aicc.spe.1$delta[1],LOOCV.mse.PModel.aicc.spe.2$delta[1],LOOCV.mse.PModel.aicc.spe.3$delta[1], LOOCV.mse.PModel.aicc.spe.avg )
cbind(Forecast.PModel.aicc.spe$delta,Forecast.PModel.aicc.spe.1$delta,Forecast.PModel.aicc.spe.2$delta,Forecast.PModel.aicc.spe.3$delta,Forecast.mse.PModel.aicc.spe.avg)

###Best model----------------------------------------------------

##AR

summary(LModel.aicc.ar.simple)
summary(PModel.aicc.ar.simple)

###Fitted plot
plot(Year, tc.ar, ylab = "Fitted TC counts", main = "AR", ylim = c(0,25), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.ar.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.ar.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.ar, ylab = "Hindcast TC counts", main = "AR", ylim = c(0,25), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.ar$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.ar$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.ar[44:48], ylab = "TC counts", main = "AR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.ar$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.ar$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

##ARW

summary(LModel.aicc.arw.simple)
summary(PModel.aicc.arw.simple)

###Fitted plot
plot(Year, tc.arw, ylab = "Fitted TC counts", main = "ARW", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.arw.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.arw.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.arw, ylab = "Hindcast TC counts", main = "ARW", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.arw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.arw$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.arw[44:48], ylab = "TC counts", main = "ARW", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.arw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.arw$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

##ARE

summary(LModel.aicc.are.simple)
summary(PModel.aicc.are.simple)

###Fitted plot
plot(Year, tc.are, ylab = "Fitted TC counts", main = "ARE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.are.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.are.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.are, ylab = "Hindcast TC counts", main = "ARE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.are$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.are$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.are[44:48], ylab = "TC counts", main = "ARE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.are$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.are$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

##SP

summary(LModel.aicc.sp.simple)
summary(PModel.aicc.sp.simple)

###Fitted plot
plot(Year, tc.sp, ylab = "Fitted TC counts", main = "SP", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.sp.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.sp.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.sp, ylab = "Hindcast TC counts", main = "SP", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.sp$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.sp$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.sp[44:48], ylab = "TC counts", main = "SP", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.sp$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.sp$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

##SPW

summary(LModel.aicc.spw.simple)
summary(PModel.aicc.spw.simple)

###Fitted plot
plot(Year, tc.spw, ylab = "Fitted TC counts", main = "SPW", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.spw.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.spw.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.spw, ylab = "Hindcast TC counts", main = "SPW", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.spw$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.spw$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.spw[44:48], ylab = "TC counts", main = "SPW", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.spw$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.spw$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

##SPE

summary(LModel.aicc.spe.simple)
summary(PModel.aicc.spe.simple)

###Fitted plot
plot(Year, tc.spe, ylab = "Fitted TC counts", main = "SPE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue" )
lg <- c("Counts", "Linear",  "Poisson")
lines(Year, fitted(LModel.aicc.spe.simple), col = c("red"), lwd = 2)
lines(Year, fitted(PModel.aicc.spe.simple), col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Hindcast plot
plot(Year, tc.spe, ylab = "Hindcast TC counts", main = "SPE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson" )
lines(Year, LOOCV.mse.LModel.aicc.spe$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.spe$Hindcast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[44:48], tc.spe[44:48], ylab = "TC counts", main = "SPE", ylim = c(0,15), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue")
lg <- c("Counts", "Linear",  "Poisson")
lines(Year[44:48], Forecast.LModel.aicc.spe$Forecast, col = c("red"), lwd = 2)
lines(Year[44:48], Forecast.PModel.aicc.spe$Forecast, col = c("blue"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###Table for AR

AR.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.ar$delta[1])/sd(tc.ar), sqrt(LOOCV.mse.LModel.aicc.ar$delta[1])/sd(tc.ar), sqrt(LOOCV.mse.PModel.aicc.ar$delta[1])/sd(tc.ar))
AR.table.mae <- cbind(LOOCV.mae.NULModel.ar$delta[1], LOOCV.mae.LModel.aicc.ar$delta[1],LOOCV.mae.PModel.aicc.ar$delta[1])
AR.table <- rbind(AR.table.nmse, AR.table.mae)
colnames(AR.table) <- c("Climatic", "AR/L", "AR/P")
rownames(AR.table) <- c("nMSE","MAE")
AR.table

ARW.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.arw$delta[1])/sd(tc.arw), sqrt(LOOCV.mse.LModel.aicc.arw$delta[1])/sd(tc.arw), sqrt(LOOCV.mse.PModel.aicc.arw$delta[1])/sd(tc.arw))
ARW.table.mae <- cbind(LOOCV.mae.NULModel.arw$delta[1],LOOCV.mae.LModel.aicc.arw$delta[1],LOOCV.mae.PModel.aicc.arw$delta[1])
ARW.table <- rbind(ARW.table.nmse, ARW.table.mae)
colnames(ARW.table) <- c("Climatic","ARW/L", "ARW/P")
rownames(ARW.table) <- c("nMSE","MAE")
ARW.table

ARE.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.are$delta[1])/sd(tc.are), sqrt(LOOCV.mse.LModel.aicc.are$delta[1])/sd(tc.are), sqrt(LOOCV.mse.PModel.aicc.are$delta[1])/sd(tc.are))
ARE.table.mae <- cbind(LOOCV.mae.NULModel.are$delta[1],LOOCV.mae.LModel.aicc.are$delta[1],LOOCV.mae.PModel.aicc.are$delta[1])
ARE.table <- rbind(ARE.table.nmse, ARE.table.mae)
colnames(ARE.table) <- c("Climatic","ARE/L", "ARE/P")
rownames(ARE.table) <- c("nMSE","MAE")
ARE.table

### Table for SPO

SP.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.sp$delta[1])/sd(tc.sp), sqrt(LOOCV.mse.LModel.aicc.sp$delta[1])/sd(tc.sp), sqrt(LOOCV.mse.PModel.aicc.sp$delta[1])/sd(tc.sp))
SP.table.mae <- cbind(LOOCV.mae.NULModel.sp$delta[1], LOOCV.mae.LModel.aicc.sp$delta[1],LOOCV.mae.PModel.aicc.sp$delta[1])
SP.table <- rbind(SP.table.nmse, SP.table.mae)
colnames(SP.table) <- c("Climatic", "SP/L", "SP/P")
rownames(SP.table) <- c("nMSE","MAE")
SP.table

SPW.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.spw$delta[1])/sd(tc.spw), sqrt(LOOCV.mse.LModel.aicc.spw$delta[1])/sd(tc.spw), sqrt(LOOCV.mse.PModel.aicc.spw$delta[1])/sd(tc.spw))
SPW.table.mae <- cbind(LOOCV.mae.NULModel.spw$delta[1],LOOCV.mae.LModel.aicc.spw$delta[1],LOOCV.mae.PModel.aicc.spw$delta[1])
SPW.table <- rbind(SPW.table.nmse, SPW.table.mae)
colnames(SPW.table) <- c("Climatic","SPW/L", "SPW/P")
rownames(SPW.table) <- c("nMSE","MAE")
SPW.table

SPE.table.nmse <- cbind(sqrt(LOOCV.mse.NULModel.spe$delta[1])/sd(tc.spe), sqrt(LOOCV.mse.LModel.aicc.spe$delta[1])/sd(tc.spe), sqrt(LOOCV.mse.PModel.aicc.spe$delta[1])/sd(tc.spe))
SPE.table.mae <- cbind(LOOCV.mae.NULModel.spe$delta[1],LOOCV.mae.LModel.aicc.spe$delta[1],LOOCV.mae.PModel.aicc.spe$delta[1])
SPE.table <- rbind(SPE.table.nmse, SPE.table.mae)
colnames(SPE.table) <- c("Climatic","SPE/L", "SPE/P")
rownames(SPE.table) <- c("nMSE","MAE")
ARE.table

