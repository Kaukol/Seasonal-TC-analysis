#Library
library(dr)
library(glmnet)
library(SISIR)
library(plsRglm)
library(AICcmodavg)
library(ggplot2)
library(doParallel)

###Set seed-------------------------------------------------
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


###Search program-----------------------------------------------------
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
#x.r <- cbind(x1[,3:14], x2[,3:14], x3[,3:14])
x.r <- cbind(x1[,3:14], x2[,3:14])

##Remark on SWIO from 1983-now
tc.swio <- tc.counts[14:48,]$SWIO
tc.swio1 <- tc.counts[14:48,]$SWIO.1
tc.swio2 <- tc.counts[14:48,]$SWIO.2
tc.swio3 <- tc.counts[14:48,]$SWIO.3
x.r.swio <- x.r[14:48,]

Year <-tc.counts$Year[14:48]
n <- dim(x.r.swio)[1]
p <- dim(x.r.swio)[2]

###Random start model -----------------------------------------------------
Multiplier <- 5
Seq_length <- 300
Cut_length <- 100
Loop <- 200
threshold_num <- 0.8

###Linear Model -------------------------------------------------------------------------------------------
distribution <- gaussian

###AIC-------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.aic.swio <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.swio <- SelectedLModel.aic.swio[,colSums(SelectedLModel.aic.swio)!=1]

LModel.aic.swio.table <- apply(SelectedLModel.aic.swio,2,AICselect)
table(LModel.aic.swio.table)

LModel.aic.swio <- glm(y~., data = z[,SelectedLModel.aic.swio[,which.min(LModel.aic.swio.table)]], family = distribution)
summary(LModel.aic.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.aic.swio1 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.swio1 <- SelectedLModel.aic.swio1[,colSums(SelectedLModel.aic.swio1)!=1]

LModel.aic.swio1.table <- apply(SelectedLModel.aic.swio1,2,AICselect)
table(LModel.aic.swio1.table)

LModel.aic.swio1 <- glm(y~., data = z[,SelectedLModel.aic.swio1[,which.min(LModel.aic.swio1.table)]], family = distribution)
summary(LModel.aic.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.aic.swio2 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.swio2 <- SelectedLModel.aic.swio2[,colSums(SelectedLModel.aic.swio2)!=1]

LModel.aic.swio2.table <- apply(SelectedLModel.aic.swio2,2,AICselect)
table(LModel.aic.swio2.table)

LModel.aic.swio2 <- glm(y~., data = z[,SelectedLModel.aic.swio2[,which.min(LModel.aic.swio2.table)]], family = distribution)
summary(LModel.aic.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.aic.swio3 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedLModel.aic.swio3 <- SelectedLModel.aic.swio3[,colSums(SelectedLModel.aic.swio3)!=1]

LModel.aic.swio3.table <- apply(SelectedLModel.aic.swio3,2,AICselect)
table(LModel.aic.swio3.table)

LModel.aic.swio3 <- glm(y~., data = z[,SelectedLModel.aic.swio3[,which.min(LModel.aic.swio3.table)]], family = distribution)
summary(LModel.aic.swio3)
stopCluster(cl)

###BIC -------------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.bic.swio <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.swio <- SelectedLModel.bic.swio[,colSums(SelectedLModel.bic.swio)!=1]

LModel.bic.swio.table <- apply(SelectedLModel.bic.swio,2,AICselect)
table(LModel.bic.swio.table)

LModel.bic.swio <- glm(y~., data = z[,SelectedLModel.bic.swio[,which.min(LModel.bic.swio.table)]], family = distribution)
summary(LModel.bic.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.bic.swio1 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.swio1 <- SelectedLModel.bic.swio1[,colSums(SelectedLModel.bic.swio1)!=1]

LModel.bic.swio1.table <- apply(SelectedLModel.bic.swio1,2,AICselect)
table(LModel.bic.swio1.table)

LModel.bic.swio1 <- glm(y~., data = z[,SelectedLModel.bic.swio1[,which.min(LModel.bic.swio1.table)]], family = distribution)
summary(LModel.bic.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.bic.swio2 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.swio2 <- SelectedLModel.bic.swio2[,colSums(SelectedLModel.bic.swio2)!=1]

LModel.bic.swio2.table <- apply(SelectedLModel.bic.swio2,2,AICselect)
table(LModel.bic.swio2.table)

LModel.bic.swio2 <- glm(y~., data = z[,SelectedLModel.bic.swio2[,which.min(LModel.bic.swio2.table)]], family = distribution)
summary(LModel.bic.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedLModel.bic.swio3 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedLModel.bic.swio3 <- SelectedLModel.bic.swio3[,colSums(SelectedLModel.bic.swio3)!=1]

LModel.bic.swio3.table <- apply(SelectedLModel.bic.swio3,2,AICselect)
table(LModel.bic.swio3.table)

LModel.bic.swio3 <- glm(y~., data = z[,SelectedLModel.bic.swio3[,which.min(LModel.bic.swio3.table)]], family = distribution)
summary(LModel.bic.swio3)
stopCluster(cl)

###AICc--------------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedLModel.aicc.swio <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.swio <- SelectedLModel.aicc.swio[,colSums(SelectedLModel.aicc.swio)!=1]

LModel.aicc.swio.table <- apply(SelectedLModel.aicc.swio,2,AICselect)
table(LModel.aicc.swio.table)

LModel.aicc.swio <- glm(y~., data = z[,SelectedLModel.aicc.swio[,which.min(LModel.aicc.swio.table)]], family = distribution)
summary(LModel.aicc.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedLModel.aicc.swio1 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.swio1 <- SelectedLModel.aicc.swio1[,colSums(SelectedLModel.aicc.swio1)!=1]

LModel.aicc.swio1.table <- apply(SelectedLModel.aicc.swio1,2,AICselect)
table(LModel.aicc.swio1.table)

LModel.aicc.swio1 <- glm(y~., data = z[,SelectedLModel.aicc.swio1[,which.min(LModel.aicc.swio1.table)]], family = distribution)
summary(LModel.aicc.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedLModel.aicc.swio2 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.swio2 <- SelectedLModel.aicc.swio2[,colSums(SelectedLModel.aicc.swio2)!=1]

LModel.aicc.swio2.table <- apply(SelectedLModel.aicc.swio2,2,AICselect)
table(LModel.aicc.swio2.table)

LModel.aicc.swio2 <- glm(y~., data = z[,SelectedLModel.aicc.swio2[,which.min(LModel.aicc.swio2.table)]], family = distribution)
summary(LModel.aicc.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedLModel.aicc.swio3 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.swio3 <- SelectedLModel.aicc.swio3[,colSums(SelectedLModel.aicc.swio3)!=1]

LModel.aicc.swio3.table <- apply(SelectedLModel.aicc.swio3,2,AICselect)
table(LModel.aicc.swio3.table)

LModel.aicc.swio3 <- glm(y~., data = z[,SelectedLModel.aicc.swio3[,which.min(LModel.aicc.swio3.table)]], family = distribution)
summary(LModel.aicc.swio3)
stopCluster(cl)


###Poisson Model------------------------------------------------------
distribution <- poisson

###AIC-------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.aic.swio <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.swio <- SelectedPModel.aic.swio[,colSums(SelectedPModel.aic.swio)!=1]

PModel.aic.swio.table <- apply(SelectedPModel.aic.swio,2,AICselect)
table(PModel.aic.swio.table)

PModel.aic.swio <- glm(y~., data = z[,SelectedPModel.aic.swio[,which.min(PModel.aic.swio.table)]], family = distribution)
summary(PModel.aic.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.aic.swio1 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.swio1 <- SelectedPModel.aic.swio1[,colSums(SelectedPModel.aic.swio1)!=1]

PModel.aic.swio1.table <- apply(SelectedPModel.aic.swio1,2,AICselect)
table(PModel.aic.swio1.table)

PModel.aic.swio1 <- glm(y~., data = z[,SelectedPModel.aic.swio1[,which.min(PModel.aic.swio1.table)]], family = distribution)
summary(PModel.aic.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.aic.swio2 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.swio2 <- SelectedPModel.aic.swio2[,colSums(SelectedPModel.aic.swio2)!=1]

PModel.aic.swio2.table <- apply(SelectedPModel.aic.swio2,2,AICselect)
table(PModel.aic.swio2.table)

PModel.aic.swio2 <- glm(y~., data = z[,SelectedPModel.aic.swio2[,which.min(PModel.aic.swio2.table)]], family = distribution)
summary(PModel.aic.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.aic.swio3 <- parSapply(cl,StartPoint,fastsearch.aic)
##Exclude null model
SelectedPModel.aic.swio3 <- SelectedPModel.aic.swio3[,colSums(SelectedPModel.aic.swio3)!=1]

PModel.aic.swio3.table <- apply(SelectedPModel.aic.swio3,2,AICselect)
table(PModel.aic.swio3.table)

PModel.aic.swio3 <- glm(y~., data = z[,SelectedPModel.aic.swio3[,which.min(PModel.aic.swio3.table)]], family = distribution)
summary(PModel.aic.swio3)
stopCluster(cl)

###BIC -------------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.bic.swio <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.swio <- SelectedPModel.bic.swio[,colSums(SelectedPModel.bic.swio)!=1]

PModel.bic.swio.table <- apply(SelectedPModel.bic.swio,2,BICselect)
table(PModel.bic.swio.table)

###NULL Model
PModel.bic.swio <- glm(y~., data = z[,SelectedPModel.bic.swio[,which.min(PModel.bic.swio.table)]], family = distribution)
summary(PModel.bic.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.bic.swio1 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.swio1 <- SelectedPModel.bic.swio1[,colSums(SelectedPModel.bic.swio1)!=1]

PModel.bic.swio1.table <- apply(SelectedPModel.bic.swio1,2,BICselect)
table(PModel.bic.swio1.table)

### NULL Model 
PModel.bic.swio1 <- glm(y~., data = z[,SelectedPModel.bic.swio1[,which.min(PModel.bic.swio1.table)]], family = distribution)
summary(PModel.bic.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.bic.swio2 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.swio2 <- SelectedPModel.bic.swio2[,colSums(SelectedPModel.bic.swio2)!=1]

PModel.bic.swio2.table <- apply(SelectedPModel.bic.swio2,2,BICselect)
table(PModel.bic.swio2.table)

PModel.bic.swio2 <- glm(y~., data = z[,SelectedPModel.bic.swio2[,which.min(PModel.bic.swio2.table)]], family = distribution)
summary(PModel.bic.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.5), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
SelectedPModel.bic.swio3 <- parSapply(cl,StartPoint,fastsearch.bic)
##Exclude null model
SelectedPModel.bic.swio3 <- SelectedPModel.bic.swio3[,colSums(SelectedPModel.bic.swio3)!=1]

PModel.bic.swio3.table <- apply(SelectedPModel.bic.swio3,2,BICselect)
table(PModel.bic.swio3.table)

PModel.bic.swio3 <- glm(y~., data = z[,SelectedPModel.bic.swio3[,which.min(PModel.bic.swio3.table)]], family = distribution)
summary(PModel.bic.swio3)
stopCluster(cl)

###AICc--------------------------------------------------------------------------------
###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
z <- as.data.frame(cbind(y,x.r.swio))
###Setting up start point
i <- 1
StartPoint <- list()
for(i in 1:Loop){
  StartPoint <- c(StartPoint, list(list(StartModelIndex = rbinom(p,1,0.2), Multiplier = Multiplier , Seq = Seq_length, Cut = Cut_length, distribution = distribution, database = z, threshold = threshold_num)))
  i <- i+1
}

###Searching
cl <- makeCluster(ncores)
clusterEvalQ(cl, library(AICcmodavg))
SelectedPModel.aicc.swio <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.swio <- SelectedPModel.aicc.swio[,colSums(SelectedPModel.aicc.swio)!=1]

PModel.aicc.swio.table <- apply(SelectedPModel.aicc.swio,2,AICcselect)
table(PModel.aicc.swio.table)

PModel.aicc.swio <- glm(y~., data = z[,SelectedPModel.aicc.swio[,which.min(PModel.aicc.swio.table)]], family = distribution)
summary(PModel.aicc.swio)
stopCluster(cl)

###SWIO1 
y <- tc.swio1
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedPModel.aicc.swio1 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.swio1 <- SelectedPModel.aicc.swio1[,colSums(SelectedPModel.aicc.swio1)!=1]

PModel.aicc.swio1.table <- apply(SelectedPModel.aicc.swio1,2,AICcselect)
table(PModel.aicc.swio1.table)

PModel.aicc.swio1 <- glm(y~., data = z[,SelectedPModel.aicc.swio1[,which.min(PModel.aicc.swio1.table)]], family = distribution)
summary(PModel.aicc.swio1)
stopCluster(cl)

###SWIO2
y <- tc.swio2
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedPModel.aicc.swio2 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.swio2 <- SelectedPModel.aicc.swio2[,colSums(SelectedPModel.aicc.swio2)!=1]

PModel.aicc.swio2.table <- apply(SelectedPModel.aicc.swio2,2,AICcselect)
table(PModel.aicc.swio2.table)

PModel.aicc.swio2 <- glm(y~., data = z[,SelectedPModel.aicc.swio2[,which.min(PModel.aicc.swio2.table)]], family = distribution)
summary(PModel.aicc.swio2)
stopCluster(cl)

###SWIO3 
y <- tc.swio3
z <- as.data.frame(cbind(y,x.r.swio))
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
SelectedPModel.aicc.swio3 <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedPModel.aicc.swio3 <- SelectedPModel.aicc.swio3[,colSums(SelectedPModel.aicc.swio3)!=1]

PModel.aicc.swio3.table <- apply(SelectedPModel.aicc.swio3,2,AICcselect)
table(PModel.aicc.swio3.table)

PModel.aicc.swio3 <- glm(y~., data = z[,SelectedPModel.aicc.swio3[,which.min(PModel.aicc.swio3.table)]], family = distribution)
summary(PModel.aicc.swio3)
stopCluster(cl)

###Summary---------------------------------------------------------
###SWIO--------------------------------------------

table(LModel.aic.swio.table)
summary(LModel.aic.swio)
table(LModel.bic.swio.table)
summary(LModel.bic.swio)
table(LModel.aicc.swio.table)
summary(LModel.aicc.swio)

table(PModel.aic.swio.table)
summary(PModel.aic.swio)
table(PModel.bic.swio.table)
summary(PModel.bic.swio)
table(PModel.aicc.swio.table)
summary(PModel.aicc.swio)

LInfo.swio <- rbind(c(AIC(LModel.aic.swio),AIC(LModel.bic.swio),AIC(LModel.aicc.swio),AIC(PModel.aic.swio),AIC(PModel.bic.swio),AIC(PModel.aicc.swio)), 
                    c(BIC(LModel.aic.swio),BIC(LModel.bic.swio),BIC(LModel.aicc.swio),BIC(PModel.aic.swio),BIC(PModel.bic.swio),BIC(PModel.aicc.swio)),
                    c(AICc(LModel.aic.swio),AICc(LModel.bic.swio),AICc(LModel.aicc.swio),AICc(PModel.aic.swio),AICc(PModel.bic.swio),AICc(PModel.aicc.swio)))

rownames(LInfo.swio) <- c("AIC","BIC","AICc")
colnames(LInfo.swio) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.swio

###Fitting plot
plot(Year, tc.swio, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.swio
z <- cbind(y, x.r.swio)
NULModel.swio <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.swio <- cv.glm1(z, NULModel.swio)
LOOCV.mae.NULModel.swio <- cv.glm1(z, NULModel.swio, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.swio <- recentHindcast(z, NULModel.swio)

z1 <- z[, SelectedLModel.aic.swio[, which.min(LModel.aic.swio.table)]]
LModel.aic.swio.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.swio <- cv.glm1(z1, LModel.aic.swio.simple)
LOOCV.mae.LModel.aic.swio <- cv.glm1(z1, LModel.aic.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.swio <- recentHindcast(z1, LModel.aic.swio.simple)

z1 <- z[, SelectedLModel.bic.swio[, which.min(LModel.bic.swio.table)]]
LModel.bic.swio.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.swio <- cv.glm1(z1, LModel.bic.swio.simple)
LOOCV.mae.LModel.bic.swio <- cv.glm1(z1, LModel.bic.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.swio <- recentHindcast(z1, LModel.bic.swio.simple)

z1 <- z[, SelectedLModel.aicc.swio[, which.min(LModel.aicc.swio.table)]]
LModel.aicc.swio.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.swio <- cv.glm1(z1, LModel.aicc.swio.simple)
LOOCV.mae.LModel.aicc.swio <- cv.glm1(z1, LModel.aicc.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.swio <- recentHindcast(z1, LModel.aicc.swio.simple)

z1 <- z[, SelectedPModel.aic.swio[, which.min(PModel.aic.swio.table)]]
PModel.aic.swio.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.swio <- cv.glm1(z1, PModel.aic.swio.simple)
LOOCV.mae.PModel.aic.swio <- cv.glm1(z1, PModel.aic.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.swio <- recentHindcast(z1, PModel.aic.swio.simple)

z1 <- z[, SelectedPModel.bic.swio[, which.min(PModel.bic.swio.table)]]
PModel.bic.swio.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.swio <- cv.glm1(z1, PModel.bic.swio.simple)
LOOCV.mae.PModel.bic.swio <- cv.glm1(z1, PModel.bic.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.swio <- recentHindcast(z1, PModel.bic.swio.simple)

z1 <- z[, SelectedPModel.aicc.swio[, which.min(PModel.aicc.swio.table)]]
PModel.aicc.swio.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.swio <- cv.glm1(z1, PModel.aicc.swio.simple)
LOOCV.mae.PModel.aicc.swio <- cv.glm1(z1, PModel.aicc.swio.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.swio <- recentHindcast(z1, PModel.aicc.swio.simple)

Hindcast.mse.swio <- cbind(LOOCV.mse.NULModel.swio$delta, LOOCV.mse.LModel.aic.swio$delta,LOOCV.mse.LModel.bic.swio$delta,LOOCV.mse.LModel.aicc.swio$delta,
                           LOOCV.mse.PModel.aic.swio$delta,LOOCV.mse.PModel.bic.swio$delta,LOOCV.mse.PModel.aicc.swio$delta)
rownames(Hindcast.mse.swio) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.swio) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.swio

Hindcast.mae.swio <- cbind(LOOCV.mae.NULModel.swio$delta, LOOCV.mae.LModel.aic.swio$delta,LOOCV.mae.LModel.bic.swio$delta,LOOCV.mae.LModel.aicc.swio$delta,
                           LOOCV.mae.PModel.aic.swio$delta,LOOCV.mae.PModel.bic.swio$delta,LOOCV.mae.PModel.aicc.swio$delta)
rownames(Hindcast.mae.swio) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.swio) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.swio

Forecast.swio <- cbind(Forecast.NULModel.swio$delta, Forecast.LModel.aic.swio$delta,Forecast.LModel.bic.swio$delta,Forecast.LModel.aicc.swio$delta,
                       Forecast.PModel.aic.swio$delta,Forecast.PModel.bic.swio$delta,Forecast.PModel.aicc.swio$delta)
rownames(Forecast.swio) <- c("MSE")
colnames(Forecast.swio) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.swio

###Forecast plot
plot(Year, tc.swio, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.swio$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.swio$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.swio$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.swio$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.swio$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.swio$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[31:35], tc.swio[31:35], ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.LModel.aic.swio$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.LModel.bic.swio$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.LModel.aicc.swio$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[31:35], tc.swio[31:35], ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.PModel.aic.swio$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.PModel.bic.swio$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.PModel.aicc.swio$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###SWIO1--------------------------------------------

PModel.bic.swio1 <- glm(tc.swio1~1, family = distribution)

table(LModel.aic.swio1.table)
summary(LModel.aic.swio1)
table(LModel.bic.swio1.table)
summary(LModel.bic.swio1)
table(LModel.aicc.swio1.table)
summary(LModel.aicc.swio1)

table(PModel.aic.swio1.table)
summary(PModel.aic.swio1)
table(PModel.bic.swio1.table)
summary(PModel.bic.swio1)
table(PModel.aicc.swio1.table)
summary(PModel.aicc.swio1)

LInfo.swio1 <- rbind(c(AIC(LModel.aic.swio1),AIC(LModel.bic.swio1),AIC(LModel.aicc.swio1),AIC(PModel.aic.swio1),AIC(PModel.bic.swio1),AIC(PModel.aicc.swio1)), 
                     c(BIC(LModel.aic.swio1),BIC(LModel.bic.swio1),BIC(LModel.aicc.swio1),BIC(PModel.aic.swio1),BIC(PModel.bic.swio1),BIC(PModel.aicc.swio1)),
                     c(AICc(LModel.aic.swio1),AICc(LModel.bic.swio1),AICc(LModel.aicc.swio1),AICc(PModel.aic.swio1),AICc(PModel.bic.swio1),AICc(PModel.aicc.swio1)))

rownames(LInfo.swio1) <- c("AIC","BIC","AICc")
colnames(LInfo.swio1) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.swio1

###Fitting plot
plot(Year, tc.swio1, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio1), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio1), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio1), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio1, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio1), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio1), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio1), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.swio1
z <- cbind(y, x.r.swio)
NULModel.swio1 <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.swio1 <- cv.glm1(z, NULModel.swio1)
LOOCV.mae.NULModel.swio1 <- cv.glm1(z, NULModel.swio1, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.swio1 <- recentHindcast(z, NULModel.swio1)

z1 <- z[, SelectedLModel.aic.swio1[, which.min(LModel.aic.swio1.table)]]
LModel.aic.swio1.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.swio1 <- cv.glm1(z1, LModel.aic.swio1.simple)
LOOCV.mae.LModel.aic.swio1 <- cv.glm1(z1, LModel.aic.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.swio1 <- recentHindcast(z1, LModel.aic.swio1.simple)

z1 <- z[, SelectedLModel.bic.swio1[, which.min(LModel.bic.swio1.table)]]
LModel.bic.swio1.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.swio1 <- cv.glm1(z1, LModel.bic.swio1.simple)
LOOCV.mae.LModel.bic.swio1 <- cv.glm1(z1, LModel.bic.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.swio1 <- recentHindcast(z1, LModel.bic.swio1.simple)

z1 <- z[, SelectedLModel.aicc.swio1[, which.min(LModel.aicc.swio1.table)]]
LModel.aicc.swio1.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.swio1 <- cv.glm1(z1, LModel.aicc.swio1.simple)
LOOCV.mae.LModel.aicc.swio1 <- cv.glm1(z1, LModel.aicc.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.swio1 <- recentHindcast(z1, LModel.aicc.swio1.simple)

z1 <- z[, SelectedPModel.aic.swio1[, which.min(PModel.aic.swio1.table)]]
PModel.aic.swio1.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.swio1 <- cv.glm1(z1, PModel.aic.swio1.simple)
LOOCV.mae.PModel.aic.swio1 <- cv.glm1(z1, PModel.aic.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.swio1 <- recentHindcast(z1, PModel.aic.swio1.simple)

#z1 <- z[, SelectedPModel.bic.swio1[, which.min(PModel.bic.swio1.table)]]
#PModel.bic.swio1.simple <- glm(formula = y ~ ., family = poisson, data = z1)
#LOOCV.mse.PModel.bic.swio1 <- cv.glm1(z1, PModel.bic.swio1.simple)
#LOOCV.mae.PModel.bic.swio1 <- cv.glm1(z1, PModel.bic.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
#Forecast.PModel.bic.swio1 <- recentHindcast(z1, PModel.bic.swio1.simple)

z1 <- z
PModel.bic.swio1.simple <- glm(formula = y ~ 1, family = poisson, data = z1)
LOOCV.mse.PModel.bic.swio1 <- cv.glm1(z1, PModel.bic.swio1.simple)
LOOCV.mae.PModel.bic.swio1 <- cv.glm1(z1, PModel.bic.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.swio1 <- recentHindcast(z1, PModel.bic.swio1.simple)

z1 <- z[, SelectedPModel.aicc.swio1[, which.min(PModel.aicc.swio1.table)]]
PModel.aicc.swio1.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.swio1 <- cv.glm1(z1, PModel.aicc.swio1.simple)
LOOCV.mae.PModel.aicc.swio1 <- cv.glm1(z1, PModel.aicc.swio1.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.swio1 <- recentHindcast(z1, PModel.aicc.swio1.simple)

Hindcast.mse.swio1 <- cbind(LOOCV.mse.NULModel.swio1$delta, LOOCV.mse.LModel.aic.swio1$delta,LOOCV.mse.LModel.bic.swio1$delta,LOOCV.mse.LModel.aicc.swio1$delta,
                            LOOCV.mse.PModel.aic.swio1$delta,LOOCV.mse.PModel.bic.swio1$delta,LOOCV.mse.PModel.aicc.swio1$delta)
rownames(Hindcast.mse.swio1) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.swio1) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.swio1

Hindcast.mae.swio1 <- cbind(LOOCV.mae.NULModel.swio1$delta, LOOCV.mae.LModel.aic.swio1$delta,LOOCV.mae.LModel.bic.swio1$delta,LOOCV.mae.LModel.aicc.swio1$delta,
                            LOOCV.mae.PModel.aic.swio1$delta,LOOCV.mae.PModel.bic.swio1$delta,LOOCV.mae.PModel.aicc.swio1$delta)
rownames(Hindcast.mae.swio1) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.swio1) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.swio1

Forecast.swio1 <- cbind(Forecast.NULModel.swio1$delta, Forecast.LModel.aic.swio1$delta,Forecast.LModel.bic.swio1$delta,Forecast.LModel.aicc.swio1$delta,
                        Forecast.PModel.aic.swio1$delta,Forecast.PModel.bic.swio1$delta,Forecast.PModel.aicc.swio1$delta)
rownames(Forecast.swio1) <- c("MSE")
colnames(Forecast.swio1) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.swio1

###Forecast plot
plot(Year, tc.swio1, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.swio1$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.swio1$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.swio1$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio1, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.swio1$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.swio1$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.swio1$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[31:35], tc.swio1[31:35], ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.LModel.aic.swio1$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.LModel.bic.swio1$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.LModel.aicc.swio1$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[31:35], tc.swio1[31:35], ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.PModel.aic.swio1$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.PModel.bic.swio1$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.PModel.aicc.swio1$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)


###SWIO2--------------------------------------------

table(LModel.aic.swio2.table)
summary(LModel.aic.swio2)
table(LModel.bic.swio2.table)
summary(LModel.bic.swio2)
table(LModel.aicc.swio2.table)
summary(LModel.aicc.swio2)

table(PModel.aic.swio2.table)
summary(PModel.aic.swio2)
table(PModel.bic.swio2.table)
summary(PModel.bic.swio2)
table(PModel.aicc.swio2.table)
summary(PModel.aicc.swio2)

LInfo.swio2 <- rbind(c(AIC(LModel.aic.swio2),AIC(LModel.bic.swio2),AIC(LModel.aicc.swio2),AIC(PModel.aic.swio2),AIC(PModel.bic.swio2),AIC(PModel.aicc.swio2)), 
                     c(BIC(LModel.aic.swio2),BIC(LModel.bic.swio2),BIC(LModel.aicc.swio2),BIC(PModel.aic.swio2),BIC(PModel.bic.swio2),BIC(PModel.aicc.swio2)),
                     c(AICc(LModel.aic.swio2),AICc(LModel.bic.swio2),AICc(LModel.aicc.swio2),AICc(PModel.aic.swio2),AICc(PModel.bic.swio2),AICc(PModel.aicc.swio2)))

rownames(LInfo.swio2) <- c("AIC","BIC","AICc")
colnames(LInfo.swio2) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.swio2

###Fitting plot
plot(Year, tc.swio2, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio2), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio2), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio2), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio2, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio2), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio2), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio2), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.swio2
z <- cbind(y, x.r.swio)
NULModel.swio2 <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.swio2 <- cv.glm1(z, NULModel.swio2)
LOOCV.mae.NULModel.swio2 <- cv.glm1(z, NULModel.swio2, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.swio2 <- recentHindcast(z, NULModel.swio2)

z1 <- z[, SelectedLModel.aic.swio2[, which.min(LModel.aic.swio2.table)]]
LModel.aic.swio2.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.swio2 <- cv.glm1(z1, LModel.aic.swio2.simple)
LOOCV.mae.LModel.aic.swio2 <- cv.glm1(z1, LModel.aic.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.swio2 <- recentHindcast(z1, LModel.aic.swio2.simple)

z1 <- z[, SelectedLModel.bic.swio2[, which.min(LModel.bic.swio2.table)]]
LModel.bic.swio2.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.swio2 <- cv.glm1(z1, LModel.bic.swio2.simple)
LOOCV.mae.LModel.bic.swio2 <- cv.glm1(z1, LModel.bic.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.swio2 <- recentHindcast(z1, LModel.bic.swio2.simple)

z1 <- z[, SelectedLModel.aicc.swio2[, which.min(LModel.aicc.swio2.table)]]
LModel.aicc.swio2.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.swio2 <- cv.glm1(z1, LModel.aicc.swio2.simple)
LOOCV.mae.LModel.aicc.swio2 <- cv.glm1(z1, LModel.aicc.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.swio2 <- recentHindcast(z1, LModel.aicc.swio2.simple)

z1 <- z[, SelectedPModel.aic.swio2[, which.min(PModel.aic.swio2.table)]]
PModel.aic.swio2.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.swio2 <- cv.glm1(z1, PModel.aic.swio2.simple)
LOOCV.mae.PModel.aic.swio2 <- cv.glm1(z1, PModel.aic.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.swio2 <- recentHindcast(z1, PModel.aic.swio2.simple)

z1 <- z[, SelectedPModel.bic.swio2[, which.min(PModel.bic.swio2.table)]]
PModel.bic.swio2.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.swio2 <- cv.glm1(z1, PModel.bic.swio2.simple)
LOOCV.mae.PModel.bic.swio2 <- cv.glm1(z1, PModel.bic.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.swio2 <- recentHindcast(z1, PModel.bic.swio2.simple)

#z1 <- z
#PModel.bic.swio2.simple <- glm(formula = y ~ 1, family = poisson, data = z1)
#LOOCV.mse.PModel.bic.swio2 <- cv.glm1(z1, PModel.bic.swio2.simple)
#LOOCV.mae.PModel.bic.swio2 <- cv.glm1(z1, PModel.bic.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
#Forecast.PModel.bic.swio2 <- recentHindcast(z1, PModel.bic.swio2.simple)

z1 <- z[, SelectedPModel.aicc.swio2[, which.min(PModel.aicc.swio2.table)]]
PModel.aicc.swio2.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.swio2 <- cv.glm1(z1, PModel.aicc.swio2.simple)
LOOCV.mae.PModel.aicc.swio2 <- cv.glm1(z1, PModel.aicc.swio2.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.swio2 <- recentHindcast(z1, PModel.aicc.swio2.simple)

Hindcast.mse.swio2 <- cbind(LOOCV.mse.NULModel.swio2$delta, LOOCV.mse.LModel.aic.swio2$delta,LOOCV.mse.LModel.bic.swio2$delta,LOOCV.mse.LModel.aicc.swio2$delta,
                            LOOCV.mse.PModel.aic.swio2$delta,LOOCV.mse.PModel.bic.swio2$delta,LOOCV.mse.PModel.aicc.swio2$delta)
rownames(Hindcast.mse.swio2) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.swio2) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.swio2

Hindcast.mae.swio2 <- cbind(LOOCV.mae.NULModel.swio2$delta, LOOCV.mae.LModel.aic.swio2$delta,LOOCV.mae.LModel.bic.swio2$delta,LOOCV.mae.LModel.aicc.swio2$delta,
                            LOOCV.mae.PModel.aic.swio2$delta,LOOCV.mae.PModel.bic.swio2$delta,LOOCV.mae.PModel.aicc.swio2$delta)
rownames(Hindcast.mae.swio2) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.swio2) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.swio2

Forecast.swio2 <- cbind(Forecast.NULModel.swio2$delta, Forecast.LModel.aic.swio2$delta,Forecast.LModel.bic.swio2$delta,Forecast.LModel.aicc.swio2$delta,
                        Forecast.PModel.aic.swio2$delta,Forecast.PModel.bic.swio2$delta,Forecast.PModel.aicc.swio2$delta)
rownames(Forecast.swio2) <- c("MSE")
colnames(Forecast.swio2) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.swio2

###Forecast plot
plot(Year, tc.swio2, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.swio2$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.swio2$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.swio2$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio2, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.swio2$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.swio2$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.swio2$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[31:35], tc.swio2[31:35], ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.LModel.aic.swio2$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.LModel.bic.swio2$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.LModel.aicc.swio2$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[31:35], tc.swio2[31:35], ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.PModel.aic.swio2$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.PModel.bic.swio2$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.PModel.aicc.swio2$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)


###SWIO3--------------------------------------------

table(LModel.aic.swio3.table)
summary(LModel.aic.swio3)
table(LModel.bic.swio3.table)
summary(LModel.bic.swio3)
table(LModel.aicc.swio3.table)
summary(LModel.aicc.swio3)

table(PModel.aic.swio3.table)
summary(PModel.aic.swio3)
table(PModel.bic.swio3.table)
summary(PModel.bic.swio3)
table(PModel.aicc.swio3.table)
summary(PModel.aicc.swio3)

LInfo.swio3 <- rbind(c(AIC(LModel.aic.swio3),AIC(LModel.bic.swio3),AIC(LModel.aicc.swio3),AIC(PModel.aic.swio3),AIC(PModel.bic.swio3),AIC(PModel.aicc.swio3)), 
                     c(BIC(LModel.aic.swio3),BIC(LModel.bic.swio3),BIC(LModel.aicc.swio3),BIC(PModel.aic.swio3),BIC(PModel.bic.swio3),BIC(PModel.aicc.swio3)),
                     c(AICc(LModel.aic.swio3),AICc(LModel.bic.swio3),AICc(LModel.aicc.swio3),AICc(PModel.aic.swio3),AICc(PModel.bic.swio3),AICc(PModel.aicc.swio3)))

rownames(LInfo.swio3) <- c("AIC","BIC","AICc")
colnames(LInfo.swio3) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
LInfo.swio3

###Fitting plot
plot(Year, tc.swio3, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio3), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio3), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio3), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio3, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio3), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio3), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio3), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

### LOCCV(Hindcasts)
y <- tc.swio3
z <- cbind(y, x.r.swio)
NULModel.swio3 <- glm(y~1, data = z, family = gaussian)
LOOCV.mse.NULModel.swio3 <- cv.glm1(z, NULModel.swio3)
LOOCV.mae.NULModel.swio3 <- cv.glm1(z, NULModel.swio3, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.NULModel.swio3 <- recentHindcast(z, NULModel.swio3)

z1 <- z[, SelectedLModel.aic.swio3[, which.min(LModel.aic.swio3.table)]]
LModel.aic.swio3.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aic.swio3 <- cv.glm1(z1, LModel.aic.swio3.simple)
LOOCV.mae.LModel.aic.swio3 <- cv.glm1(z1, LModel.aic.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aic.swio3 <- recentHindcast(z1, LModel.aic.swio3.simple)

z1 <- z[, SelectedLModel.bic.swio3[, which.min(LModel.bic.swio3.table)]]
LModel.bic.swio3.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.bic.swio3 <- cv.glm1(z1, LModel.bic.swio3.simple)
LOOCV.mae.LModel.bic.swio3 <- cv.glm1(z1, LModel.bic.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.bic.swio3 <- recentHindcast(z1, LModel.bic.swio3.simple)

z1 <- z[, SelectedLModel.aicc.swio3[, which.min(LModel.aicc.swio3.table)]]
LModel.aicc.swio3.simple <- glm(formula = y ~ ., family = gaussian, data = z1)
LOOCV.mse.LModel.aicc.swio3 <- cv.glm1(z1, LModel.aicc.swio3.simple)
LOOCV.mae.LModel.aicc.swio3 <- cv.glm1(z1, LModel.aicc.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.LModel.aicc.swio3 <- recentHindcast(z1, LModel.aicc.swio3.simple)

z1 <- z[, SelectedPModel.aic.swio3[, which.min(PModel.aic.swio3.table)]]
PModel.aic.swio3.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aic.swio3 <- cv.glm1(z1, PModel.aic.swio3.simple)
LOOCV.mae.PModel.aic.swio3 <- cv.glm1(z1, PModel.aic.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aic.swio3 <- recentHindcast(z1, PModel.aic.swio3.simple)

z1 <- z[, SelectedPModel.bic.swio3[, which.min(PModel.bic.swio3.table)]]
PModel.bic.swio3.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.bic.swio3 <- cv.glm1(z1, PModel.bic.swio3.simple)
LOOCV.mae.PModel.bic.swio3 <- cv.glm1(z1, PModel.bic.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.bic.swio3 <- recentHindcast(z1, PModel.bic.swio3.simple)

#z1 <- z
#PModel.bic.swio3.simple <- glm(formula = y ~ 1, family = poisson, data = z1)
#LOOCV.mse.PModel.bic.swio3 <- cv.glm1(z1, PModel.bic.swio3.simple)
#LOOCV.mae.PModel.bic.swio3 <- cv.glm1(z1, PModel.bic.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
#Forecast.PModel.bic.swio3 <- recentHindcast(z1, PModel.bic.swio3.simple)

z1 <- z[, SelectedPModel.aicc.swio3[, which.min(PModel.aicc.swio3.table)]]
PModel.aicc.swio3.simple <- glm(formula = y ~ ., family = poisson, data = z1)
LOOCV.mse.PModel.aicc.swio3 <- cv.glm1(z1, PModel.aicc.swio3.simple)
LOOCV.mae.PModel.aicc.swio3 <- cv.glm1(z1, PModel.aicc.swio3.simple, cost = function(y, yhat) mean(abs(y-yhat)))
Forecast.PModel.aicc.swio3 <- recentHindcast(z1, PModel.aicc.swio3.simple)

Hindcast.mse.swio3 <- cbind(LOOCV.mse.NULModel.swio3$delta, LOOCV.mse.LModel.aic.swio3$delta,LOOCV.mse.LModel.bic.swio3$delta,LOOCV.mse.LModel.aicc.swio3$delta,
                            LOOCV.mse.PModel.aic.swio3$delta,LOOCV.mse.PModel.bic.swio3$delta,LOOCV.mse.PModel.aicc.swio3$delta)
rownames(Hindcast.mse.swio3) <- c("MSE","MSE-adj")
colnames(Hindcast.mse.swio3) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mse.swio3

Hindcast.mae.swio3 <- cbind(LOOCV.mae.NULModel.swio3$delta, LOOCV.mae.LModel.aic.swio3$delta,LOOCV.mae.LModel.bic.swio3$delta,LOOCV.mae.LModel.aicc.swio3$delta,
                            LOOCV.mae.PModel.aic.swio3$delta,LOOCV.mae.PModel.bic.swio3$delta,LOOCV.mae.PModel.aicc.swio3$delta)
rownames(Hindcast.mae.swio3) <- c("MAE","MAE-adj")
colnames(Hindcast.mae.swio3) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Hindcast.mae.swio3

Forecast.swio3 <- cbind(Forecast.NULModel.swio3$delta, Forecast.LModel.aic.swio3$delta,Forecast.LModel.bic.swio3$delta,Forecast.LModel.aicc.swio3$delta,
                        Forecast.PModel.aic.swio3$delta,Forecast.PModel.bic.swio3$delta,Forecast.PModel.aicc.swio3$delta)
rownames(Forecast.swio3) <- c("MSE")
colnames(Forecast.swio3) <- c("Climatic", "Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
Forecast.swio3

###Forecast plot
plot(Year, tc.swio3, ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.LModel.aic.swio3$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.LModel.bic.swio3$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.LModel.aicc.swio3$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio3, ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, LOOCV.mse.PModel.aic.swio3$Hindcast, col = c("red"), lwd = 2)
lines(Year, LOOCV.mse.PModel.bic.swio3$Hindcast, col = c("blue"), lwd = 2)
lines(Year, LOOCV.mse.PModel.aicc.swio3$Hindcast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###5 Year forecast
plot(Year[31:35], tc.swio3[31:35], ylab = "TC counts", main = "SWIO LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.LModel.aic.swio3$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.LModel.bic.swio3$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.LModel.aicc.swio3$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year[31:35], tc.swio3[31:35], ylab = "TC counts", main = "SWIO PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year[31:35], Forecast.PModel.aic.swio3$Forecast, col = c("red"), lwd = 2)
lines(Year[31:35], Forecast.PModel.bic.swio3$Forecast, col = c("blue"), lwd = 2)
lines(Year[31:35], Forecast.PModel.aicc.swio3$Forecast, col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

