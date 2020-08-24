#Library
library(dr)
library(glmnet)
library(SISIR)
library(plsRglm)
library(AICcmodavg)
library(ggplot2)
library(doParallel)

###Set seed
set.seed(1)

###Search program
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

##Remark on SWIO from 1983-now
tc.swio <- tc.counts[14:48,]$SWIO
tc.swio1 <- tc.counts[14:48,]$SWIO.1
tc.swio2 <- tc.counts[14:48,]$SWIO.2
tc.swio3 <- tc.counts[14:48,]$SWIO.3
x.r.swio <- x.r[14:48,]

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

###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
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
SelectedLModel.aic.swio <- parSapply(cl,StartPoint,fastsearch.aic)

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

LModel.aic.swio3.table <- apply(SelectedLModel.aic.swio3,2,AICselect)
table(LModel.aic.swio3.table)

LModel.aic.swio3 <- glm(y~., data = z[,SelectedLModel.aic.swio3[,which.min(LModel.aic.swio3.table)]], family = distribution)
summary(LModel.aic.swio3)
stopCluster(cl)

###BIC ----------------------------------------------------------------------------------------------------------

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

###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
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
SelectedLModel.bic.swio <- parSapply(cl,StartPoint,fastsearch.bic)

LModel.bic.swio.table <- apply(SelectedLModel.bic.swio,2,BICselect)
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

LModel.bic.swio1.table <- apply(SelectedLModel.bic.swio1,2,BICselect)
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

LModel.bic.swio2.table <- apply(SelectedLModel.bic.swio2,2,BICselect)
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

LModel.bic.swio3.table <- apply(SelectedLModel.bic.swio3,2,BICselect)
table(LModel.bic.swio3.table)

LModel.bic.swio3 <- glm(y~., data = z[,SelectedLModel.bic.swio3[,which.min(LModel.bic.swio3.table)]], family = distribution)
summary(LModel.bic.swio3)
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
SelectedLModel.aicc.spe <- parSapply(cl,StartPoint,fastsearch.aicc)
##Exclude null model
SelectedLModel.aicc.spe <- SelectedLModel.aicc.spe[,colSums(SelectedLModel.aicc.spe)!=1]

LModel.aicc.spe.table <- apply(SelectedLModel.aicc.spe,2,AICcselect)
table(LModel.aicc.spe.table)

LModel.aicc.spe <- glm(y~., data = z[,SelectedLModel.aicc.spe[,which.min(LModel.aicc.spe.table)]], family = distribution)
summary(LModel.aicc.spe)
stopCluster(cl)

###SWIO ----------------------------------------------------------------------------------
y <- tc.swio
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
SelectedLModel.aicc.swio <- parSapply(cl,StartPoint,fastsearch.aicc)

LModel.aicc.swio.table <- apply(SelectedLModel.aicc.swio,2,AICcselect)
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
SelectedLModel.aicc.swio1 <- parSapply(cl,StartPoint,fastsearch.aicc)

LModel.aicc.swio1.table <- apply(SelectedLModel.aicc.swio1,2,AICcselect)
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
SelectedLModel.aicc.swio2 <- parSapply(cl,StartPoint,fastsearch.aicc)

LModel.aicc.swio2.table <- apply(SelectedLModel.aicc.swio2,2,AICcselect)
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
SelectedLModel.aicc.swio3 <- parSapply(cl,StartPoint,fastsearch.aicc)

LModel.aicc.swio3.table <- apply(SelectedLModel.aicc.swio3,2,AICcselect)
table(LModel.aicc.swio3.table)

LModel.aicc.swio3 <- glm(y~., data = z[,SelectedLModel.aicc.swio3[,which.min(LModel.aicc.swio3.table)]], family = distribution)
summary(LModel.aicc.swio3)
stopCluster(cl)

####Poisson Model ----------------------------------------------------------------------------------------------------------

distribution <- poisson


###Summary   -----------------------------------------------------------------------------------

###AR--------------------------------------------

summary(LModel.aic.ar)
summary(LModel.bic.ar)
summary(LModel.aicc.ar)

summary(PModel.aic.ar)
summary(PModel.bic.ar)
summary(PModel.aicc.ar)

LInfo.ar <- rbind(c(AIC(LModel.aic.ar),AIC(LModel.bic.ar),AIC(LModel.aicc.ar),AIC(PModel.aic.ar),AIC(PModel.bic.ar),AIC(PModel.aicc.ar)), 
                  c(BIC(LModel.aic.ar),BIC(LModel.bic.ar),BIC(LModel.aicc.ar),BIC(PModel.aic.ar),BIC(PModel.bic.ar),BIC(PModel.aicc.ar)),
                  c(AICc(LModel.aic.ar),AICc(LModel.bic.ar),AICc(LModel.aicc.ar),AICc(PModel.aic.ar),AICc(PModel.bic.ar),AICc(PModel.aicc.ar)))

rownames(LInfo.ar) <- c("AIC","BIC","AICc")
colnames(LInfo.ar) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.ar)

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

###ARW

summary(LModel.aic.arw)
summary(LModel.bic.arw)
summary(LModel.aicc.arw)

summary(PModel.aic.arw)
summary(PModel.bic.arw)
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

###ARE

summary(LModel.aic.are)
summary(LModel.bic.are)
summary(LModel.aicc.are)

summary(PModel.aic.are)
summary(PModel.bic.are)
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

###SPW

summary(LModel.aic.spw)
summary(LModel.bic.spw)
summary(LModel.aicc.spw)

summary(PModel.aic.spw)
summary(PModel.bic.spw)
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

###SPE

summary(LModel.aic.spe)
summary(LModel.bic.spe)
summary(LModel.aicc.spe)

summary(PModel.aic.spe)
summary(PModel.bic.spe)
summary(PModel.aicc.spe)

LInfo.spe <- rbind(c(AIC(LModel.aic.spe),AIC(LModel.bic.spe),AIC(LModel.aicc.spe),AIC(PModel.aic.spe),AIC(PModel.bic.spe),AIC(PModel.aicc.spe)), 
                   c(BIC(LModel.aic.spe),BIC(LModel.bic.spe),BIC(LModel.aicc.spe),BIC(PModel.aic.spe),BIC(PModel.bic.spe),BIC(PModel.aicc.spe)),
                   c(AICc(LModel.aic.spe),AICc(LModel.bic.spe),AICc(LModel.aicc.spe),AICc(PModel.aic.spe),AICc(PModel.bic.spe),AICc(PModel.aicc.spe)))

rownames(LInfo.spe) <- c("AIC","BIC","AICc")
colnames(LInfo.spe) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.spe)

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

###SWIO ----------------------------------------------------------------------------

summary(LModel.aic.swio)
summary(LModel.bic.swio)
summary(LModel.aicc.swio)

summary(PModel.aic.swio)
summary(PModel.bic.swio)
summary(PModel.aicc.swio)

LInfo.swio <- rbind(c(AIC(LModel.aic.swio),AIC(LModel.bic.swio),AIC(LModel.aicc.swio),AIC(PModel.aic.swio),AIC(PModel.bic.swio),AIC(PModel.aicc.swio)), 
                    c(BIC(LModel.aic.swio),BIC(LModel.bic.swio),BIC(LModel.aicc.swio),BIC(PModel.aic.swio),BIC(PModel.bic.swio),BIC(PModel.aicc.swio)),
                    c(AICc(LModel.aic.swio),AICc(LModel.bic.swio),AICc(LModel.aicc.swio),AICc(PModel.aic.swio),AICc(PModel.bic.swio),AICc(PModel.aicc.swio)))

rownames(LInfo.swio) <- c("AIC","BIC","AICc")
colnames(LInfo.swio) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.swio)

plot(Year, tc.swio, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###SWIO.1

summary(LModel.aic.swio1)
summary(LModel.bic.swio1)
summary(LModel.aicc.swio1)

summary(PModel.aic.swio1)
summary(PModel.bic.swio1)
summary(PModel.aicc.swio1)

LInfo.swio1 <- rbind(c(AIC(LModel.aic.swio1),AIC(LModel.bic.swio1),AIC(LModel.aicc.swio1),AIC(PModel.aic.swio1),AIC(PModel.bic.swio1),AIC(PModel.aicc.swio1)), 
                     c(BIC(LModel.aic.swio1),BIC(LModel.bic.swio1),BIC(LModel.aicc.swio1),BIC(PModel.aic.swio1),BIC(PModel.bic.swio1),BIC(PModel.aicc.swio1)),
                     c(AICc(LModel.aic.swio1),AICc(LModel.bic.swio1),AICc(LModel.aicc.swio1),AICc(PModel.aic.swio1),AICc(PModel.bic.swio1),AICc(PModel.aicc.swio1)))

rownames(LInfo.swio1) <- c("AIC","BIC","AICc")
colnames(LInfo.swio1) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.swio1)

plot(Year, tc.swio1, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio1), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio1), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio1), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio1, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio1), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio1), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio1), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

###SWIO.2

summary(LModel.aic.swio2)
summary(LModel.bic.swio2)
summary(LModel.aicc.swio2)

summary(PModel.aic.swio2)
summary(PModel.bic.swio2)
summary(PModel.aicc.swio2)

LInfo.swio2 <- rbind(c(AIC(LModel.aic.swio2),AIC(LModel.bic.swio2),AIC(LModel.aicc.swio2),AIC(PModel.aic.swio2),AIC(PModel.bic.swio2),AIC(PModel.aicc.swio2)), 
                     c(BIC(LModel.aic.swio2),BIC(LModel.bic.swio2),BIC(LModel.aicc.swio2),BIC(PModel.aic.swio2),BIC(PModel.bic.swio2),BIC(PModel.aicc.swio2)),
                     c(AICc(LModel.aic.swio2),AICc(LModel.bic.swio2),AICc(LModel.aicc.swio2),AICc(PModel.aic.swio2),AICc(PModel.bic.swio2),AICc(PModel.aicc.swio2)))

rownames(LInfo.swio2) <- c("AIC","BIC","AICc")
colnames(LInfo.swio2) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.swio2)

plot(Year, tc.swio2, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio2), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio2), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio2), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio2, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio2), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio2), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio2), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)


###SWIO.3

summary(LModel.aic.swio3)
summary(LModel.bic.swio3)
summary(LModel.aicc.swio3)

summary(PModel.aic.swio3)
summary(PModel.bic.swio3)
summary(PModel.aicc.swio3)

LInfo.swio3 <- rbind(c(AIC(LModel.aic.swio3),AIC(LModel.bic.swio3),AIC(LModel.aicc.swio3),AIC(PModel.aic.swio3),AIC(PModel.bic.swio3),AIC(PModel.aicc.swio3)), 
                     c(BIC(LModel.aic.swio3),BIC(LModel.bic.swio3),BIC(LModel.aicc.swio3),BIC(PModel.aic.swio3),BIC(PModel.bic.swio3),BIC(PModel.aicc.swio3)),
                     c(AICc(LModel.aic.swio3),AICc(LModel.bic.swio3),AICc(LModel.aicc.swio3),AICc(PModel.aic.swio3),AICc(PModel.bic.swio3),AICc(PModel.aicc.swio3)))

rownames(LInfo.swio3) <- c("AIC","BIC","AICc")
colnames(LInfo.swio3) <- c("Linear.aic","Linear.bic","Linear.aicc","Poisson.aic","Poisson.bic","Poisson.aicc")
t(LInfo.swio3)

plot(Year, tc.swio3, ylab = "TC counts", main = "AR LR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(LModel.aic.swio3), col = c("red"), lwd = 2)
lines(Year, fitted.values(LModel.bic.swio3), col = c("blue"), lwd = 2)
lines(Year, fitted.values(LModel.aicc.swio3), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

plot(Year, tc.swio3, ylab = "TC counts", main = "AR PR", ylim = c(0,20), xlab = "Years", col ="grey", type = "h", lwd = 1)
cls <- c("grey", "red",  "blue", "black" )
lg <- c("Counts", "AIC",  "BIC", "AICc")
lines(Year, fitted.values(PModel.aic.swio3), col = c("red"), lwd = 2)
lines(Year, fitted.values(PModel.bic.swio3), col = c("blue"), lwd = 2)
lines(Year, fitted.values(PModel.aicc.swio3), col = c("black"), lwd = 2)
legend("topright", lty =1, lwd = 2, col = cls, legend = lg)

####------------------------------------------------------------------------------------------------------------------