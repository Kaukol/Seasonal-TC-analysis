#Library
library(AICcmodavg)
library(doParallel)

#Sigmoid
Sigmoid <- function(x,InverseT){ return(1/(1+exp(-InverseT*x)))}

#Information criteria
InfoCriteria <- function(mod_1, mod_2, n, info = c("AIC", "BIC", "AICc")){
  if(info == "AIC"){
    return(Sigmoid(AIC(mod_2)-AIC(mod_1),n))
  }
  else if(info == "BIC"){
    return(Sigmoid(BIC(mod_2)-BIC(mod_1),n))
  }
  else if(info == "AICc"){
    return(Sigmoid(AICc(mod_2)-AICc(mod_1),n))
  }
  else {return("Wrong criteria")}
}


#Gibbs sampler
GibbsSampler_Search <- function(Y, X, StartModel, Seq_Length, InverseT, info = c("AIC", "BIC", "AICc"), family = c("gaussian","poisson", "binomial")){
  
  #data
  y <- Y
  x <- X
  z <- as.data.frame(cbind(y,x))
  
  #Start
  n <- InverseT
  StartIndex <- StartModel
  StartModel <- c(1, StartIndex)
  ModelSeletionMatrix <- StartModel
  
  #Loop
  j <- 0
  while(j < Seq_Length){
    for(i in 1:(ncol(z) - 1)){
      TempIndex <- StartIndex
      
      #print(StartIndex)
      #print(TempIndex)
      
      TempDual <- TempIndex
      TempDual[i] <- 1 - TempIndex[i]
      
      #print(TempIndex)
      #print(TempDual)
      
      TempCheck <- c(1,TempIndex) == 1
      TempCheckDual <- c(1, TempDual) == 1
      
      if(sum(!(TempIndex == 0)) == 0){
        Model_One <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_One <- glm(y ~ .,family = family, data = z[,TempCheck])
      }
      
      if(sum(!(TempDual == 0)) == 0){
        Model_Two <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_Two <- glm(y ~ .,family = family, data = z[,TempCheckDual])
      }
      
      #Random jump
      BinomProb <- InfoCriteria(Model_One, Model_Two, n, info = info)
      
      #print(c(BIC(Model_One), BIC(Model_Two)))
      #print(BinomProb)
      
      TossUp <- rbinom(1,1,BinomProb)
      if(TossUp == 1){
        TempIndex[i] <- TempIndex[i]
      }
      else{
        TempIndex[i] <- 1 - TempIndex[i]
      }
      StartIndex <- TempIndex
      #print(StartIndex)
    }
    ModelSeletionMatrix <- rbind(ModelSeletionMatrix, c(1,StartIndex))
    j = j + 1
  }
  return(ModelSeletionMatrix)
}

#Sigmoid
Sigmoid_MH <- function(x,InverseT){ 
  t1 <- 1/exp(-InverseT*x)
  return(min(c(1,t1)))
}

#Information criteria
InfoCriteria_MH <- function(mod_1, mod_2, n, info = c("AIC", "BIC", "AICc")){
  if(info == "AIC"){
    return(Sigmoid_MH(-AIC(mod_2)+ AIC(mod_1),n))
  }
  else if(info == "BIC"){
    return(Sigmoid_MH(-BIC(mod_2)+ BIC(mod_1),n))
  }
  else if(info == "AICc"){
    return(Sigmoid_MH(-AICc(mod_2)+ AICc(mod_1),n))
  }
  else {return("Wrong criteria")}
}

#Gibbs sampler
GibbsSampler_Search_MH <- function(Y, X, StartModel, Seq_Length, InverseT, info = c("AIC", "BIC", "AICc"), family = c("gaussian","poisson", "binomial")){
  
  #data
  y <- Y
  x <- X
  z <- as.data.frame(cbind(y,x))
  
  #Start
  n <- InverseT
  StartIndex <- StartModel
  StartModel <- c(1, StartIndex)
  ModelSeletionMatrix <- StartModel
  
  #Loop
  j <- 0
  while(j < Seq_Length){
    for(i in 1:(ncol(z) - 1)){
      TempIndex <- StartIndex
      TempDual <- TempIndex
      TempDual[i] <- 1 - TempIndex[i]
      
      TempCheck <- c(1,TempIndex) == 1
      TempCheckDual <- c(1, TempDual) == 1
      
      if(sum(!(TempIndex == 0)) == 0){
        Model_One <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_One <- glm(y ~ .,family = family, data = z[,TempCheck])
      }
      
      if(sum(!(TempDual == 0)) == 0){
        Model_Two <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_Two <- glm(y ~ .,family = family, data = z[,TempCheckDual])
      }
      
      #Random jump
      BinomProb <- InfoCriteria_MH(Model_One, Model_Two, n, info = info)
      
      #print(c(BIC(Model_One), BIC(Model_Two)))
      #print(BinomProb)
      
      TossUp <- rbinom(1,1,BinomProb)
      if(TossUp == 0){
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

#Model selection 
BurnSeq <- function(ModelSeletionMatrix, Cut_Length){
  CutModelMatrix <- ModelSeletionMatrix[(nrow(ModelSeletionMatrix) - Cut_Length + 1):nrow(ModelSeletionMatrix),]
  return(CutModelMatrix)
}


####Model selection by
InfoCritForModel <- function(SelectedModel, info = c("AIC", "BIC", "AICc"), family = c("gaussian","poisson", "binomial")){
  if(info == "AIC"){
    return(AIC(glm(y~., data = z[,SelectedModel], family = family)))
  }
  else if(info == "BIC"){
    return(BIC(glm(y~., data = z[,SelectedModel], family = family)))
  }
  else if(info == "AICc"){
    return(AICc(glm(y~., data = z[,SelectedModel], family = family)))
  }
  else {return("Wrong criteria")}
}

#Gibbs sampler simualted annealing
GibbsSampler_Search_SA <- function(Y, X, StartModel, Seq_Length, InverseT = c("Log", "Geo"), info = c("AIC", "BIC", "AICc"), family = c("gaussian","poisson", "binomial")){
  
  #data
  y <- Y
  x <- X
  z <- as.data.frame(cbind(y,x))
  
  #Start
  #n <- InverseT
  StartIndex <- StartModel
  StartModel <- c(1, StartIndex)
  ModelSeletionMatrix <- StartModel
  
  #Loop
  j <- 0
  while(j < Seq_Length){
    n <- log(j+1)/2
    for(i in 1:(ncol(z) - 1)){
      TempIndex <- StartIndex
      
      #print(StartIndex)
      #print(TempIndex)
      
      TempDual <- TempIndex
      TempDual[i] <- 1 - TempIndex[i]
      
      #print(TempIndex)
      #print(TempDual)
      
      TempCheck <- c(1,TempIndex) == 1
      TempCheckDual <- c(1, TempDual) == 1
      
      if(sum(!(TempIndex == 0)) == 0){
        Model_One <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_One <- glm(y ~ .,family = family, data = z[,TempCheck])
      }
      
      if(sum(!(TempDual == 0)) == 0){
        Model_Two <- glm(y ~ 1,family = family) #Null model
      }else{
        Model_Two <- glm(y ~ .,family = family, data = z[,TempCheckDual])
      }
      
      #Random jump
      BinomProb <- InfoCriteria(Model_One, Model_Two, n, info = info)
      
      #print(c(BIC(Model_One), BIC(Model_Two)))
      #print(BinomProb)
      
      TossUp <- rbinom(1,1,BinomProb)
      if(TossUp == 1){
        TempIndex[i] <- TempIndex[i]
      }
      else{
        TempIndex[i] <- 1 - TempIndex[i]
      }
      StartIndex <- TempIndex
      #print(StartIndex)
    }
    ModelSeletionMatrix <- rbind(ModelSeletionMatrix, c(1,StartIndex))
    j = j + 1
  }
  return(ModelSeletionMatrix)
}
