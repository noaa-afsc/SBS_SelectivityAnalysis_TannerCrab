#----------------------------------------------------------Hybrid metaheuristic functions--------------------------
ModellEpit=function(egyed, 
                    X,        #--dataframe or matrix of covariates
                    Y,        #--response vector
                    csalad,   #--model family
                    faktorok, #--string vector containing the names of the dummy variables representing factor in X
                    konkurv_strict, # An int: indicates whether the concurvity constraint is based on the (1) pessimistic or the (2) observed measures from mgcv
                    magok     # An int: determines the number of CPU cores the algorithm can use for parallel computation of GAMs. It is advisable to use the number of available cores inus 1 here, not to overload your CPU. The version of the algorithm in the *HybridFunctions.R* file applies this parameter for the parallel computation of the parameter estimates of a single GAM. The version in *HybridFunctions_Parallelized.R* applies this parameter for the simultaneous computation of GAMs corresponding to each individual in the current harmony memory.
                    ){
  genszam=length(egyed)
  
  alapnevek<-names(as.data.frame(X))
  nevek<-character(sum(egyed))
  szamlalo=1
  for(i in 1:genszam) {
    if(egyed[i]==1){
      nevek[szamlalo]<-alapnevek[i]
      szamlalo=szamlalo+1
    }
  }
  
  adatok=as.data.frame(matrix(nrow=nrow(X),ncol = sum(egyed)))
  szamlalo=1
  for(i in 1:genszam){
    if(egyed[i]==1){
      adatok[,szamlalo]=X[,i]
      szamlalo=szamlalo+1
    }
  }
  colnames(adatok)<-nevek
  adatok$target<-Y
  
  #library(parallel)
  #cl <- makeCluster(magok)
  
  SpecStatus<-FALSE
  felsokorlat<-1
  kvektor<-rep(10,sum(egyed))
  while ((!SpecStatus)&(felsokorlat<=2)) {
    for (i in 1:length(kvektor)) {
      if (length(unique(adatok[,i]))<10) {
        kvektor[i]=length(unique(adatok[,i]))
      }
    }
    
    modelstring<-"target~"
    szamlalo=1
    for(i in 1:genszam) {
      if(egyed[i]==1){
        if (szamlalo==length(kvektor)) {
          if (alapnevek[i] %in% faktorok) {
            modelstring<-paste(modelstring, alapnevek[i], sep="")
          } else {
            modelstring<-paste(modelstring, "s(",alapnevek[i],",k=",kvektor[szamlalo],")", sep="")
          }
        } else {
          if (alapnevek[i] %in% faktorok) {
            modelstring<-paste(modelstring,alapnevek[i],"+", sep="")
          } else {
            modelstring<-paste(modelstring,"s(",alapnevek[i],",k=",kvektor[szamlalo],")+", sep="")
          }
        }
        szamlalo=szamlalo+1
      }
    }
    
    library(mgcv)
    gam.mod<-bam(as.formula(modelstring),family=csalad,data=adatok, method="fREML")
    ellen<-mgcv:::k.check(gam.mod)
    
    SpecStatus<-TRUE
    szamlalo<-1
    for (i in 1:length(kvektor)) {
      if (!(alapnevek[i]%in% faktorok)) {
        if (!is.null(ellen[szamlalo,1])&!is.null(ellen[szamlalo,2])&!is.null(ellen[szamlalo,4])) {
          tryCatch({
            if (((ellen[szamlalo,1]-ellen[szamlalo,2])<3)&(ellen[szamlalo,4]<0.05)) {
              SpecStatus<-FALSE
              kvektor[i]=kvektor[i]+5;#--??
            }
          },
          error=function(cond) {
            message(cond)
          },
          warning=function(cond) {
            message(cond)
          },
          finally={
            
          })
        }
        szamlalo=szamlalo+1
      }
    }#--i in 1:length(kvektor)
    felsokorlat=felsokorlat+1
  } #--while ((!SpecStatus)&(felsokorlat<=2))
  
  szumma<-summary(gam.mod)
  
  akaike<-szumma$r.sq
  
  szignif<-TRUE
  if (length(szumma$p.pv)>0) {
    for(i in 1:length(szumma$p.pv)){
      if(szumma$p.pv[i]>0.05){szignif<-FALSE
      break}
    }
  }
  if (length(szumma$s.pv)>0) {
    for(i in 1:length(szumma$s.pv)){
      if(szumma$s.pv[i]>0.05){szignif<-FALSE
      break}
    }
  }
  
  vifre<-TRUE
  tryCatch({
    konkurv<-concurvity(gam.mod)[konkurv_strict,]
    for(i in 1:length(konkurv)){
      if(konkurv[i]>0.5) {vifre<-FALSE
      break}
    }
  }, warning = function(w) {
    
  }, error = function(e) {
    vifre<-FALSE
  }, finally = {
    
  })
  
  #stopCluster(cl)
  return(c(akaike,szignif, vifre))
}#--ModellEpit

#' @param genszam:     An *int*, that gives the number of possible feature variables in the current dataset.
#' @param pop_meret:   An *int*, that determines the size of the population / harmony memory.
#' @param maxlepes:    An *int*, that determines the maximum number of iterations (or generations) the algorithm can run.
#' @param mutacio:     A *double*, determines the initial mutation (*bw*) probability.
#' @param HMCR:        A *double*, that determines the initial *HMCR* probability.
#' @param vegmutacio:  A *double*, that determines the mutation (*bw*) probability in the last generation (the last generation is determined in the *maxlepes* parameter).
#' @param vegHMCR:     A *double*, that determines the *HMCR* probability in the last generation (the last generation is determined in the *maxlepes* parameter).
#' @param konvergKrit: An *int*, that determines the early stopping criterion. If the best solution does not change for the number of generations given here, the algorithm stops.
#' @param X:           A *dataframe* or a *named matrix* object, that contains the realized values of all the feature variables on the training set. Important: feature variables of *factor* type have to be represented by dummy variables in this object!
#' @param Y:           A *vector*, that contains the realized values of the target variable on the training set (in an order matching with that of the table given in the *X* parameter).
#' @param csalad:      A *string*, that gives the distribution of the target variable. A list of acceptable values can be found in the <a href="https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/family.mgcv" target="_blank">documentation</a> for the *family* parameter of the *bam* function in the *mgcv* package.
#' @param faktorok:    A *vector of strings*, that contains the names of the dummy variables representing *factor*s in the table given in the *X* parameter.
#' @param konkurv_strict: An *int*, that controls whether the *concurvity* constraint should consider the pessimistic or the observed concurvity measure from the *mgcv* package. Pessimistic measure = 1; Observed measure = 2.
#' @param magok:          An *int*, that determines the number of CPU cores the algorithm can use for parallel computation of GAMs. It is advisable to use the number of available cores inus 1 here, not to overload your CPU. The version of the algorithm in the *HybridFunctions.R* file applies this parameter for the parallel computation of the parameter estimates of a single GAM. The version in *HybridFunctions_Parallelized.R* applies this parameter for the simultaneous computation of GAMs corresponding to each individual in the current harmony memory.
Hibrid=function(genszam,  #  An *int*, that gives the number of possible feature variables in the current dataset.
                pop_meret,# An *int*, that determines the size of the populaton / harmony memory.
                maxlepes, # An *int*, that determines the maximum number of iterations (or generations) the algorithm can run.
                mutacio,  # A *double*, determines the initial mutation (*bw*) probability.
                HMCR,     # A *double*, that determines the inital *HMCR* probability.
                vegmutacio,  # A *double*, that determines the mutation (*bw*) probability in the last generation (the last generation is determined in the *maxlepes* parameter).
                vegHMCR,     # A *double*, that determines the *HMCR* probability in the last generation (the last generation is determined in the *maxlepes* parameter).
                konvergKrit, # An *int*, that determines the early stopping criterion. If the best solution does not change for the number of generations given here, the algorithm stops.
                X,           # A *dataframe* or a *named matrix* object, that contains the realized values of all the feature variables on the training set. Important: feature variables of *factor* type have to be represented by dummy variables in this object!
                Y,           # A *vector*, that contains the realized values of the target variable on the training set (in an order matching with that of the table given in the *X* parameter).
                csalad,      # A *string*, that gives the distribution of the target variable. A list of acceptable values can be found in the <a href="https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/family.mgcv" target="_blank">documentation</a> for the *family* parameter of the *bam* function in the *mgcv* package.
                faktorok,       # A *vector of strings*, that contains the names of the dummy variables representing *factor*s in the table given in the *X* parameter.
                konkurv_strict, # An *int*, that controls whether the *concurvity* constraint should consider the pessimistic or the observed concurvity measure from the *mgcv* package. Pessimistic measure = 1; Observed measure = 2.
                magok           # An *int*, that determines the number of CPU cores the algorithm can use for parallel computation of GAMs. It is advisable to use the number of available cores inus 1 here, not to overload your CPU. The version of the algorithm in the *HybridFunctions.R* file applies this parameter for the parallel computation of the parameter estimates of a single GAM. The version in *HybridFunctions_Parallelized.R* applies this parameter for the simultaneous computation of GAMs corresponding to each individual in the current harmony memory.
                ){
  
  library(foreach)
  library(parallel)
  library(doParallel)
  
  populacio<-matrix(nrow = pop_meret,ncol = 4)
  
  cl <- makeCluster(magok)
  registerDoParallel(cl)
  
  populacio <- foreach (i = 1:pop_meret, .combine='rbind', .export="ModellEpit") %dopar% {
    egyed.akt=rbinom(genszam,1,1/2)#--random selection of model components
    #print(egyed.akt)
    #flush.console()
    if(sum(egyed.akt)==0){
      c(toString(egyed.akt),-200000,0,0)
    } else{
      #--fit model defined by egyed.akt
      mod.out<-ModellEpit(egyed.akt,X,Y, csalad, faktorok, konkurv_strict, magok)
      c(toString(egyed.akt),
        mod.out[1],#--AIC
        mod.out[2],#--significance
        mod.out[3])#--concurvity check (FALSE=> concurvity>0.5)
    }
  }
  stopCluster(cl)
  
  populacio<-populacio[order(populacio[,2],decreasing=TRUE),]
  
  fitneszeddig=0
  konvergszamlalo=0
  HMCRszorzo=(vegHMCR/HMCR)^(1/(maxlepes-1))
  mutacioszoro=(vegmutacio/mutacio)^(1/(maxlepes-1))
  
  for(index in 1:maxlepes){
    
    legjobb<-populacio[1,]
    
    if(as.numeric(legjobb[2])==fitneszeddig){
      konvergszamlalo=konvergszamlalo+1
      if(konvergszamlalo==konvergKrit){
        break
      }
    } else{
      fitneszeddig=as.numeric(legjobb[2])
      konvergszamlalo=0
    }
    
    if(sum(as.numeric(populacio[,3]))==0){
      atlagfit<-mean(as.numeric(populacio[,2]))
    } else{
      atlagfit<-sum(as.numeric(populacio[,2])*as.numeric(populacio[,3])*as.numeric(populacio[,4]))/sum(as.numeric(populacio[,3])*as.numeric(populacio[,4]))
    }
    ujpop<-matrix(nrow = pop_meret,ncol = 4)
    
    cl <- makeCluster(magok)
    registerDoParallel(cl)
    ujpop <- foreach (j = 1:pop_meret, .combine='rbind') %dopar% {
      if(as.numeric(populacio[j,2])>atlagfit && as.numeric(populacio[j,3])==1 && as.numeric(populacio[j,4])==1){
        c(populacio[j,])
      } else{
        c(NA,NA,NA,NA)
      }
    }
    stopCluster(cl)
    
    ujpop <- ujpop[order(ujpop[,4],decreasing=TRUE),]
    
    startIndex=min(which(is.na(ujpop)))
    
    cl <- makeCluster(magok)
    registerDoParallel(cl)
    ujpop[startIndex:pop_meret,] <- foreach (j = startIndex:pop_meret, .combine='rbind', .export="ModellEpit") %dopar% {
      if(runif(1,0,1)<HMCR){
        valasztott=floor(runif(1, 1,startIndex))
        if(startIndex==1){
          ujegyed=legjobb[1]
        } else {
          ujegyed=populacio[valasztott,1]
        }
        ujegyed<-strsplit(ujegyed,",")[[1]]
        for(k in 1:genszam){
          if(runif(1,0,1)<mutacio){
            if(ujegyed[k]==0){ujegyed[k]<-1}else{ujegyed[k]<-0}
          }
        }
        ujegyed=toString(ujegyed)
      } else{
        ujegyed=toString(rbinom(genszam,1,1/2))
      }
      egyed.akt=as.numeric(strsplit(ujegyed,",")[[1]])
      #print(j)
      #print(egyed.akt)
      #flush.console()
      if(sum(egyed.akt)==0){
        c(toString(egyed.akt),-200000,0,0)
      } else{
        mod.out<-ModellEpit(egyed.akt,X,Y, csalad, faktorok, konkurv_strict, magok)
        c(toString(egyed.akt),mod.out[1],mod.out[2],mod.out[3])
      }
    }
    stopCluster(cl)
    
    populacio<-ujpop[order(ujpop[,4],ujpop[,3],ujpop[,2],decreasing=TRUE),]
    
    print(index)
    flush.console()
    
    HMCR=HMCR*HMCRszorzo
    mutacio=mutacio*mutacioszoro
  }
  return(c(legjobb, konvergszamlalo))
}#--Hibrid

ModellEpit_B=function(egyed, X, Y, csalad, faktorok, magok){
  genszam=length(egyed)
  
  kvektor<-rep(10,sum(egyed))
  SpecStatus<-FALSE
  felsokorlat<-1
  
  alapnevek<-names(as.data.frame(X))
  nevek<-character(sum(egyed))
  szamlalo=1
  for(i in 1:genszam) {
    if(egyed[i]==1){
      nevek[szamlalo]<-alapnevek[i]
      szamlalo=szamlalo+1
    }
  }
  
  adatok=as.data.frame(matrix(nrow=nrow(X),ncol = sum(egyed)))
  szamlalo=1
  for(i in 1:genszam){
    if(egyed[i]==1){
      adatok[,szamlalo]=X[,i]
      szamlalo=szamlalo+1
    }
  }
  colnames(adatok)<-nevek
  adatok$target<-Y
  
  library(parallel)
  cl <- makeCluster(magok)
  
  while ((!SpecStatus)&(felsokorlat<=2)) {
    for (i in 1:length(kvektor)) {
      if (length(unique(adatok[,i]))<10) {
        kvektor[i]=length(unique(adatok[,i]))
      }
    }
    
    modelstring<-"target~"
    szamlalo=1
    for(i in 1:genszam) {
      if(egyed[i]==1){
        if (szamlalo==length(kvektor)) {
          if (alapnevek[i] %in% faktorok) {
            modelstring<-paste(modelstring, alapnevek[i], sep="")
          } else {
            modelstring<-paste(modelstring, "s(",alapnevek[i],",k=",kvektor[szamlalo],")", sep="")
          }
        } else {
          if (alapnevek[i] %in% faktorok) {
            modelstring<-paste(modelstring,alapnevek[i],"+", sep="")
          } else {
            modelstring<-paste(modelstring,"s(",alapnevek[i],",k=",kvektor[szamlalo],")+", sep="")
          }
        }
        szamlalo=szamlalo+1
      }
    }
    
    library(mgcv)
    gam.mod<-bam(as.formula(modelstring),family=csalad,data=adatok, method="fREML", cluster=cl)
    ellen<-mgcv:::k.check(gam.mod)
    
    SpecStatus<-TRUE
    szamlalo<-1
    for (i in 1:length(kvektor)) {
      if (!(alapnevek[i]%in% faktorok)) {
        if (!is.null(ellen[szamlalo,1])&!is.null(ellen[szamlalo,2])&!is.null(ellen[szamlalo,4])) {
          tryCatch({
            if (((ellen[szamlalo,1]-ellen[szamlalo,2])<3)&(ellen[szamlalo,4]<0.05)) {
              SpecStatus<-FALSE
              kvektor[i]=kvektor[i]+5
            }
          },
          error=function(cond) {
            message(cond)
          },
          warning=function(cond) {
            message(cond)
          },
          finally={
            
          })
        }
        szamlalo=szamlalo+1
      }
    }
    felsokorlat=felsokorlat+1
  }
  
  stopCluster(cl)
  return(gam.mod)
}