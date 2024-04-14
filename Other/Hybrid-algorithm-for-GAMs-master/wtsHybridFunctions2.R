require(foreach);
require(parallel);
require(doParallel);
source("wtsHybridFunctions1.R");
#' @param numAllCovars:  An *int*, gives the number of possible feature variables in the current dataset.              Was `genszam`.
#' @param pop_size:      An *int*, determines the size of the population / harmony memory.                             Was `pop_meret`
#' @param maxGens:       An *int*, determines the maximum number of iterations (or generations) the algorithm can run. Was `maxlepes`.
#' @param prbMut_init:   A *double*, determines the initial mutation (*bw*) probability.                               Was `mutacio`.
#' @param prbHMCR_init:  A *double*, determines the initial *prbHMCR_init* probability.                                Was `HMCR`.
#' @param prbMut_last:   A *double*, determines the mutation (*bw*) probability in the last generation (the last generation is determined in the *maxGens* parameter).
#' @param prbHMCR_last:  A *double*, determines the *prbHMCR_init* probability in the last generation (the last generation is determined in the *maxGens* parameter).
#' @param earlyStopCrit: An *int*, determines the early stopping criterion. If the best solution does not change for the number of generations given here, the algorithm stops.
#' @param X:             A *dataframe* or a *named matrix* object, that contains the realized values of all the feature variables on the training set. Important: feature variables of *factor* type have to be represented by dummy variables in this object!
#' @param Y:             A *vector*, that contains the realized values of the target variable on the training set (in an order matching with that of the table given in the *X* parameter).
#' @param mdlfam:        A *string*, that gives the distribution of the target variable. A list of acceptable values can be found in the <a href="https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/family.mgcv" target="_blank">documentation</a> for the *family* parameter of the *bam* function in the *mgcv* package.
#' @param facnms:        A *vector of strings*, that contains the names of the dummy variables representing *factor*s in the table given in the *X* parameter. Was `faktorok`.
#' @param concrv_opt:    An *int*, that controls whether the *concurvity* constraint should consider the pessimistic (1) or the observed (2) concurvity measure from the *mgcv* package. Was `konkurv_strict`.
#' @param ncores:        An *int*, that determines the number of CPU cores the algorithm can use for parallel computation of GAMs. Was `magoak`.
#' @details For `ncores`, it is advisable to use the number of available cores minus 1 here, not to overload your CPU. 
#' 
#' The version used here (based on *HybridFunctions_Parallelized.R*) applies this parameter for the simultaneous computation of 
#' GAMs corresponding to each individual in the current harmony memory. 
#' 
Hybrid=function(mdl, # [mgcv](bam) model object
                pop_size=20,      # An *int*, that determines the size of the population / harmony memory.                             Was `pop_meret`.
                maxGens=20,       # An *int*, that determines the maximum number of iterations (or generations) the algorithm can run. Was `maxlepes`.
                prbMut_init=0.9,  # A *double*, determines the initial mutation (*bw*) probability.                                    Was `mutacio`.
                prbHMCR_init=0.05,# A *double*, that determines the initial *prbHMCR_init* probability.                                Was `HMCR`.
                prbMut_last=0.1,  # A *double*, that determines the mutation (*bw*) probability in the last generation (the last generation is determined in the *maxGens* parameter). Was `vegmutacio`.
                prbHMCR_last=0.35,# A *double*, that determines the *prbHMCR_init* probability in the last generation (the last generation is determined in the *maxGens* parameter). Was `vegHMCR`.
                earlyStopCrit=5,  # An *int*, that determines the early stopping criterion. If the best solution does not change for the number of generations given here, the algorithm stops. Was `konvergKrit`.
                X,                # A *dataframe* or a *named matrix* object, that contains the realized values of all the feature variables on the training set. Important: feature variables of *factor* type have to be represented by dummy variables in this object!
                Y,                # A *vector*, that contains the realized values of the target variable on the training set (in an order matching with that of the table given in the *X* parameter).
                mdlfam,           # the *family* parameter for the *bam* function in the *mgcv* package.
                facnms,           # A *vector of strings*, that contains the names of the dummy variables representing *factor*s in the table given in the *X* parameter.
                concrv_opt=2,     # An *int*, that controls whether the *concurvity* constraint should consider the pessimistic or the observed concurvity measure from the *mgcv* package. Pessimistic measure = 1; Observed measure = 2.
                ncores=parallel::detectCores()-1 # An *int*, that determines the number of CPU cores the algorithm can use for parallel computation of GAMs.
                ){
  
  #--get model info
  mdlfmly   = family(mdl);                  #--model family object
  mdlvars   = gratia::model_vars(mdl);      #--model covariates
  mdlsmths  = gratia::smooths(mdl);         #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);      #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);    #--list, by smooth, of variables involved in each smooth
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  
  #--want to 
  numAllCovars = length(mdlvars);# number of total number of potential model covariates (Was `genszam`, an input)
  
  population<-matrix(nrow = pop_size,ncol = 4);#--was `populacio`
  
  cl <- makeCluster(ncores);
  registerDoParallel(cl);
  
  population <- foreach (i = 1:pop_size, .combine='rbind', .export="evalModel1") %dopar% {
    ids_sel_covars = rbinom(numAllCovars,1,1/2);#--random selection of model components (Was `egyed.akt`)
    nms_sel_covars = mdlvars[id_sel_covars];
    #print(egyed.akt)
    #flush.console()
    if(sum(egyed.akt)==0){
      c(toString(egyed.akt),-200000,0,0)
    } else{
      #--fit model defined by egyed.akt
      #----if initial k's too small, increase and try again
      mod.out<-evalModel1(egyed.akt,X,Y, mdlfam, facnms, concrv_opt, ncores);#--was ModellEpit
      c(toString(egyed.akt),
        mod.out[1],#--AIC
        mod.out[2],#--significance
        mod.out[3])#--concurvity check (FALSE=> concurvity>0.5)
    }
  }
  stopCluster(cl)
  
  population<-population[order(population[,2],decreasing=TRUE),]
  
  fitneszeddig=0
  konvergszamlalo=0
  HMCRszorzo=(prbHMCR_last/prbHMCR_init)^(1/(maxGens-1))
  mutacioszoro=(prbMut_last/prbMut_init)^(1/(maxGens-1))
  
  for(index in 1:maxGens){
    
    legjobb<-population[1,]
    
    if(as.numeric(legjobb[2])==fitneszeddig){
      konvergszamlalo=konvergszamlalo+1
      if(konvergszamlalo==earlyStopCrit){
        break
      }
    } else{
      fitneszeddig=as.numeric(legjobb[2])
      konvergszamlalo=0
    }
    
    if(sum(as.numeric(population[,3]))==0){
      atlagfit<-mean(as.numeric(population[,2]))
    } else{
      atlagfit<-sum(as.numeric(population[,2])*as.numeric(population[,3])*as.numeric(population[,4]))/sum(as.numeric(population[,3])*as.numeric(population[,4]))
    }
    ujpop<-matrix(nrow = pop_size,ncol = 4)
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    ujpop <- foreach (j = 1:pop_size, .combine='rbind') %dopar% {
      if(as.numeric(population[j,2])>atlagfit && as.numeric(population[j,3])==1 && as.numeric(population[j,4])==1){
        c(population[j,])
      } else{
        c(NA,NA,NA,NA)
      }
    }
    stopCluster(cl)
    
    ujpop <- ujpop[order(ujpop[,4],decreasing=TRUE),]
    
    startIndex=min(which(is.na(ujpop)))
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    ujpop[startIndex:pop_size,] <- foreach (j = startIndex:pop_size, .combine='rbind', .export="evalModel1") %dopar% {
      if(runif(1,0,1)<prbHMCR_init){
        valasztott=floor(runif(1, 1,startIndex))
        if(startIndex==1){
          ujegyed=legjobb[1]
        } else {
          ujegyed=population[valasztott,1]
        }
        ujegyed<-strsplit(ujegyed,",")[[1]]
        for(k in 1:numAllCovars){
          if(runif(1,0,1)<prbMut_init){
            #--flip selected covar
            if(ujegyed[k]==0){ujegyed[k]<-1}else{ujegyed[k]<-0}
          }
        }
        ujegyed=toString(ujegyed);
      } else{
        #--try new combination of all covars
        ujegyed=toString(rbinom(numAllCovars,1,1/2));
      }
      egyed.akt=as.numeric(strsplit(ujegyed,",")[[1]])
      #print(j)
      #print(egyed.akt)
      #flush.console()
      if(sum(egyed.akt)==0){
        c(toString(egyed.akt),-200000,0,0)
      } else {
        mod.out<-evalModel1(egyed.akt,X,Y, mdlfam, facnms, concrv_opt, ncores)
        c(toString(egyed.akt),mod.out[1],mod.out[2],mod.out[3])
      }
    }#--ujpop loop
    stopCluster(cl)
    
    population<-ujpop[order(ujpop[,4],ujpop[,3],ujpop[,2],decreasing=TRUE),]
    
    print(index)
    flush.console()
    
    prbHMCR_init=prbHMCR_init*HMCRszorzo
    prbMut_init=prbMut_init*mutacioszoro
  }
  return(c(legjobb, konvergszamlalo))
}#--Hibrid


