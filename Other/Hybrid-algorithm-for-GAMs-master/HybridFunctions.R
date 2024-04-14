#----------------------------------------------------------Hybrid metaheuristic functions--------------------------
ModellEpit=function(egyed, X, Y, csalad, faktorok, konkurv_strict, magok){
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
  
  stopCluster(cl)
  return(c(akaike,szignif, vifre))
}

#' 
#' @title Hibrid
#' 
#' @param genszam: An *int*, that gives the number of possible feature variables in the current dataset.
#' @param pop_meret**: An *int*, that determines the size of the populaton / harmony memory.
#' @param maxlepes**: An *int*, that determines the maximum number of iterations (or generations) the algorithm can run.
#' @param mutacio**: A *double*, determines the initial mutation (*bw*) probability.
#' @param HMCR**: A *double*, that determines the inital *HMCR* probability.
#' @param vegmutacio**: A *double*, that determines the mutation (*bw*) probability in the last generation (the last generation is determined in the *maxlepes* parameter).
#' @param vegHMCR**: A *double*, that determines the *HMCR* probability in the last generation (the last generation is determined in the *maxlepes* parameter).
#' @param konvergKrit**: An *int*, that determines the early stopping criterion. If the best solution does not change for the number of generations given here, the algorithm stops.
#' @param X**: A *dataframe* or a *named matrix* object, that contains the realized values of all the feature variables on the training set. **Important:** feature variables of *factor* type have to be represented by dummy variables in this object!
#' @param Y**: A *vector*, that contains the realized values of the target variable on the training set (in an order matching with that of the table given in the *X* parameter).
#' @param csalad**: A *string*, that gives the distribution of the target variable. A list of acceptable values can be found in the <a href="https://www.rdocumentation.org/packages/mgcv/versions/1.8-31/topics/family.mgcv" target="_blank">documentation</a> for the *family* parameter of the *bam* function in the *mgcv* package.
#' @param faktorok**: A *vector of strings*, that contains the names of the dummy variables representing *factor*s in the table given in the *X* parameter.
#' @param konkurv_strict**: An *int*, that controls whether the *concurvity* constraint should consider the pessimistic or the observed concurvity measure from the *mgcv* package. Pessimistic measure = 1; Observed measure = 2.
#' @param magok**: An *int*, that determines the number of CPU cores the algorithm can use for parallel computation of GAMs. It is advisable to use the number of available cores inus 1 here, not to overload your CPU. The version of the algorithm in the *HybridFunctions.R* file applies this parameter for the parallel computation of the parameter estimates of a single GAM. The version in *HybridFunctions_Parallelized.R* applies this parameter for the simultaneous computation of GAMs corresponding to each individual in the current harmony memory.
#' 
#' @return a list with two elements:
#'   1. **best**: a list with four elements that contains the parameters of the *best individual* in the last generation.
#'      1. A *string* describing the binary representation of the individual/solution.
#'      2. A *double*, that gives the pseudo R-squared value of the GAM corresponding to the best individual in the last generation.
#'      3. A *logical*, that describes whether the GAM corresponding to the best individual in the last generation satisfies the significance constraint, *S<sub>i</sub>*.
#'      4. A *logical*, that describes whether the GAM corresponding to the best individual in the last generation satisfies the concurvity constraint, *C<sub>i</sub>*.
#'   2. **konvergszamlalo**: An *int*, that indicates how many generations the best individual in the last generation spent in the memory/population before the algorithm stopped.
#' 
Hibrid=function(genszam,
                pop_meret,
                maxlepes,
                mutacio, 
                HMCR, 
                vegmutacio, 
                vegHMCR, 
                konvergKrit, 
                X, 
                Y, 
                csalad, 
                faktorok, 
                konkurv_strict, 
                magok){
  
  populacio<-matrix(nrow = pop_meret,ncol = 4)
  
  for (i in 1:pop_meret){
    egyed.akt=rbinom(genszam,1,1/2)
    #print(egyed.akt)
    #flush.console()
    if(sum(egyed.akt)==0){
      populacio[i,1]=toString(egyed.akt)
      populacio[i,2]=-200000
      populacio[i,3]=0
      populacio[i,4]=0
    } else{
      mod.out<-ModellEpit(egyed.akt,X,Y, csalad, faktorok, konkurv_strict, magok)
      populacio[i,1]=toString(egyed.akt)
      populacio[i,2]=mod.out[1]
      populacio[i,3]=mod.out[2]
      populacio[i,4]=mod.out[3]
    }
  }
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
    szamlalo=1
    for(j in 1:pop_meret){
      if(as.numeric(populacio[j,2])>atlagfit && as.numeric(populacio[j,3])==1 && as.numeric(populacio[j,4])==1){
        ujpop[szamlalo,1]=populacio[j,1]
        ujpop[szamlalo,2]=populacio[j,2]
        ujpop[szamlalo,3]=populacio[j,3]
        ujpop[szamlalo,4]=populacio[j,4]
        szamlalo=szamlalo+1
      }
    }
    
    startIndex=szamlalo
    
    for(j in startIndex:pop_meret){
      if(runif(1,0,1)<HMCR){
        valasztott=floor(runif(1, 1,startIndex))
        if(startIndex==1){
          ujpop[j,1]=legjobb[1]
        } else {
          ujpop[j,1]=populacio[valasztott,1]
        }
        ujegyed<-strsplit(ujpop[j,1],",")[[1]]
        for(k in 1:genszam){
          if(runif(1,0,1)<mutacio){
            if(ujegyed[k]==0){ujegyed[k]<-1}else{ujegyed[k]<-0}
          }
        }
        ujpop[j,1]=toString(ujegyed)
      } else{
        ujpop[j,1]=toString(rbinom(genszam,1,1/2))
      }
      egyed.akt=as.numeric(strsplit(ujpop[j,1],",")[[1]])
      #print(j)
      #print(egyed.akt)
      #flush.console()
      if(sum(egyed.akt)==0){
        ujpop[j,2]=-200000
        ujpop[j,3]=0
        ujpop[j,4]=0
      } else{
        mod.out<-ModellEpit(egyed.akt,X,Y, csalad, faktorok, konkurv_strict, magok)
        ujpop[j,2]=mod.out[1]
        ujpop[j,3]=mod.out[2]
        ujpop[j,4]=mod.out[3]
      }
    }
    
    populacio<-ujpop[order(ujpop[,4],ujpop[,3],ujpop[,2],decreasing=TRUE),]
    
    print(index)
    flush.console()
    
    HMCR=HMCR*HMCRszorzo
    mutacio=mutacio*mutacioszoro
  }
  return(c(legjobb, konvergszamlalo))
}#--Hibrid

#' 
#' 
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
