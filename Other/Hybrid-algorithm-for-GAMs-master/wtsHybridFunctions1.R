getSmoothsMarginalBasisTypes<-function(frmla){
  rsp = as.character(frmla)[2];
  fps = strsplit(as.character(frmla)[3]," + ",fixed=TRUE)[[1]];
  bss = rep(NA_character_,length(fps));
  for (i in 1:length(fps)){
    bsp = stringr::str_extract_all(fps[i],  "(?<=bs =).+?(?=,)")[[1]];
    if (length(bsp)>0) bss[i] = bsp;
  }
  return(bss);
}#--getSmoothsMarginalBasisTypes
if (FALSE){
  mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;
  bss = getSmoothsMarginalBasisTypes(formula(mdl));
}

#' 
#' @param resp - response name
#' @param smths - vector of smooths names 
#' @param smtrms - list (by smooth) of smooth terms 
#' @param bss  - vector of marginal basis types (from getSmoothMarginalBasisTypes)
#' @param k1d - max k for 1-d smooths 
#' @param k2d - max k for 2-d smooths
createModelFormula<-function(resp,smths,smtrms,bss,k1d,k2d){
  n = length(smths);
  str = paste0(resp,"~");
  for (i in 1:n){
    #--testing: i = 7;
    sm  = smths[i];
    smt = smtrms[[i]];
    bs  = bss[i];
    smfcn = strsplit(sm,"(",fixed=TRUE)[[1]][1];
    strp = paste0(smfcn[1],"(",paste0(smt,collapse=","),",bs=",bs,",k=c(");
    if (length(smt)==1) {
      strp = paste0(strp,k1d,")) ");
    } else {
      strp = paste0(strp,paste(c(k1d,k2d),collapse=","),")) ");
    }
    if (i<n) str = paste(str,strp,"+");
  }
  str = paste(str,strp);
  return(as.formula(str));
}#--createModelFormula
if (FALSE){
  mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;
  mdlvars   = gratia::model_vars(mdl); #--model covariates
  mdlsmths  = gratia::smooths(mdl);    #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl); #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl); #--list, by smooth, of variables involved in each smooth
  bss = getSmoothsMarginalBasisTypes(formula(mdl));
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  mdlfmly   = family(mdl);#--model family object
  createModelFormula(smths,smtrms,bss,k1d=10,k2d=5)
}
  
  
#----------------------------------------------------------Hybrid metaheuristic functions--------------------------
#' @param egeyed
#' @param X - dataframe or matrix of covariates
#' @param Y - response vector
#' @param mdlfam - model family
#' @param facnms - string vector containing the names of the dummy variables representing factor in X
#' @param concrv_opt - An int: indicates whether the concurvity constraint is based on the (1) pessimistic or the (2) observed measures from mgcv
#' @param ncores - An int: determines the number of CPU cores the algorithm can use for parallel computation of GAMs. 
#' 
#' @details For `ncores`, it is advisable to use the number of available cores minus 1, not to overload your CPU. 
#' The version of the algorithm here (based on *HybridFunctions_Parallelized.R*) applies this parameter 
#' for the simultaneous computation of GAMs corresponding to each individual in the current harmony memory.
#' 
evalModel1<-function(egyed,      #--selection vector for covariates TODO: replace by selected model terms
                     X,          #--dataframe or matrix of covariates
                     Y,          #--response vector
                     mdlfam,     #--model family
                     facnms=character(0),#--string vector containing the names of the dummy variables representing factors in X
                     concrv_opt=2,       #--An int: indicates whether the concurvity constraint is based on the (1) pessimistic or the (2) observed measures from mgcv
                     ncores              #--An int: determines the number of CPU cores the algorithm can use for parallel computation of GAMs.
                    ){
  nTermsFull=length(egyed); #number of (all) covariates; formerly genzsam
  
  #--determine model covariates (subset of all covariates)
  covars_all<-names(as.data.frame(X));#--names of (all) covariates; was alapnevek
  covars_sel<-character(sum(egyed));  #--vector of selected covariate names; was nevek
  tmpctr=1; #--temporary counter; was szamlalo
  for(i in 1:nTermsFull) {
    if(egyed[i]==1){#--covariate is selected
      covars_sel[tmpctr]<-covars_all[i]; #--add covariate term to selected
      tmpctr=tmpctr+1;                   #--increment counter
    }
  }
  rm(tmpctr);
  
  #--create data matrix for selected model covariates
  dfrSelCovarsData=as.data.frame(matrix(nrow=nrow(X),ncol = sum(egyed)));#--matrix of selected covars; was adatok
  tmpctr=1; #--temporary counter; was szamlalo
  for(i in 1:nTermsFull){
    if(egyed[i]==1){#--covariate is selected
      dfrSelCovarsData[,tmpctr]=X[,i];
      tmpctr=tmpctr+1;
    }
  }
  rm(tmpctr);
  colnames(dfrSelCovarsData)<-covars_sel;
  dfrSelCovarsData$target<-Y
  
  #library(parallel)
  #cl <- makeCluster(ncores)
  
  SpecStatus<-FALSE
  ctr1<-1; #--was felsokorlat
  #--set default k's: probably need to set separately for 1d/2d smooths
  kvektor<-rep(10,sum(egyed)); #--set kvektor[i] = 10 as default for 1d smooth
  while ((!SpecStatus)&(ctr1<=2)) {
    for (i in 1:length(kvektor)) {
      if (length(unique(dfrSelCovarsData[,i]))<10) {
        kvektor[i]=length(unique(dfrSelCovarsData[,i])) ;#--adjusting
      }
    }
    
    modelstring<-"target~"
    tmpctr=1;#--temporary counter; was szamlalo
    for(i in 1:nTermsFull) {
      if(egyed[i]==1){
        if (tmpctr==length(kvektor)) {
          if (covars_all[i] %in% facnms) {
            modelstring<-paste(modelstring, covars_all[i], sep="")
          } else {
            modelstring<-paste(modelstring, "s(",covars_all[i],",k=",kvektor[tmpctr],")", sep="")
          }
        } else {
          if (covars_all[i] %in% facnms) {
            modelstring<-paste(modelstring,covars_all[i],"+", sep="")
          } else {
            modelstring<-paste(modelstring,"s(",covars_all[i],",k=",kvektor[tmpctr],")+", sep="")
          }
        }
        tmpctr=tmpctr+1
      }
    }
    rm(tmpctr);
    
    library(mgcv)
    gam.mod<-bam(as.formula(modelstring),family=mdlfam,data=dfrSelCovarsData, method="fREML")
    ellen<-mgcv:::k.check(gam.mod)
    
    SpecStatus<-TRUE
    tmpctr<-1
    for (i in 1:length(kvektor)) {
      if (!(covars_all[i]%in% facnms)) {
        if (!is.null(ellen[tmpctr,1])&!is.null(ellen[tmpctr,2])&!is.null(ellen[tmpctr,4])) {
          tryCatch({
            if (((ellen[tmpctr,1]-ellen[tmpctr,2])<3)&(ellen[tmpctr,4]<0.05)) {
              SpecStatus<-FALSE
              kvektor[i]=kvektor[i]+5;#--k.check failed for this component. try increased value of k next time.
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
        tmpctr=tmpctr+1
      }
    }#--i in 1:length(kvektor)
    rm(tmpctr);
    ctr1=ctr1+1
  } #--while ((!SpecStatus)&(ctr1<=2))
  
  mdl_smry<-summary(gam.mod);#--was szumma
  
  akaike<-mdl_smry$r.sq
  
  signif<-TRUE;#--was szignif
  if (length(mdl_smry$p.pv)>0) {
    for(i in 1:length(mdl_smry$p.pv)){
      if(mdl_smry$p.pv[i]>0.05){
        signif<-FALSE;
        break;
      }
    }
  }
  if (length(mdl_smry$s.pv)>0) {
    for(i in 1:length(mdl_smry$s.pv)){
      if(mdl_smry$s.pv[i]>0.05){
        signif<-FALSE;
        break;
      }
    }#--i
  }
  
  concrv_tst<-TRUE;#--was vifre
  tryCatch(
    { concrv<-concurvity(gam.mod)[concrv_opt,];
      for(i in 1:length(concrv)){
        if(concrv[i]>0.5) {
          concrv_tst<-FALSE;
          break;
        }
      }#--i
    }, 
    warning = function(w) {}, 
    error = function(e) {
              concrv_tst<-FALSE
            }, 
    finally = {}
  )#--tryCatch
  
  #stopCluster(cl)
  return(c(akaike,signif, concrv_tst))
}#--evalModel; was ModellEpit

