# 
require(gratia);
require(mgcv);
require(tibble);
require(foreach);
require(doParallel);
require(parallel);


dirPrj = rstudioapi::getActiveProject();
#source(file.path(dirPrj,"Other/Hybrid-algorithm-for-GAMs-master","wtsHybridFunctions1.R"));
#source(file.path(dirPrj,"Other/Hybrid-algorithm-for-GAMs-master","HybridFunctions.R"));

#---------------------------------------------
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

#---------------------------------------------
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
  
#---------------------------------------------
getCombs<-function(mdlsmths){
  smbase = mdlsmths[1]; # required (base) terms
  smxtra = mdlsmths[2:length(mdlsmths)];# potential (extra) terms to evaluate
  n_xtra = length(smxtra);# number of potential (extra)  terms to evaluate
  cmbs = tibble::tibble(x=1); names(cmbs)[ncol(cmbs)] = smbase;
  for (i in 1:n_xtra) {
    cmbs = tidyr::expand_grid(cmbs,tibble::tibble(x=0:1))
    names(cmbs)[ncol(cmbs)] = smxtra[i];
  }
  ncmbs = nrow(cmbs);
  cat("Total number of combinations:",ncmbs,"\n");
  return(cmbs);
}

#---------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param i - index of "best" model from [evalAllModels] 
#' 
#' @return fitted gam model representing "best" combination of terms from the input model `mdl`.
#' 
evalBestModel<-function(mdl,
                        i){
  resp = as.character(formula(mdl))[2];     #--model response
  mdldata   = mdl$model;                    #--data for model fit
  mdlfmly   = family(mdl);                  #--model family object
  offset_   = model.offset(mdl);  #--model offsets, if any
  weights_  = model.weights(mdl); #--model weights, if any
  mdlvars   = gratia::model_vars(mdl);      #--model covariates
  mdlsmths  = gratia::smooths(mdl);         #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);      #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);    #--list, by smooth, of variables involved in each smooth
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  mdlbss    = getSmoothsMarginalBasisTypes(formula(mdl));#--basis types for smooths
  
  #--set default values for offset_, weights_ as necessary
  if (is.null(offset_)) 
    offset_ = rep(0,nrow(mdldata));
  if (is.null(weights_)) 
    weights_ = rep(1,nrow(mdldata));
  
  #--set up models to run
  cmbs = getCombs(mdlsmths);
  
  idx = as.logical(as.vector(t(cmbs[i,])));
  if (debug) cat(paste0("\n\nEvaluating best model (",i,"): ",paste(mdlsmths[idx],collapse=" + "),"\n"));
  frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],10,8);
  best   = mgcv::gam(formula=frmla,
                     family=mdlfmly,
                     data=mdldata,
                     offset=offset_,
                     weights=weights_,
                     method="REML");
return(best);
}

#---------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' 
evalAllModels<-function(mdl,
                        concrv_opt=2,
                        doParallel=FALSE,
                        ncores=parallel::detectCores()-1,
                        max_ncmbs=NULL,
                        icmbs_=NULL,
                        debug=FALSE){
  resp = as.character(formula(mdl))[2];     #--model response
  mdldata   = mdl$model;                    #--data for model fit
  mdlfmly   = family(mdl);                  #--model family object
  mdlvars   = gratia::model_vars(mdl);      #--model covariates
  mdlsmths  = gratia::smooths(mdl);         #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);      #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);    #--list, by smooth, of variables involved in each smooth
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  mdlbss    = getSmoothsMarginalBasisTypes(formula(mdl));#--basis types for smooths
  offset_   = model.offset(mdl);  #--model offsets, if any
  weights_  = model.weights(mdl); #--model weights, if any
  
  #--set default values for offset_, weights_ as necessary
  if (is.null(offset_)) 
    offset_ = rep(0,nrow(mdldata));
  if (is.null(weights_)) 
    weights_ = rep(1,nrow(mdldata));
  
  #--set up models to run
  cmbs = getCombs(mdlsmths);
  ncmbs = nrow(cmbs);
  if (is.null(max_ncmbs)) max_ncmbs = ncmbs;
  cat("Evaluating",max_ncmbs,"combinations.\n");
  icmbs = 1:max_ncmbs;
  if (!is.null(icmbs_)) icmbs = icmbs_;

  #--set up parallel processing
  if (doParallel){
    cl <- makeCluster(ncores);
    registerDoParallel(cl)
  }
  
  if (debug) writeLines(c(""), "log.txt");
  
  lstAll<-foreach (i = icmbs,
                   .inorder=debug,
                   .export=c("createModelFormula"),
                   .packages=c("mgcv","tibble")) %dopar% {
    if (debug) sink("log.txt", append=TRUE);
    idx = as.logical(as.vector(t(cmbs[i,])));
    if (debug) cat(paste0("\n\nEvaluating model ",i,": ",paste(mdlsmths[idx],collapse=" + "),"\n"));
    frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],10,8);
    #cat("offset_:",offset_,"\n");
    #cat("weights_:",weights_,"\n");
    mdlres = mgcv::gam(formula=frmla,
                       family=mdlfmly,
                       data=mdldata,
                       offset=offset_,
                       weights=weights_,
                       method="REML");
    kchk = k.check(mdlres);#--not really doing anything with k.check result here
    if (debug) {cat("\nk.check:\n"); print(kchk);}
    mdl_smry<-summary(mdlres);
    if (is.null(mdl_smry)){
      if (debug) cat("\nNULL: model summary failed in model",i,"\n");
      rsqr    = -1;
      aic     = Inf;
      signifp = FALSE;#--default
      sigidxp = -1;
      sigvalp = -1;
      signifs = FALSE;#--default
      sigidxs = -1;
      sigvals = -1;
    } else {
      if (debug) {cat("\nSummary:\n"); print(mdl_smry);}
      rsqr      = mdl_smry$r.sq;
      aic      = AIC(mdlres);
      #--check for non-significant terms
      signifp = TRUE;#--default
      sigidxp = -1;
      sigvalp = -1;
      if (length(mdl_smry$p.pv)>0) {
        for(j in 1:length(mdl_smry$p.pv)){
          if(mdl_smry$p.pv[j]>0.05){
            signifp = FALSE;
            sigidxp = j;
            sigvalp = mdl_smry$p.pv[j];
            break;
          }
        }
      }
      signifs = TRUE;#--default
      sigidxs = -1;
      sigvals = -1;
      if (length(mdl_smry$s.pv)>0) {
        for(j in 1:length(mdl_smry$s.pv)){
          if(mdl_smry$s.pv[j]>0.05){
            signifs = FALSE;
            sigidxs = j;
            sigvals = mdl_smry$s.pv[j];
            break;
          }
        }#--i
      }
    }#--is.null(mdl_smry)
    
    #--test for concurvity
    concrv_tst = TRUE;
    concrv_idx = -1;
    concrv_val = -1;
    tryCatch(
      { concrv<-mgcv::concurvity(mdlres);
        if (is.null(concrv)){
          if (debug) cat("\nNULL: concurvity calculation failed in model",i,"\n");
        } else {
          concrv = concrv[concrv_opt,];
          if (debug) {cat("\nConcurvity:\n"); print(concrv);}
          for(j in 1:length(concrv)){
            concrv_val = max(concrv_val,concrv[j]);
            if(concrv[j]>0.5) {
              concrv_tst = FALSE;
              concrv_idx = j;
              break;
            }
          }#--i
        }#--is.null(concrv)
      }, 
      warning = function(w) {}, 
      error   = function(e) {
                              concrv_tst<-FALSE;
                              if (debug) cat("Caught error in concurvity check\n");
                            }, 
      finally = {}
    )#--tryCatch
    if (debug) cat("#################################\n\n");
    
    tbl = tibble::tibble(i=i,
                         frmla=as.character(frmla),
                         smths=paste(mdlsmths[idx],collapse="+"),
                         rsqr = rsqr,
                         aic=aic,
                         signifp=signifp,
                         sigdxp=sigidxp,
                         sigvalp=sigvalp,
                         signifs=signifs,
                         sigdxs=sigidxs,
                         sigvals=sigvals,
                         concrv_tst=concrv_tst,
                         concrv_idx=concrv_idx,
                         concrv_val=concrv_val);
    tbl
  }#--foreach(i = 1:ncmbs)
  if (doParallel) stopCluster(cl);
  
  dfrAll = dplyr::bind_rows(lstAll);
  return(dfrAll);
}

#--load full model to evaluate
mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;

#--evaluate all models
dfrAll = evalAll(mdl,
                  concrv_opt=2,
                  doParallel=TRUE,
                  ncores=parallel::detectCores()-1,
                  max_ncmbs=NULL,
                  icmbs=NULL,
                  debug=FALSE);
dfrAllp = dfrAll |> dplyr::arrange(dplyr::desc(concrv_tst),aic);

#--evaluate
best_idx = dfrAllp$i[1];

################################################################################
#--for debugging
if(FALSE){
  resp = as.character(formula(mdl))[2];     #--model response
  mdldata   = mdl$model;                    #--data for model fit
  mdlfmly   = family(mdl);                  #--model family object
  mdlvars   = gratia::model_vars(mdl);      #--model covariates
  mdlsmths  = gratia::smooths(mdl);         #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);      #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);    #--list, by smooth, of variables involved in each smooth
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  mdlbss    = getSmoothsMarginalBasisTypes(formula(mdl));#--basis types for smooths
  offset_   = model.offset(mdl);  #--model offsets, if any
  weights_  = model.weights(mdl); #--model weights, if any
  cmbs = getCombs(mdlsmths);
  ncmbs = nrow(cmbs);
  
  i = 3;
    idx = as.logical(as.vector(t(cmbs[i,])));
    cat(paste0("\n\nEvaluating model ",i,": ",paste(mdlsmths[idx],collapse=" + "),"\n"));
    frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],10,8);
    #cat("offset_:",offset_,"\n");
    #cat("weights_:",weights_,"\n");
    mdlres = mgcv::gam(formula=frmla,
                       family=mdlfmly,
                       data=mdldata,
                       offset=offset_,
                       weights=weights_,
                       method="REML");
    kchk = k.check(mdlres);#--not really doing anything with k.check result here
    {cat("\nk.check:\n"); print(kchk);}
    mdl_smry<-summary(mdlres);
    if (is.null(mdl_smry)){
      cat("\nNULL: model summary failed in model",i,"\n");
      signifp = FALSE;#--default
      sigidxp = -1;
      sigvalp = -1;
      signifs = FALSE;#--default
      sigidxs = -1;
      sigvals = -1;
    } else {
      {cat("\nSummary:\n"); print(mdl_smry);}
      rsqr      = mdl_smry$r.sq;
      aic      = AIC(mdlres);
      #--check for non-significant terms
      signifp = TRUE;#--default
      sigidxp = -1;
      sigvalp = -1;
      if (length(mdl_smry$p.pv)>0) {
        for(j in 1:length(mdl_smry$p.pv)){
          if(mdl_smry$p.pv[j]>0.05){
            signifp = FALSE;
            sigidxp = j;
            sigvalp = mdl_smry$p.pv[j];
            break;
          }
        }
      }
      signifs = TRUE;#--default
      sigidxs = -1;
      sigvals = -1;
      if (length(mdl_smry$s.pv)>0) {
        for(j in 1:length(mdl_smry$s.pv)){
          if(mdl_smry$s.pv[j]>0.05){
            signifs = FALSE;
            sigidxs = j;
            sigvals = mdl_smry$s.pv[j];
            break;
          }
        }#--i
      }
    }#--is.null(mdl_smry)
    
    #--test for concurvity
    concrv_tst = TRUE;
    concrv_idx = -1;
    concrv_val = -1;
    tryCatch(
      { concrv<-mgcv::concurvity(mdlres);
        if (is.null(concrv)){
          cat("\nNULL: concurvity calculation failed in model",i,"\n");
        } else {
          concrv = concrv[concrv_opt,];
          {cat("\nConcurvity:\n"); print(concrv);}
          for(j in 1:length(concrv)){
            concrv_val = max(concrv_val,concrv[j]);
            if(concrv[j]>0.5) {
              concrv_tst = FALSE;
              concrv_idx = j;
              break;
            }
          }#--i
        }#--is.null(concrv)
      }, 
      warning = function(w) {}, 
      error   = function(e) {
                              concrv_tst<-FALSE;
                              if (debug) cat("Caught error in concurvity check\n");
                            }, 
      finally = {}
    )#--tryCatch
    cat("#################################\n\n");
}

