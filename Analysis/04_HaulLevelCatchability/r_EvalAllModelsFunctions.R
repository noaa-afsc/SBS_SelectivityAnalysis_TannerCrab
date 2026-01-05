#--functions to evaluate "all" models----
require(mgcv);
require(tibble);
require(foreach);
require(doParallel);
require(parallel);
#--getSmoothsMarginalBasisTypes code-------------------------------------------
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
  #--for testing
  dirPrj = rstudioapi::getActiveProject();
  mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;
  bss = getSmoothsMarginalBasisTypes(formula(mdl));
}

#--createModelFormula code-------------------------------------------
#' 
#' @param resp - response name
#' @param smths - vector of smooths names 
#' @param smtrms - list (by smooth) of smooth terms 
#' @param bss  - vector of marginal basis types (from getSmoothMarginalBasisTypes)
#' @param ks - max k's for smooths 
#' 
#' @return character string representing a model formula 
#'
createModelFormula<-function(resp,smths,smtrms,bss,ks,
                             env=parent.frame()){
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
      strp = paste0(strp,ks[1],")) ");
    } else {
      strp = paste0(strp,paste(ks[1:2],collapse=","),")) ");
    }
    if (i<n) str = paste(str,strp,"+");
  }
  str = paste(str,strp);
  return(as.formula(str,env=env));
}#--createModelFormula

#--test createModelFormula----
if (FALSE){
  #--for testing
  dirPrj = rstudioapi::getActiveProject();
  mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;
  resp      = as.character(formula(mdl))[2];     #--model response
  mdlvars   = gratia::model_vars(mdl);           #--model covariates
  mdlsmths  = gratia::smooths(mdl);              #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);           #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);         #--list, by smooth, of variables involved in each smooth
  bss = getSmoothsMarginalBasisTypes(formula(mdl));
  mdlprtms  = gratia::parametric_terms(mdl);     #--vector of parametric terms names
  mdlfmly   = mdl$family;                        #--model family object
  createModelFormula(resp,smths,smtrms,bss,k1d=10,k2d=5)
}
  
#--getCombs code-------------------------------------------
#' 
#' @param mdlsmths - vector of model smooths 
#' @param verbose - flag to print details
#' 
#' @return dataframe with all possible combinations of 
#' smooth terms involving the first term. Each column represents a 
#' separate smooth term. Each row represents a 
#' potential model, with included smooth terms having a value of 1. 
#' Excluded terms have a vlaue of 0.
#' 
getCombs<-function(mdlsmths,verbose=FALSE){
  smbase = mdlsmths[1]; # required (base) terms
  smxtra = mdlsmths[2:length(mdlsmths)];# potential (extra) terms to evaluate
  n_xtra = length(smxtra);# number of potential (extra)  terms to evaluate
  cmbs = tibble::tibble(x=1); names(cmbs)[ncol(cmbs)] = smbase;
  for (i in 1:n_xtra) {
    cmbs = tidyr::expand_grid(cmbs,tibble::tibble(x=0:1))
    names(cmbs)[ncol(cmbs)] = smxtra[i];
  }
  ncmbs = nrow(cmbs);
  if (verbose) cat("Total number of combinations:",ncmbs,"\n");
  return(cmbs);
}

#--evalAllModels code-------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param dfrTrain - dataframe with training data
#' @param dfrTest - dataframe for testing out-of-sample predictive capability
#' @param col_obs - column name in dfrTest to use as observed values for testing predicting capability
#' @param col_offsets - column name in dfrTest to use as link-scale offsets (or NULL if no offsets)
#' @param col_wgts - column name in dfrTest to use as weights (or NULL if no weights)
#' @param col_link - column name in dfrTest to use as values for testing predicting capability
#' @param concrv_opt - option (1=worst, 2=observed) for type to use for concurvity criterion
#' @param doParallel - flag to do parallel processing
#' @param ncores - number of cores to use for parallel processing
#' @param max_ncmbs - max number of potential models to run (NULL results in running all)
#' @param icmbs_ - vector of indices from model combinations dataframe for specific models to evaluate 
#' @param logfile - NULL, or name of file to write logging or debugging info to
#' @param debug - flag to print debugging info to 
#' 
#' @return unordered dataframe with model evaluations 
#' 
#' @details Returned dataframe has columns
#' \itemize{
#'  \item{fold - fold index}       
#'  \item{i - model index}          
#'  \item{frmla - model formula}      
#'  \item{smths - smooths in model}      
#'  \item{rsqr - adjusted R-sqr from model summary}       
#'  \item{aic - AIC}        
#'  \item{rsqr_prd - R-sqr prediction error (covariance) on the link scale}   
#'  \item{mspe_prd - mean square absolute prediction error on the link scale}   
#'  \item{mase_prd - mean absolute scaled error on the link scale}   
#'  \item{signifp - flag if all parametric terms are significant} 
#'  \item{sigdxp - index of first(?) non-significant parametric term}     
#'  \item{sigvalp - p-value of first(?) non-significant parametric term}    
#'  \item{signifs - flag if all smooth terms are significant}    
#'  \item{sigdxs - index of first(?) non-significant smooth term}     
#'  \item{sigvals - p-value of first(?) non-significant smooth term}    
#'  \item{concrv_tst - flag if all concurvity tests pass} 
#'  \item{concrv_idx - index of first(?) concurvity test that failed} 
#'  \item{concrv_val - value of first(?) concurvity test that failed}
#'  \item{test - list}
#' }
#' 
evalAllModels<-function(mdl,
                        ks,
                        dfrTrain,
                        dfrTest,
                        col_link="lnR",
                        concrv_opt=2,
                        doParallel=FALSE,
                        ncores=parallel::detectCores()-1,
                        max_ncmbs=NULL,
                        icmbs_=NULL,
                        logfile=NULL,
                        verbose=TRUE,
                        debug=FALSE){
  cat("in evalAllModels\n");
  t1 = Sys.time();
  
  #--get information for fitting models with training dataset
  mdl_frmla = formula(mdl);
  mdl_env   = environment(mdl_frmla);
  resp = as.character(mdl_frmla)[2];        #--model response
  # offset_   = mdl$offset;                   #--vector of model offsets, if any
  # weights_  = mdl$w;                        #--vector of model weights, if any
  mdlfmly   = family(mdl);                  #--model family object
  mdlvars   = gratia::model_vars(mdl);      #--model covariates
  mdlsmths  = gratia::smooths(mdl);         #--vector of the names of the smooths
  mdlsmdims = gratia::smooth_dim(mdl);      #--vector of the number of dimensions for each smooth
  mdlsmtrms = gratia::smooth_terms(mdl);    #--list, by smooth, of variables involved in each smooth
  mdlprtms  = gratia::parametric_terms(mdl);#--vector of parametric terms names
  mdlbss    = getSmoothsMarginalBasisTypes(mdl_frmla);#--basis types for smooths
  #--column names for training dataset
  #----need to get column name for offsets and weights from the model call
  col_offsets = as.character(mdl$cl["offset"]);  #if (col_offsets=="NULL") col_offsets = NULL;
  col_weights = as.character(mdl$cl["weights"]); #if (col_weights=="NULL") col_weights = NULL;
  sym_col_offsets = sym(col_offsets);
  sym_col_weights = sym(col_weights);
  if (verbose) cat("offsets, weights in columns",col_offsets,col_weights,"\n");
  #--gam options (same as "full" model)
  str_method = mdl_env$str_method;
  if (verbose) cat("str_method =",str_method,"\n")
  lgl_select = mdl_env$lgl_select;
  if (verbose) cat("lgl_select =",lgl_select,"\n")
  
  #--get column names for testing dataset
  col_respons = resp;#--same as for training dataset

  #--set up models to run
  cmbs = getCombs(mdlsmths);
  ncmbs = nrow(cmbs);
  if (is.null(max_ncmbs)) max_ncmbs = ncmbs;
  cat("Evaluating",max_ncmbs,"combinations.\n");
  icmbs = 1:max_ncmbs;
  if (!is.null(icmbs_)) icmbs = icmbs_;

  #--set up parallel processing
  if (doParallel){
    cat("setting up parallel procesing\n");
    cl <- makeCluster(ncores);
    registerDoParallel(cl)
  }
  
  #--set up logging
  log = FALSE;
  if (!is.null(logfile)) log = TRUE;
  if (debug&&is.null(logfile)) logfile = "log.txt";
  if (log||debug) {
    writeLines(paste0("Started evalAllModels at ",t1), logfile);#--creates or overwrites logfile
    if (!doParallel) cat(paste0("Started evalAllModels at ",t1,"\n"));
  }
  
  lstAll<-foreach(i = icmbs,
                 .errorhandling="remove",
                 .inorder=debug,
                 .export=c("createModelFormula","calcLogLike.binomial","calcLogLike.tw"),
                 .packages=c("mgcv","tibble")) %dopar% {
    if (doParallel&&(log||debug)) sink(logfile, append=TRUE);
    if (log||debug) cat("\n\nicmbs for next model:",paste(as.vector(t(cmbs[i,])),collapse=" "),"\n");
    idx = as.logical(as.vector(t(cmbs[i,])));
    if (log||debug) cat(paste0("Evaluating model ",i,": ",paste(mdlsmths[idx],collapse=" + "),"\n"));
    frmla  = createModelFormula(resp,
                                mdlsmths[idx],
                                mdlsmtrms[idx],
                                mdlbss[idx],
                                ks,
                                env=environment());
    #--fit candidate model with training dataset
    mdlres = NULL;
    tryCatch(
      {mdlres = rlang::inject(
                  mgcv::gam(formula=frmla,
                            family=mdlfmly,
                            data=dfrTrain,
                            offset=!!sym_col_offsets,
                            weights=!!sym_col_weights,
                            method=str_method,
                            select=lgl_select,
                            na.action=na.omit,
                            fit=TRUE)
               );},
      warning = function(w) {cat("Caught warning calculating mdlres\n");}, 
      error   = function(e) {
                              cat("Caught error calculating mdlres\n");
                              cat(as.character(e));
                            }, 
      finally = function(c){cat("Caught condition calculating mdlres.\n");
                            cat(conditionMessage(c));}
    );
    if (is.null(mdlres)){
      if (log||debug) cat("model",i,"result is NULL.\n"); 
      tbl = NULL;
    } else {
      if (log||debug) {cat("model",i,"result:\n"); print(mdlres);}
      mdl_smry<-summary(mdlres);
      if (is.null(mdl_smry)){
        if (log||debug) cat(paste("\nNULL: model summary failed in model",i,"\n"));
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
        rsqr = mdl_smry$r.sq;
        aic  = AIC(mdlres);
        prd_rsp = as.vector(unname(predict(mdlres,dfrTrain,type="response")[]));
        obs_rsp = dfrTrain[[col_respons]];
        if (mdlres$family$family=="binomial"){
          llsTrn = calcLogLike.binomial(obs_rsp,dfrTrain[[col_weights]],prd_rsp,family=mdlres$family,offsets=dfrTrain[[col_offsets]]);
        } else if (stringr::str_starts(mdlres$family$family,"Tweedie")){
          if (debug) cat("calculating llsTrn for Tweedie family. scale = ",mdlres$scale,"\n")
          llsTrn = calcLogLike.tw(obs_rsp,prd_rsp,family=mdlres$family,rho=log(mdlres$scale),offsets=dfrTrain[[col_offsets]]);
          if (debug) cat("calculated llsTrn for Tweedie family\n")
        } else {
          cat( paste0("evalAllModels failure: model family '",mdlres$family$family,"'is not recognized.\n"));
          stop(paste0("evalAllModels failure: model family '",mdlres$family$family,"'is not recognized."));
        }
        llMdl  = as.numeric(logLik(mdlres));
        llTrn  = sum(llsTrn);
        scrTrn = sum(llsTrn)/length(prd_rsp);
        tblTrn = tibble::tibble(obs_rsp=obs_rsp,
                                prd_rsp=prd_rsp,
                                obs_lnk=dfrTrain[[col_link]],
                                prd_lnk=mdlres$family$linkfun(prd_rsp),
                                offsets=dfrTrain[[col_offsets]],
                                ll=llsTrn);
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
            if (log||debug) cat(paste("\nNULL: concurvity calculation failed in model",i,"\n"));
          } else {
            concrv = concrv[concrv_opt,];
            if (log||debug) {cat(paste("\nConcurvity:\n")); print(concrv);}
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
                                if (log||debug) cat("Caught error in concurvity check\n");
                              }, 
        finally = {}
      )#--tryCatch
      prd_rsp  = NULL;
      obs_rsp  = NULL;
      rsqr_prd = NA_real_;
      mspe_prd = NA_real_;
      mase_prd = NA_real_;
      llTst    = NA_real_;
      scrTst   = NA_real_
      tblTst   = NULL;
      if (!is.null(dfrTest)) {
        prd_rsp = as.vector(unname(predict(mdlres,dfrTest,type="response")[]));
        obs_rsp = dfrTest[[col_respons]];
        if (mdlres$family$family=="binomial"){
          llsTst = calcLogLike.binomial(obs_rsp,dfrTest[[col_weights]],prd_rsp,family=mdlres$family,offsets=dfrTest[[col_offsets]]);
        } else if (stringr::str_starts(mdlres$family$family,"Tweedie")){
          if (debug) cat("calculating llsTst for Tweedie family\n")
          llsTst = calcLogLike.tw(obs_rsp,prd_rsp,family=mdlres$family,rho=log(mdlres$scale),offsets=dfrTest[[col_offsets]]);
          if (debug) cat("calculated llsTst for Tweedie family\n")
        } else {
          cat( paste0("evalAllModels failure: model family '",mdlres$family$family,"'is not recognized.\n"));
          stop(paste0("evalAllModels failure: model family '",mdlres$family$family,"'is not recognized."));
        }
        llTst  = sum(llsTst);
        scrTst = llTst/length(prd_rsp);
        tblTst = tibble::tibble(y=dfrTest$y,h=dfrTest$h,
                                d=dfrTest$d,t=dfrTest$t,f=dfrTest$f,s=dfrTest$s,
                                z=dfrTest$z,x=dfrTest$x,
                                n=dfrTest$n,p=dfrTest$p,
                                asr=dfrTest$asr,q=dfrTest$q,lnq=dfrTest$lnq,
                                obs_rsp=obs_rsp,
                                prd_rsp=prd_rsp,
                                obs_lnk=dfrTest[[col_link]],
                                prd_lnk=mdlres$family$linkfun(prd_rsp),
                                offsets=dfrTest[[col_offsets]],
                                ll=llsTst);
        #--R-sqr prediction error
        use = is.finite(prd_rsp)&is.finite(obs_rsp);
        prd = prd_rsp[use]; obs = obs_rsp[use];
        rsqr_prd = cor(prd,obs,use="na.or.complete")^2;
        #--mean square prediction error
        mspe_prd = mean((obs-prd)^2,na.rm=TRUE);
        #--mean absolute scaled error
        mase_prd = mean(abs(obs-prd),na.rm=TRUE)/mean(abs(obs-mean(obs,na.rm=TRUE)),na.rm=TRUE);
      }
      str = as.character(frmla);
      tbl = tibble::tibble(i=i,
                           frmla=paste(str[2],"~",str[3]),
                           smths=paste(mdlsmths[idx],collapse="+"),
                           rsqr = rsqr,
                           aic=aic,
                           rsqr_prd=rsqr_prd,
                           mspe_prd=mspe_prd,
                           mase_prd=mase_prd,
                           signifp=signifp,
                           sigdxp=sigidxp,
                           sigvalp=sigvalp,
                           signifs=signifs,
                           sigdxs=sigidxs,
                           sigvals=sigvals,
                           concrv_tst=concrv_tst,
                           concrv_idx=concrv_idx,
                           concrv_val=concrv_val,
                           llMdl=llMdl,
                           llTrn=llTrn,
                           scrTrn=scrTrn,
                           llTst=llTst,
                           scrTst=scrTst,
                           lstMdl=list(mdl=mdlres),
                           lstTrn=list(tblTrn=tblTrn),
                           lstTst=list(tblTst=tblTst));
      if (log||debug) cat(paste(paste0("\nResults for model ",i,": ",paste(str[2],"~",str[3])),rsqr,aic,rsqr_prd,signifp,signifs,concrv_tst,llTrn,scrTrn,llTst,scrTst,"\n",collapse=", "));
      print(head(tbl));
      if (log||debug) cat("#################################\n\n");
      if (doParallel&&(log||debug)) sink();
    }#--mdlres is NOT NULL
    tbl;#--return tbl as list element to lstAll
  }#--foreach(i = 1:ncmbs)
  if (doParallel) stopCluster(cl);
  
  if (log||debug) {
    t2 = Sys.time();
    cat(paste0("Finished evalAllModels logging at ",t2,"\n"), file=logfile, append=TRUE);
    if (!doParallel && verbose) cat(paste0("Elapsed time was ",t2-t1,"\n"),file=logfile, append=TRUE);
  }

  
  #--create dataframe
  dfrAll = dplyr::bind_rows(lstAll);
  return(dfrAll);
}

#--test: eval all models-------------------------------------------------
if (FALSE){
  dirPrj = rstudioapi::getActiveProject();
  
  #--load full model to evaluate
  mdl = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/Binomial_Males_CensoredData1/rda_mdlNrml_ZE2D.RData"))$model;
  
  #--evaluate all models
  dfrAll = evalAllModels(mdl,
                         ks=c(10,8),
                         concrv_opt=2,
                         doParallel=TRUE,
                         ncores=parallel::detectCores()-1,
                         max_ncmbs=NULL,
                         icmbs=NULL,
                         debug=FALSE);
}

#--For debugging evalAllModels-------------------------------------
if(FALSE){
  ks = c(10,8); #--max k values
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
    frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],ks);
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
