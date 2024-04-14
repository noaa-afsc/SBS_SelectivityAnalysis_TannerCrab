#--Functions to evaluate "best" model in a model hierarchy
#----criteria include CONCURVITY and AIC

require(ggplot2);
require(gratia);
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
#' 
#' @return dataframe with all possible combinations of 
#' smooth terms involving the first term. Each column represents a 
#' separate smooth term. Each row represents a 
#' potential model, with included smooth terms having a value of 1. 
#' Excluded terms have a vlaue of 0.
#' 
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

#--evalBestModel code-------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param i - index of "best" model from [evalAllModels] 
#' 
#' @return fitted mgcv::gam model representing "best" combination of terms 
#' from the input model `mdl`.
#' 
evalBestModel<-function(mdl,
                        ks,
                        i){
  #--data for model fit
  if (class(mdl)=="gam.prefit"){
    mdldata   = mdl$mf;
  } else {
    mdldata = mdl$model;
  }
  resp = as.character(formula(mdl))[2];     #--model response
  mdlfmly   = mdl$family;                   #--model family object
  offset_   = mdl$offset;                   #--model offsets, if any
  weights_  = mdl$w;                        #--model weights, if any
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
  cat(paste0("\n\nEvaluating best model (",i,"): ",paste(mdlsmths[idx],collapse=" + "),"\n"));
  frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],ks,
                              env=environment());
  best   = mgcv::gam(formula=frmla,
                     family=mdlfmly,
                     data=mdldata,
                     offset=offset_,
                     weights=weights_,
                     method="REML");
  return(best);
}

#--evalAllModels code-------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param dfrTest - dataframe for testing out-of-sample predictive capability
#' @param link - column name in dfrTest to use as values for testing predicting capability
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
evalAllModels<-function(mdl,
                        ks,
                        dfrTest,
                        link="lnR",
                        concrv_opt=2,
                        doParallel=FALSE,
                        ncores=parallel::detectCores()-1,
                        max_ncmbs=NULL,
                        icmbs_=NULL,
                        logfile=NULL,
                        debug=FALSE){
  t1 = Sys.time();
  
  if (class(mdl)=="gam.prefit"){
    mdldata   = mdl$mf;
  } else {
    mdldata = mdl$model;
  }
  resp = as.character(formula(mdl))[2];     #--model response
  offset_   = mdl$offset;  #--model offsets, if any
  weights_  = mdl$w;       #--model weights, if any
  mdlfmly   = family(mdl);                  #--model family object
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
  
  #--set up logging
  log = FALSE;
  if (!is.null(logfile)) log = TRUE;
  if (debug&&is.null(logfile)) logfile = "log.txt";
  if (log||debug) writeLines(paste0("Started evalAllModels at ",t1), logfile);#--creates or overwrites logfile
  
  lstAll<-foreach (i = icmbs,
                   .inorder=debug,
                   .export=c("createModelFormula"),
                   .packages=c("mgcv","tibble")) %dopar% {
    if (log||debug) sink(logfile, append=TRUE);
    idx = as.logical(as.vector(t(cmbs[i,])));
    if (log||debug) cat(paste0("\n\nEvaluating model ",i,": ",paste(mdlsmths[idx],collapse=" + "),"\n"));
    frmla  = createModelFormula(resp,mdlsmths[idx],mdlsmtrms[idx],mdlbss[idx],ks,
                                env=environment());
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
    rsqr_prd = 0;
    if (!is.null(dfrTest)) {
      prd = predict(mdlres,dfrTest,type="link");
      obs = dfrTest[[link]];
      rsqr_prd = cor(prd,obs)^2;
      mspe_prd = mean((obs-prd)^2);                          #--mean square prediction error
      mase_prd = mean(abs(obs-prd))/mean(abs(obs-mean(obs)));#--mean absolute scaled error
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
                         concrv_val=concrv_val);
    if (log) cat(i,paste(str[2],"~",str[3]),rsqr,aic,rsqr_prd,signifp,signifs,concrv_tst,"\n",sep=", ");
    if (debug) cat("#################################\n\n");
    if (log||debug) sink();
    tbl
  }#--foreach(i = 1:ncmbs)
  if (doParallel) stopCluster(cl);
  
  if (log||debug) {
    t2 = Sys.time();
    cat(paste0("Finished evalAllModels logging at ",t2,"\n"), file=logfile, append=TRUE);
    cat(paste0("Elapsed time was ",t2-t1,"\n"),               file=logfile, append=TRUE);
  }

  
  #--create dataframe and order by 
  #----1. passing concurvity test
  #----2. AIC value
  dfrAll = dplyr::bind_rows(lstAll) |> 
             dplyr::arrange(dplyr::desc(concrv_tst),mase_prd);
  return(dfrAll);
}

#--runCrossValidation code-------------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param dfrData - complete dataframe
#' @param link - column name in dfrTest to use as values for testing predicting capability
#' @param numFolds - number of cross validation folds to evaluate
#' @param concrv_opt - option (1=worst, 2=observed) for type to use for concurvity criterion
#' @param doParallel - flag to do parallel processing
#' @param ncores - number of cores to use for parallel processing
#' @param max_ncmbs - max number of potential models to run (NULL results in running all)
#' @param icmbs_ - vector of indices from model combinations dataframe for specific models to evaluate 
#' @param log - flag to write to log 
#' @param debug - flag to print 
#' 
#' @return unordered dataframe 
#' 
runCrossValidation<-function(mdl,
                             ks,
                             dfrData,
                             link="lnR",
                             numFolds=10,
                             concrv_opt=2,
                             doParallel=TRUE,
                             ncores=parallel::detectCores()-1,
                             max_ncmbs=NULL,
                             icmbs=NULL,
                             log=TRUE,
                             debug=FALSE
                             ){
  #--determine unique paired hauls
  uniqHauls = dfrData |> 
                dplyr::select(h) |> 
                dplyr::distinct(h);
  nUHs = nrow(uniqHauls);
  #--randomly assign hauls to folds
  uniqHauls = uniqHauls |> 
                dplyr::arrange(floor(runif(nUHs,1,nUHs)));
  nTHs = floor(nUHs/numFolds)
  uniqHauls[["fold"]] = rep_len(1:numFolds,nUHs);
  uniqHauls = uniqHauls |> dplyr::arrange(fold,h);
  
  #--run model selection across numFolds 
  lstRes = list();
  for (fld in 1:numFolds){
    #--testing: fld = 1;
    hlsTst = (uniqHauls |> dplyr::filter(fold==fld))[["h"]];
    dfrTst = dfrData |> dplyr::filter(h %in% hlsTst);
    dfrTrn = dfrData |> dplyr::filter(!(h %in% hlsTst));
    mdlp = mgcv::gam(mdl$formula,family=family(mdl),data=dfrTrn,
                     select=TRUE,method="REML",fit=FALSE);
    logfile=NULL;
    if (log||debug) {
      logfile = paste0("log_",wtsUtilities::formatZeros(fld),".txt");
      if (file.exists(logfile)) file.remove(logfile);
    }
    dfrAll = evalAllModels(mdlp,
                           ks,
                           dfrTst,
                           link=link,
                           concrv_opt=concrv_opt,
                           doParallel=doParallel,
                           ncores=ncores,
                           max_ncmbs=max_ncmbs,
                           icmbs=icmbs,
                           logfile=logfile,
                           debug=debug);
    lstRes[[fld]] = dfrAll |> dplyr::mutate(fold=fld,.before=1);
  }
  dfrRes = dplyr::bind_rows(lstRes);
  return(dfrRes);
}

#--plot model results------------------------------------------------------
plt1D<-function(sme,prs,x_,y_){
  require(rlang);
  vx=sym(x_); vy=sym(y_);
  gratia::draw(sme) + 
        geom_point(aes(x=!!vx,y=!!vy),prs,alpha=0.2);
}
plt2D<-function(sme,prs,x_,y_,z_){
  require(rlang);
  vx=sym(x_); vy=sym(y_); vz=sym(z_);
  gratia::draw(sme) + 
        geom_point(aes(x=!!vx,y=!!vy,size=!!vz),prs,alpha=0.2) + 
        scale_size_area()
}
getModelPlots<-function(mdl,
                        subs = c("z"="size (mm CW)",
                                 "d"="depth (m)",
                                 "t"="temperature (deg C)",
                                 "f"="phi",
                                 "s"="sorting coefficient")){
  resp = as.character(formula(mdl))[2];     #--model response
  fmly = family(mdl);
  nsms  = n_smooths(mdl);
  smths = smooths(mdl);
  smtsl = smooth_terms(mdl);
  dfrDatp = mdl$model;
  allPlts = list();
  for (ism in 1:nsms){
    #--testing: ism = 1;
    smth = smths[ism];
    smts = smtsl[[ism]];
    sme = gratia::smooth_estimates(mdl,select=smth,
                                   unconditional=TRUE,
                                   overall_uncertainty=TRUE);
    sme= add_confint(sme,0.80);
    prs = add_partial_residuals(dfrDatp,mdl,select=smth);
    plts = list();
    if (length(smts)==1){
      p = plt1D(sme,prs,smts[1],smth) + 
                    labs(x=subs[smts[1]],y="Partial effect on link scale");
      cap = paste0("Partial effect of ",smts[1]," for the link-scale response ",
                    resp," assuming a ",fmly$family," error distribution and a ",
                    fmly$link," link.")
      plts[[1]] = list(p=p,cap=cap);
      if (smth=="ti(z)") {
        p = gratia::draw(transform_fun(sme,fun=exp,constant=model_constant(mdl))) + 
                      labs(x=subs[smts[1]],y="Base Selectivity");
        cap = paste0("Estimated base selectivity (independent of covariate effects).");
        plts[[2]] = list(p=p,cap=cap);
      } else {
        # plts[[2]] = gratia::draw(transform_fun(sme,fun=inv_link(mdl))) + 
        #               labs(x=subs[smts[1]],y="Response Scale Effect");
        # cap = paste0("Estimated selectivity (independent of covariate effects.");
        # plts[[2]] = list(p=p,cap=cap);
      }
    } else if (length(smts==2)){
      p = plt2D(sme,prs,smts[1],smts[2],smth) + 
                    labs(x=subs[smts[1]],y=subs[smts[2]],size="Partial\nresiduals");
      cap = paste0("Partial effect of ",smts[1]," and ",smths[2]," for the link-scale response ",
                    resp," assuming a ",fmly$family," error distribution and a ",
                    fmly$link," link.");
      plts[[1]] = list(p=p,cap=cap);
    }
    allPlts[[smth]] = list(plts);
    rm(plts);
  }
  return(allPlts)
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

#--test: eval best model-------------------------------------------------
  #--evaluate best model
  best_idx = dfrAll$i[1];#--index of best model in evaluated combinations
  best_mdl = evalBestModel(mdl,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
  
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



