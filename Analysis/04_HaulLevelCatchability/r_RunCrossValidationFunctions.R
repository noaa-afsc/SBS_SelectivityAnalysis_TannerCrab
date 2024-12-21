#--Functions to evaluate "best" model in a model hierarchy
#----criteria include CONCURVITY and AIC
require(ggplot2);
require(gratia);
require(mgcv);
require(tibble);
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

#--test: eval best model-------------------------------------------------
if (FALSE){
  #--evaluate best model
  best_idx = dfrAll$i[1];#--index of best model in evaluated combinations
  best_mdl = evalBestModel(mdl,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
}

#--runCrossValidation code-------------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param dfrData - complete dataframe
#' @param link - column name in dfrTest to use as values for testing predicting capability
#' @param numFolds - number of cross validation folds to evaluate
#' @param selectBy - method to select testing data within a fold ("by_haul" or "by_obs")
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
                             col_link="lnR",
                             numFolds=10,
                             selectBy="by_haul",
                             concrv_opt=2,
                             doParallel=TRUE,
                             ncores=parallel::detectCores()-1,
                             max_ncmbs=NULL,
                             icmbs=NULL,
                             log=TRUE,
                             debug=FALSE
                             ){
  if (selectBy=="by_haul"){
    #--create testing/training splits for folds by selecting hauls
    #--determine unique paired hauls
    uniqHauls = dfrData |> 
                  dplyr::select(h) |> 
                  dplyr::distinct(h);
    nUHs = nrow(uniqHauls);
    #--randomly assign hauls to folds
    uniqHauls = uniqHauls |> 
                  dplyr::arrange(floor(runif(nUHs,1,nUHs)));
    nTHs = floor(nUHs/numFolds);#--number of testing hauls/fold
    uniqHauls[["fold"]] = rep_len(1:numFolds,nUHs);
    uniqHauls = uniqHauls |> dplyr::arrange(fold,h);
  } else if (selectBy=="by_obs"){
    #--create testing/training splits for folds by selecting individual observations
    uniqObs = dfrData |> 
                dplyr::select(y,h,z) |> 
                dplyr::distinct(y,h,z);
    nUOs = nrow(uniqObs);
    #--randomly assign individuals to folds
    uniqObs = uniqObs |> 
                  dplyr::arrange(floor(runif(nUOs,1,nUOs)));
    nTOs = floor(nUOs/numFolds);#--number of testing observations/fold
    uniqObs[["fold"]] = rep_len(1:numFolds,nUOs);
    uniqObs = uniqObs |> dplyr::arrange(fold,y,h,z);
  } else {
    stop(paste0("'",selectBy,"' is not a recognized option for splitting testing/training datasets.\n",
                "Valid options are 'by_haul' and 'by_obs'."));
  }
  
  #--column names for training datasets
  col_offsets = as.character(mdl$cl["offset"]);  if (col_offsets=="NULL") col_offsets = NULL;
  col_weights = as.character(mdl$cl["weights"]); if (col_weights=="NULL") col_weights = NULL;
  sym_col_offsets = sym(col_offsets);
  sym_col_weights = sym(col_weights);
  #--gam options
  str_method  = as.character(mdl$cl["method"]);  if (str_method=="NULL")  str_method  = NULL;
  str_select  = as.character(mdl$cl["select"]);  lgl_select = ifelse(str_select=="NULL",FALSE,as.logical(str_select));
  
  #--run model selection across numFolds 
  lstRes = list();
  for (fld in 1:numFolds){
    #--testing: fld = 1;
    cat("Running fold",fld,"\n");
    if (selectBy=="hauls"){
      #--select testing set by hauls (h) in fold, training set by obverse
      hlsTst = (uniqHauls |> dplyr::filter(fold==fld))[["h"]];
      dfrTst = dfrData |> dplyr::filter(h %in% hlsTst);
      dfrTrn = dfrData |> dplyr::filter(!(h %in% hlsTst));
    } else {
      #--select testing set by observations in fold (y,h,z is unique), training set by obverse
      iosTst = (uniqObs |> dplyr::filter(fold==fld));#--dataframe w/ unique obs by y,h,z for testing
      dfrTst = dfrData |> dplyr::semi_join(iosTst,by=c("y","h","z"));#--rows in dfrData that match iosTst, for testing
      dfrTrn = dfrData |> dplyr::anti_join(iosTst,by=c("y","h","z"));#--rows in dfrData that don't match iosTst, for training
    }
    cat("dfrTrn has",nrow(dfrTrn),"rows\n");
    cat("dfrTst has",nrow(dfrTst),"rows\n");
    frmla = frmla<-formula(paste(as.character(formula(mdl))[c(2,1,3)],collapse=" "),
                           env=environment());#--need to evaluate formula in *this* environment
    mdlp = rlang::inject(
             mgcv::gam(frmla,
                     family=family(mdl),
                     data=dfrTrn,
                     offset=!!sym_col_offsets,
                     weights=!!sym_col_weights,
                     method=str_method,
                     select=lgl_select,
                     fit=FALSE)
           );
    print(mdlp$cl);
    logfile=NULL;
    if (log||debug) {
      logfile = paste0("log_",wtsUtilities::formatZeros(fld),".txt");
      if (file.exists(logfile)) file.remove(logfile);
    }
    dfrAll = evalAllModels(mdlp,
                           ks,
                           dfrTrn,
                           dfrTst,
                           col_link=col_link,
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




