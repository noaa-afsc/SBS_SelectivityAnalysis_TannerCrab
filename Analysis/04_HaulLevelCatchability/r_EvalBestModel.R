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
                        i,
                        verbose=FALSE){
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
  if (verbose) cat(paste0("\n\nEvaluating best model (",i,"): ",paste(mdlsmths[idx],collapse=" + "),"\n"));
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

