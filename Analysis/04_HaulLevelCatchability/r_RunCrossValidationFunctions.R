#--Function to run cross validation across a set of models
#----criteria include CONCURVITY and AIC
require(ggplot2);
require(gratia);
require(mgcv);
require(tibble);
#--runCrossValidation code-------------------------------------------------
#' 
#' @param mdl - (possibly unevaluated) model with all terms 
#' @param ks - vector with max k-values for smooths (length = max dimension of smooths)
#' @param dfrData - complete dataframe
#' @param col_link - column name in dfrTest to use as values for testing predicting capability
#' @param startFold - fold to start with (default: 1; >1 to resume if analysis was prematurely interrupted were)
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
#' @return NULL (invisibly). Results are saved as "fold_`fld`.RData" in the working directory. 
#' 
runCrossValidation<-function(mdl,
                             ks,
                             dfrData,
                             col_link="lnR",
                             startFold=1,
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
  nothing = NULL;
  col_offsets = as.character(mdl$cl["offset"]);  
  if (col_offsets=="NULL") {  #--NULL doesn't work as a symbol
    col_offsets = "default_offsets";
    dfrData[[col_offsets]] = 0;
  }
  col_weights = as.character(mdl$cl["weights"]); 
  if (col_weights=="NULL") {  #--NULL doesn't work as a symbol
    col_weights = "default_weights";
    dfrData[[col_weights]] = 1;
  }
  sym_col_offsets = sym(col_offsets);
  sym_col_weights = sym(col_weights);
  #--gam options
  str_method  = as.character(mdl$cl["method"]);  if (str_method=="NULL")  str_method  = NULL;
  str_select  = as.character(mdl$cl["select"]);  lgl_select = ifelse(str_select=="NULL",FALSE,as.logical(str_select));
  
  #--run model selection across numFolds 
  lstRes = list();
  for (fld in startFold:numFolds){
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
                           debug=debug) |> 
              dplyr::mutate(fold=fld,.before=1);
    wtsUtilities::saveObj(dfrAll,paste0("fold_",fld,".RData"));
    rm(dfrAll);
  }
  return(invisible(NULL));
}





