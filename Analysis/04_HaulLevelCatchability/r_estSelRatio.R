require(magrittr)
#'
#' @title Estimate selectivity ratio
#'
#' @importFrom dplyr mutate
#' @importFrom magrittr %<>%
#' @import mgcv
#' 
#' @return list with elements `model` and `results`. `model` is the fitted model object. 
#' `results` is the input dataframe (`dfr`) with columns `fits` and `residuals` added,
#' reflecting model fitted values (via `fitted(model)`) and 
#' deviance residuals (via `residuals(model,"residuals")`).
#' 
#' @details 
#' @export
#'
estSelRatio<-function(
            dfr,
            formula,
            family=stats::binomial(),
            method="REML",
            select=FALSE,
            offset="lnq",
            scale=0){
  if (family$family=="binomial"){
    cat("Running Binomial regression model.\n")
    if (is.null(offset)){
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                      method=method,select=select,scale=scale,offset=NULL);
    } else {
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=n,
                      method=method,select=select,scale=scale,offset=lnq);
    }
  } else if (family$family=="negative binomial"){
    cat("Running Negative Binuomial regression model.\n")
    mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                    method=method,select=select,scale=scale,offset=NULL);
  } else if (family$family=="gaussian"){
    cat("Running Gaussian regression model.\n")
    mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                    method=method,select=select,scale=scale,offset=NULL);
  } else if (family$family=="Beta regression"){
    cat("Running Beta regression model.\n");
    if (is.null(offset)){
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                      method=method,select=select,scale=scale,offset=NULL);
    } else {
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                      method=method,select=select,scale=scale,offset=lnq);
    }
  } else if (family$family=="Tweedie"){
    cat("Running Tweedie model.\n");
    if (is.null(offset)){
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                      method=method,select=select,scale=scale,offset=NULL);
    } else {
      mdl = mgcv::gam(data=dfr,family=family,formula=formula,weights=NULL,
                      method=method,select=select,scale=scale,offset=lnq);
    }
  }

  #--extract fitted values and calculate associated selectivity ratio
  if (family$family=="binomial"){
    dfr %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., prNMFS)
                           residuals=residuals(mdl,type="deviance"));#--deviance residuals for qqplot
  } else if (family$family=="gaussian") {
    dfr %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., lnR = logit(prNMFS)-lnq)
                           residuals=residuals(mdl,type="deviance"));#--deviance residuals for qqplot
  } else if (family$family=="negative binomial") {
    dfr %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., lnR = logit(prNMFS)-lnq)
                           residuals=residuals(mdl,type="deviance"));#--deviance residuals for qqplot
  } else if (family$family=="Beta regression") {
    dfr %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., lnR = logit(prNMFS)-lnq)
                           residuals=residuals(mdl,type="deviance"));#--deviance residuals for qqplot
  } else if (family$family=="Tweedie") {
    dfr %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., lnR = logit(prNMFS)-lnq)
                           residuals=residuals(mdl,type="deviance"));#--deviance residuals for qqplot
  }
  return(list(model=mdl,results=dfr));
}

newDFR<-function(z=NULL,d=NULL,t=NULL,f=NULL,s=NULL){
  return(tidyr::expand_grid(z,d,t,f,s));
}

plotEsts<-function(dfrDat,showPlots=FALSE){
  ps = list();
  #--plot fitted values
  ps[["1a"]] = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=fits,colour=y)) +
                ggplot2::geom_point() + ggplot2::geom_smooth() +
                ggplot2::labs(x="size (mm CW)",y="fitted prNMFS");
  if (showPlots) print(ps[["1a"]]);
  ps[["1b"]] = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=lnR,colour=y)) +
                ggplot2::geom_point() + ggplot2::geom_smooth() +
                ggplot2::labs(x="size (mm CW)",y="predicted lnR");
  if (showPlots) print(ps[["1b"]]);
  ps[["1c"]] = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=R,colour=y)) +
                ggplot2::geom_point() + ggplot2::geom_smooth() +
                ggplot2::labs(x="size (mm CW)",y="predicted R");
  if (showPlots) print(ps[["1c"]]);
  
  #--plot residuals
  ps[["2a"]]  = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=residuals,colour=y)) +
                 ggplot2::geom_point() + ggplot2::geom_smooth() +
                 ggplot2::labs(x="size (mm CW)",y="deviance residuals");
  if (showPlots) print(ps[["2a"]]);
  ps[["2b"]]  = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=fits,y=residuals,colour=y)) +
                 ggplot2::geom_point() + ggplot2::geom_smooth() +
                 ggplot2::labs(x="fitted values",y="deviance residuals");
  if (showPlots) print(ps[["2b"]]);
  ps[["2c"]]  = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(sample=residuals)) +
                 ggplot2::geom_qq() + ggplot2::geom_qq_line() +
                 ggplot2::labs(x="fitted values",y="deviance residuals");
  if (showPlots) print(ps[["2c"]]);

  return(ps);
}

plotPrdSel<-function(dfrNew,seFactor,showPlots=FALSE){
  ps = list();

  #--plot predicted selectivity values
  ps[["3a"]] = ggplot(data=dfrNew,mapping=aes(x=z,y=prdLnR,ymin=prdLnR-seLnR,ymax=prdLnR+seLnR)) + 
                geom_line() + geom_ribbon(alpha=0.4) + geom_vline(xintercept=c(25.0,180.0),linetype=2) +
                labs(x="size (mm CW)",y="ln(selectivity ratio)",subtitle=paste0(tolower(x),"s: ","predicted ln-scale R with +/- 1 SE."));
  if (showPlots) print(ps[["3a"]]);
  ps[["3b"]] = ggplot(data=dfrNew,mapping=aes(x=z,y=prdR,ymin=lower,ymax=upper)) + 
                geom_line() + geom_ribbon(alpha=0.4) + geom_vline(xintercept=c(25.0,180.0),linetype=2) +
                labs(x="size (mm CW)",y="selectivity ratio",subtitle=paste0(tolower(x),"s: ","predicted R with +/- ",seFactor," SE."))
  if (showPlots) print(ps[["3b"]]);
  
  return(ps);
}

#'
#' @title Calculate the predicted (ln-scale) selectivity ratio across a grid of covariates
#' @description Function to calculate the predicted (ln-scale) selectivity ratio across a grid of covariates.
#'
#' @param model - [mgcv::gam()] model object
#' @param dfrNew - expanded grid tibble with covariate values at which to calculate the selectivity
#' @param type - prediction type ("link"-the linear predictor [default],"terms","iterms","response")
#' @param terms - names of terms to include in results
#' @param exclude - names of terms to exclude from results
#' 
#' @return a tibble replicating dfrNew with additional columns `prd` and `se` 
#' representing the results of `predict.gam` with type=`type` and 
#' `unconditional=TRUE`.
#' 
#' @details See [mgcv::predict.gam()].
#' 
#' @import mgcv
#' @importFrom dplyr mutate
#' @export
#'
calcPrdSel<-function(model,dfrNew,type="link"){
  prd_ = predict(model,newdata=dfrNew,se.fit=TRUE,type=type,
                newdata.guaranteed=TRUE,unconditional=TRUE);
  dfrPrd = dfrNew |>
             dplyr::mutate(prd=prd_$fit,   #--offset is lnq = 0, identically
                           se=prd_$se.fit);#--standard error
  return(dfrPrd);
}

#'
#' @title Plot predicted selectivity estimates along 1 covariate
#'
#' @param dfrPrd - dataframe with columns prdLnR, seLnR, and covariates
#' @param x_ - covariate column to plot against
#' @param xlab - label for x axis
#' @param xints - x intercepts for vertical dashed reference lines (if not NULL)
#' @param dfrDat - original data to use for rug plots
#' @param ci - two-sided confidence interval
#' @return list of two plots (lnR, R vs.x_)
#' @importFrom wtsPlots getStdTheme
#' @import ggplot2
#' @export
#'
plotPrdSel1D<-function(dfrPrd,x_,xlab="size (mm CW)",xints=NULL,dfrDat=NULL,ci=0.80){
  dci = (1-ci)/2;
  mlt = pnorm(1-dci,0,1);#se multiplier corresponding to ci
  ps = list();

  #--plot predicted selectivity values
  p = ggplot(data=dfrPrd,mapping=aes(x={{x_}},y=prdLnR,ymin=prdLnR-mlt*seLnR,ymax=prdLnR+mlt*seLnR)) +
        geom_line() + geom_ribbon(alpha=0.4) +
        geom_hline(yintercept=log(c(0.5,1)),linetype=2) +
        labs(x=xlab,y="ln(selectivity ratio)",
             subtitle=paste0(tolower(x),"s: ","predicted ln-scale R with ",100*ci,"% CI.")) +
        wtsPlots::getStdTheme();
  if (!is.null(xints))  p = p + geom_vline(xintercept=xints,linetype=2);
  if (!is.null(dfrDat)) p = p + geom_rug(data=dfrDat,mapping=aes(x={{x_}}),inherit.aes=FALSE)
  ps[["lnR"]] = p;

  p = ggplot(data=dfrPrd,mapping=aes(x={{x_}},y=exp(prdLnR+(seLnR^2)/2),
                                     ymin=exp(prdLnR-mlt*seLnR),ymax=exp(prdLnR+mlt*seLnR))) +
        geom_line() + geom_ribbon(alpha=0.4) +
        geom_hline(yintercept=c(0,0.5,1),linetype=2) +
        labs(x=xlab,y="selectivity ratio",
             subtitle=paste0(tolower(x),"s: ","predicted R with ",100*ci,"% CI.")) +
        wtsPlots::getStdTheme();
  if (!is.null(xints))  p = p + geom_vline(xintercept=xints,linetype=2);
  if (!is.null(dfrDat)) p = p + geom_rug(data=dfrDat,mapping=aes(x={{x_}}),inherit.aes=FALSE)
  ps[["R"]] = p;
  return(ps);
}

#'
#' @title Plot predicted selectivity estimates across 2 covariates
#'
#' @param dfrPrd - dataframe with columns prdLnR, seLnR, and covariates
#' @param x_ - x-axis covariate column to plot against
#' @param y_ - y-axis covariate column to plot against
#' @param xlab - label for x axis
#' @param ylab - label for y axis
#' @param xints - x intercepts for vertical dashed reference lines (if not NULL)
#' @param dfrDat - original data to use for rug plots
#' @param ci - two-sided confidence interval
#' @return list of two plots (lnR, R vs.x_)
#' @importFrom wtsPlots getStdTheme
#' @import ggplot2
#' @export
#'
plotPrdSel2D<-function(dfrPrd,x_,y_,xlab="size (mm CW)",ylab="depth (m)",xints=NULL,dfrDat=NULL,ci=0.80){
  dci = (1-ci)/2;
  mlt = pnorm(1-dci,0,1);#se multiplier corresponding to ci
  ps = list();

  #--plot predicted selectivity values
  p = ggplot(data=dfrPrd,mapping=aes(x={{x_}},y={{y_}},z=prdLnR)) +
        geom_contour_filled(breaks=log(c(0.01,0.1,0.5,1,2,10))) +
        labs(x=xlab,y=ylab,
             subtitle=paste0(tolower(x),"s: ","predicted ln-scale R")) +
        wtsPlots::getStdTheme();
  if (!is.null(xints))  p = p + geom_vline(xintercept=xints,linetype=2);
  if (!is.null(dfrDat)) {
    p = p + geom_point(data=dfrDat,mapping=aes(x={{x_}},y={{y_}}),inherit.aes=FALSE,size=0.1,colour="white")
    p = p + geom_rug(data=dfrDat,mapping=aes(x={{x_}}),inherit.aes=FALSE);
    p = p + geom_rug(data=dfrDat,mapping=aes(y={{y_}}),inherit.aes=FALSE);
  }
  ps[["lnR"]] = p;

  p = ggplot(data=dfrPrd,mapping=aes(x={{x_}},y={{y_}},z=exp(prdLnR+(seLnR^2)/2))) +
        geom_contour_filled(breaks=c(0.01,0.1,0.5,1,2,10)) + 
        labs(x=xlab,ylab=ylab,
             subtitle=paste0(tolower(x),"s: ","predicted R")) +
        wtsPlots::getStdTheme();
  if (!is.null(xints))  p = p + geom_vline(xintercept=xints,linetype=2);
  if (!is.null(dfrDat)) {
    p = p + geom_point(data=dfrDat,mapping=aes(x={{x_}},y={{y_}}),inherit.aes=FALSE,size=0.1,colour="white")
    p = p + geom_rug(data=dfrDat,mapping=aes(x={{x_}}),inherit.aes=FALSE);
    p = p + geom_rug(data=dfrDat,mapping=aes(y={{y_}}),inherit.aes=FALSE);
  }
  ps[["R"]] = p;
  return(ps);
}

plotGAM<-function(mdl,
                  se=TRUE,
                  seWithMean=TRUE,
                  residuals=TRUE,
                  rug=TRUE,
                  theta=30,
                  phi=45,
                  all.terms=TRUE,
                  shade=TRUE,
                  scheme=2
                  ){
  res = plot(mdl,se=se,seWithMean=seWithMean,residuals=residuals,rug=rug,
             theta=theta,phi=phi,all.terms=all.terms,
             shade=shade,scheme=scheme);
  return(res);
}
