#----functions to predict and plot values based on a model
#' 
#' @title Predict values based on a model 
#' @description Function to predict values based on a model. 
#' @param mdl - fitted gam model
#' @param trms - terms to include in prediction ("all" includes all terms, including the intercept)
#' @param lst - list of covariate variables to include in grid(or dataframe)
#' @param type - type of prediction (see [mgcv::predict.gam()] for options)
#' @param keep - if not NULL (default), a character vector of covariates to keep for plots 
#' @param p - value indicating size for two-sided confidence intervals (ci = 100*(1-2*p)) 
#' @return dataframe with column names corresponding to covariates, `emp_sel`, `lci`, `uci`, and `trms`.
#' @import dplyr 
#' @import wtsMGCV
#' @import tibble 
#' @import tidyselect 
#' @export
#' 
prdMod<-function(mdl,trms,lst,type="link",keep=NULL,p=0.05){
  if (!inherits(lst,"data.frame")){
    dfr = wtsMGCV::createGridTbl(lst);
  } else {
    dfr = lst;
  }
  if (any(trms=="all")){
    #--add intercept and all smooth terms
    trmsp = "(Intercept)";
    trms = wtsMGCV::getSmoothTerms(mdl);
    for (trm in trms) trmsp = c(trmsp,trm);
    trms = trmsp;
  }
  prd = dplyr::bind_cols(
            dfr,
            tibble::as_tibble(
              mgcv::predict.gam(mdl,dfr,type=type,terms=trms,se.fit=TRUE),
            ) |> 
            dplyr::mutate(type="fit",
                          lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                          uci=qnorm(p,fit,se.fit,lower.tail=FALSE),
                          terms=paste(trms,collapse=" + "))
        ) |> 
          dplyr::rename(emp_sel=fit);
  if (!is.null(keep)){
    prd = prd |> 
            dplyr::distinct(pick(tidyselect::any_of(keep),emp_sel,se.fit,lci,uci,terms));
    drop = names(lst)[!(names(lst) %in% keep)];
    for (drp in drop) prd[[drp]] = NA;
  }
  return(prd);
}
#' 
#' @title Plot predicted values by size, colored by model 
#' @description Function to plot predicted values by size, colored by model. 
#' @param ylims - y-axis limits 
#' @return ggplot2 plot object.
#' @import ggplot2 
#' @export
#' 
plotMod<-function(tmp,ylims=c(0,1.5)){
  if (("y" %in% names(tmp))&&all(is.na(tmp$y))) tmp$y = "all";
  p = ggplot(tmp,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=terms,fill=terms));
  if ("n" %in% names(tmp)){
    p = p + geom_point(aes(size=n)) + scale_size_area() + 
            geom_line();
  }
  p = p + 
         geom_ribbon(alpha=0.3) + 
         geom_line() + 
         geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
         scale_y_continuous(limits=ylims,oob=scales::squish) + 
         labs(x="size (mm CW)",y="estimated selectivity",
              colour="terms",fill="terms",size="crab\nsampled") + 
         wtsPlots::getStdTheme() + 
         theme(legend.position="inside",
               legend.position.inside=c(0.01,0.99),
               legend.justification.inside=c(0,1),
               legend.byrow=TRUE,
               legend.box="horizontal");
  return(p);
}

if (FALSE){
  #--example plot
  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  dfrPrd = prdMod(mdl_bestRE,trms=c("all"),type="response",lst=grdPrd);
  plotMod(dfrPrd)
}
