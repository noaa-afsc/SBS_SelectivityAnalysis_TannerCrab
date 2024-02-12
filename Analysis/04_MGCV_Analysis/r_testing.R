gam.check(mdlsB[["z"]]$model)
gam.check(mdlsB[["all 2-way interactions"]]$model)

estSelRatio<-function(
            dfrDat,
            formula,
            family=stats::binomial(),
            select=FALSE,
            scale=0){
  mdl = gam(data=dfrDat,family=family,formula=formula,weights=n,select=select,scale=scale,offset=lnq);
  
  #--extract fitted values and calculate associated selectivity ratio
  if (family$family=="binomial"){
    dfrDat %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., prNMFS)
                              residuals=residuals(mdl,type="deviance"), #--deviance residuals for qqplot
                              lnR=log(fits/(1-fits))-lnq,               #--corresponding values of lnR
                              R=exp(lnR));                              #--corresponding values of R
  } else { #--negative binomial
    dfrDat %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., logit(prNMFS))
                              residuals=residuals(mdl,type="deviance"), #--deviance residuals for qqplot
                              lnR=fits-lnq,                             #--corresponding values of lnR
                              R=exp(lnR));                              #--corresponding values of R
  }
  return(list(model=mdl,data=dfrDat));
}

estSelRatio<-function(
            dfrDat,
            formula,
            family=stats::binomial(),
            select=FALSE,
            scale=0,
            dfrNew=NULL,
            seFactor=2,
            showPlots=FALSE){
  mdl = gam(data=dfrDat,family=family,formula=formula,weights=n,select=select,scale=scale,offset=lnq);
  
  #--extract fitted values and calculate associated selectivity ratio
  if (family$family=="binomial"){
    dfrDat %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., prNMFS)
                              residuals=residuals(mdl,type="deviance"), #--deviance residuals for qqplot
                              lnR=log(fits/(1-fits))-lnq,               #--corresponding values of lnR
                              R=exp(lnR));                              #--corresponding values of R
  } else { #--negative binomial
    dfrDat %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., logit(prNMFS))
                              residuals=residuals(mdl,type="deviance"), #--deviance residuals for qqplot
                              lnR=fits-lnq,                             #--corresponding values of lnR
                              R=exp(lnR));                              #--corresponding values of R
  }
  #--predict new values
  if (!is.null(dfrNew)){
    #--return predictions on logit-scale w/ no offset
    prd = predict(mdl,newdata=dfrNew,se.fit=TRUE,type="link",newdata.guaranteed=TRUE);
    dfrNew %<>% dplyr::mutate(prdLnR=prd$fit,  #--offset is lnq is 0, identically
                              seLnR=prd$se.fit,#--standard error
                              prdR=exp(prdLnR),#--predicted selectivity ratio
                              lower=exp(prdLnR-seFactor*seLnR),
                              upper=exp(prdLnR+seFactor*seLnR));
  }
  
  #--make plots
  psEsts = plotEsts(dfrDat,showPlots=showPlots);
  psPrds= NULL;
  if (!is.null(dfrNew)) psPrds=plotPrdSel(dfrNew,seFactor,showPlots=showPlots)
  return(list(model=mdl,data=dfrDat,new=dfrNew,plotEsts=psEsts));
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
                labs(x="size (mm CW)",y="ln(selectivity ratio)",
                     subtitle=paste0(tolower(x),"s: ","predicted ln-scale R with +/- 1 SE."));
  if (showPlots) print(ps[["3a"]]);
  ps[["3b"]] = ggplot(data=dfrNew,mapping=aes(x=z,y=prdR,ymin=lower,ymax=upper)) + 
                geom_line() + geom_ribbon(alpha=0.4) + geom_vline(xintercept=c(25.0,180.0),linetype=2) +
                labs(x="size (mm CW)",y="selectivity ratio",
                     subtitle=paste0(tolower(x),"s: ","predicted R with +/- ",seFactor," SE."))
  if (showPlots) print(ps[["3b"]]);
  
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
