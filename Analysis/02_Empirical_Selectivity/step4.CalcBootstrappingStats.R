#--calculate statistics for empirical selectivity functions
require(mgcv);

dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

load(file=file.path(dirThs,"rda_Step3a_EmpiricalSelectivityFromBootstrapping.RData"));

#--calculate bootstrap means by sex, year, size
dfrMeans<-reshape2::dcast(dfrESs,"x+y+z~.",fun.aggregate=mean,na.rm=TRUE,value.var="emp_sel");
names(dfrMeans)[4]<-"mean";
dfrStds <-reshape2::dcast(dfrESs,"x+y+z~.",fun.aggregate=var, na.rm=TRUE,value.var="emp_sel");
dfrStds[["."]]<-sqrt(dfrStds[["."]]);
names(dfrStds)[4]<-"stddev";
cis<-wtsUtilities::calcCIs(dfrMeans$mean,sdvs=dfrStds$stddev,pdfType="normal",ci=0.8);
dfrStats<-cbind(dfrMeans,stddev=dfrStds[["stddev"]],lci=cis$lci,uci=cis$uci);
dfrStats$y<-factor(dfrStats$y)
pl<-wtsPlots::plotMDFR.XY(dfrStats,x="z",value.var="mean",plotLines=TRUE,plotPoints=FALSE,colour="y",facet_grid="x~.",showPlot=TRUE)
pl<-pl + ggplot2::geom_ribbon(mapping=aes_string(ymin="lci",ymax="uci",fill="y"),alpha=0.4,colour=NA);
pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity",colour="year",fill="year");
pl<-pl + ggplot2::geom_hline(yintercept=c(0,1),linetype=2);
pl<-pl + ggplot2::scale_y_continuous(limits=c(-1,2),oob=scales::squish);
pl<-pl+wtsPlots::getStdTheme();
print(pl);

#--define function for plotting model responses
plotPred<-function(nm,mod,zs,yfs,zfs){
  dfrPrd = tidyr::crossing(z=zs,yf=yfs,zf=zfs);
  prd<-predict(mod,dfrPrd,type="response",se.fit=TRUE);
  cis<-wtsUtilities::calcCIs(prd$fit,sdvs=prd$se.fit,pdfType="lognormal",ci=0.8);
  dfrPrdRes<-cbind(dfrPrd,fit=prd$fit,se=prd$se.fit,lci=cis$lci,uci=cis$uci);
  pl<-wtsPlots::plotMDFR.XY(dfrPrdRes,x="z",value.var="fit",colour="yf",plotLines=TRUE,plotPoints=FALSE,showPlot=FALSE)
  pl<-pl + ggplot2::geom_ribbon(mapping=aes_string(x="z",ymin="lci",ymax="uci",fill="yf"),alpha=0.4,colour=NA,inherit.aes=FALSE);
  pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity",fill="year",colour="year");
  pl<-pl + ggplot2::ggtitle(nm); 
  pl<-pl + ggplot2::geom_hline(yintercept=c(0,1),linetype=2);
  pl<-pl + wtsPlots::getStdTheme() + ggplot2::theme(panel.grid=ggplot2::element_line(colour="grey"));
  return(pl)
}

#--run models to fit smooth curves to bootstrapping results
if (FALSE){
  #----list for all models
  lstMods<-list();
  for (x_ in c("male","female")){
    #--list for models of sex x_
    lstModsX<-list();
    #--extract relevant bootstrap results
    dfrESsX = dfrESs |> 
                dplyr::filter(x==x_) |> 
                dplyr::mutate(yf = factor(y),
                              zf = factor(z),
                              nTot = n_BSFRF+n_NMFS);
    #----gaussian with log link
    nm<-paste0(x_,"s: gaussian(log link): emp_sel ~ s(z,bs='tp',id=1) + s(z,bs='re',by=yf,id=1)");
    lstModsX[[nm]]<-gam(emp_sel ~ s(z,bs="tp",id=1) + s(z,bs="re",by=yf,id=1),
                       data=dfrESsX,family=gaussian(link="log"),weights=nTot);
    #----
    nm<-paste0(x_,"s: gaussian(log link) emp_sel~ s(z,bs='tp') + s(zf,yf,bs='re')");
    lstModsX[[nm]]<-gam(emp_sel~ s(z,bs="tp") + s(zf,yf,bs="re"), 
                       data=dfrESsX,family=gaussian(link="log"),weights=nTot);
    #--
    lstMods[[x_]] = lstModsX;
  }#--x_
  save(dfrStats,lstMods,file=file.path(dirThs,"rda_Step4_EmpiricalSelectivtyFromBootstrapping.stats.RData"));
}
 if (FALSE) load(file=file.path(dirThs,"rda_Step4_EmpiricalSelectivtyFromBootstrapping.stats.RData"))
#--print/plot results for the models
zs<-list(male  =seq(27,182,by=5),
         female=seq(27,127,by=5));
for (x_ in c("female")){
  cat("\n\n");
  cat("--------------------------------------------------------------\n");
  cat("--------------------------------------------------------------\n");
  cat("results for",x_,"\n")
  #--extract model results
  lstModsX = lstMods[[x_]];
  #--extract relevant bootstrap results
  dfrESsX = dfrESs |> 
              dplyr::filter(x==x_) |> 
              dplyr::mutate(yf = factor(y),
                            zf = factor(z),
                            nTot = n_BSFRF+n_NMFS);
  for (nm in names(lstModsX)){
    cat("\n--------------------------------------------------------------\n");
    cat("results for",nm,"\n")
    summary(lstModsX[[nm]]);
    plot(lstModsX[[nm]],page=1);
    gam.check(lstModsX[[nm]]);
  }
  print(anova(lstModsX[[1]],lstModsX[[2]]));
  print(AIC(lstModsX[[1]],lstModsX[[2]]));
  nm<-names(lstModsX)[1];
  p = plotPred(nm,lstModsX[[nm]],zs=zs[[x_]],yfs=dfrESsX$yf,zfs=dfrESsX$zf);
  print(p);
  nm<-names(lstModsX)[2];
  p = plotPred(nm,lstModsX[[nm]],zs=zs[[x_]],yfs=0,zfs=0);
  print(p);
}

