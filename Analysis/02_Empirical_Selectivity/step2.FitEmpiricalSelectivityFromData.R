#--calculate statistics for empirical selectivity functions
library(dplyr)
library(ggplot2)
library(magrittr)
library(mgcv);

dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step1_EmpiricalSelectivityFromData.RData"));

#--fit with smooth models
#----function to do model diagnostics
doDiagnostics<-function(mod){
  summary(mod);
  plot(mod,page=1);
  gam.check(mod);
}

#----function to calculate estimates from model
calcEstimates<-function(mod,zs,yfs){
  dfrPrd<-data.frame(z=rep(zs,length(yfs)),yf=rep(yfs,each=length(zs)));
  prd<-predict(mod,dfrPrd,type="response",se.fit=TRUE);
  cis<-wtsUtilities::calcCIs(prd$fit,sdvs=prd$se.fit,pdfType="lognormal",ci=0.8);
  dfrPrdRes<-cbind(dfrPrd,fit=prd$fit,se=prd$se.fit,lci=cis$lci,uci=cis$uci);
  dfrPrdRes$yf<-factor(dfrPrdRes$yf);
  return(dfrPrdRes);
}

#----function to plot estimates from model
plotEstimates<-function(dfrPrdSel,nm){
  pl<-wtsPlots::plotMDFR.XY(dfrPrdSel,x="z",value.var="fit",colour="yf",facet_grid="yf~.",ylim=c(0,1.2),plotLines=TRUE,plotPoints=FALSE,showPlot=FALSE)
  pl<-pl + ggplot2::geom_ribbon(mapping=ggplot2::aes_string(ymin="lci",ymax="uci"),alpha=0.4,fill="grey");
  pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity");
  pl<-pl + ggplot2::ggtitle(nm);
  return(pl);
}

#--fit models to "observed" catchability
#----list for models
lstMods<-list();
#----fit models by sex
sexes<-c("male","female");
mxZs <-c(180,   120);     #--need one for each sex
nXs<-length(sexes);
for (x_ in 1:nXs){
  sex<-sexes[x_]; mxZ<-mxZs[x_];
  #----emp_avail: gaussian with log link + annual fits
  dfrESsp = lst$dfrESs %>% subset((x==sexes[x_])&(z<=mxZs[x_])) %>%
             dplyr::mutate(yf=factor(y),zf=factor(z),p=NMFS/(NMFS+BSFRF)) %>%
             dplyr::group_by(type,y,x) %>%
             dplyr::mutate(w=(n_BSFRF + n_NMFS)/mean(n_BSFRF + n_NMFS)) %>%
             dplyr::ungroup();
  nm<-paste0(sex,"s: annual fits");
  lstMods[[nm]]<-gam(emp_sel~ yf+s(z,bs="cs",by=yf), weights=w, data=dfrESsp,family=gaussian(link="log"));
  #----emp_avail: gaussian with log link + aggregated fit
  dfrESsp = lst$dfrESs %>% subset((x==sexes[x_])&(z<=mxZs[x_])) %>%
             dplyr::mutate(zf=factor(z),p=NMFS/(NMFS+BSFRF)) %>%
             dplyr::group_by(type,x) %>%
             dplyr::mutate(w=(n_BSFRF + n_NMFS)/mean(n_BSFRF + n_NMFS)) %>%
             dplyr::ungroup();
  nm<-paste0(sex,"s: aggregated fit");
  lstMods[[nm]]<-gam(emp_sel~ s(z,bs="cs"), weights=w, data=dfrESsp,family=gaussian(link="log"));
}
rm(x_,sex,dfrESsp,nm);

#--get model estimate of selectivity
yfs<-c(2013:2018);
dfrModEsts<-NULL;
for (x_ in 1:nXs){
  sex = sexes[x_];
  zs<-seq(25,mxZs[x_],by=5);
  #--annual model
  nm1 <-paste0(sexes[x_],"s: annual fits")
  mdl1 = lstMods[[nm1]];
  tmp<-calcEstimates(mdl1,zs,yfs);
  tmp$model_name<-nm1; 
  plt = plotEstimates(tmp,nm1);
  data=mdl1$model %>% mutate(w=`(weights)`);
  p1 = plt + ggplot2::geom_point(data=data,mapping=ggplot2::aes_string(x="z",y="emp_sel",size="w")) +
             ggplot2::labs(colour="year",size="relative\nweight",title=NULL);
  tmp$z<-tmp$z+2;#--match z's in assessment
  dfrModEsts<-rbind(dfrModEsts,tmp);
  
  #--aggregated model
  nm2 <-paste0(sexes[x_],"s: aggregated fit")
  mdl2 = lstMods[[nm2]];
  tmp<-calcEstimates(mdl2,zs," ");
  tmp$model_name<-nm2;
  plt = plotEstimates(tmp,nm2);
  data=mdl2$model %>% mutate(w=`(weights)`,yf=" ");
  p2 = plt + ggplot2::geom_point(data=data,mapping=ggplot2::aes_string(x="z",y="emp_sel",size="w")) +
             ggplot2::labs(colour="empirical",size="relative\nweight",title = NULL);
  tmp$z<-tmp$z+2;#--match z's in assessment
  dfrModEsts<-rbind(dfrModEsts,tmp);
  
  pg = ggpubr::ggarrange(p1,p2,ncol=1,heights=c(4,2));
  print(pg);
  #--ggplot2::ggsave(paste0("EstimatedEmpiricalSelectivityFromData.",sex,"s.pdf"),pg,width=6.5,height=8);
  
  modDiags = cbind(model=c(nm1,nm2),AIC(mdl1,mdl2),BIC(mdl1,mdl2)[,2],
                   "% deviance explained"=100*c(1-mdl1$deviance/mdl1$null.deviance,
                                                1-mdl2$deviance/mdl2$null.deviance));
  #--readr::write_csv(modDiags,file = paste0("ModelDiagnostics.",sex,".csv"));
  #rm(sex,zs,nm1,mdl1,tmp,plt,nm2,mdl2,pg)
}
wtsUtilities::saveObj(dfrModEsts,"rda_Step2_EmpSel_dfrModEsts.RData");
# dfrModEstsWide<-reshape2::dcast(dfrModEsts,model_name+yf~z,fun.aggregate=wtsUtilities::Sum,value.var="fit");
# dfrModEstsLong<-reshape2::melt(dfrModEsts,id.vars=c("model_name","yf","z"),measure.vars=c("fit","se","lci","uci"));
# write.csv(dfrModEstsWide,file="EstimatedEmpiricalSelectivityFromData.Wide.csv");
# write.csv(dfrModEstsLong,file="EstimatedEmpiricalSelectivityFromData.Long.csv");
# save(lstMods,dfrModEsts,file="EstimatedEmpiricalSelectivityFromData.RData");




