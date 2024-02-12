#--calculate statistics for empirical selectivity functions
load(file="EmpiricalSelectivityFromBootstrapping.RData");

#--calculate bootstrap means by sex, year, size
dfrMeans<-reshape2::dcast(dfrESs,"x+y+z~.",fun.aggregate=mean,na.rm=TRUE,value.var="emp_sel");
names(dfrMeans)[4]<-"mean";
dfrStds <-reshape2::dcast(dfrESs,"x+y+z~.",fun.aggregate=var, na.rm=TRUE,value.var="emp_sel");
dfrStds[["."]]<-sqrt(dfrStds[["."]]);
names(dfrStds)[4]<-"stddev";
cis<-wtsUtilities::calcCIs(dfrMeans$mean,sdvs=dfrStds$stddev,pdfType="normal",ci=0.8);
dfrStats<-cbind(dfrMeans,stddev=dfrStds[["stddev"]],lci=cis$lci,uci=cis$uci);
dfrStats$y<-factor(dfrStats$y)
pl<-plotMDFR.XY(dfrStats,x="z",value.var="mean",plotLines=TRUE,plotPoints=FALSE,colour="y",facet_grid="x~.",showPlot=TRUE)
pl<-pl + ggplot2::geom_ribbon(mapping=aes_string(ymin="lci",ymax="uci",fill="y"),alpha=0.4);
pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity",colour="year",fill="year");
print(pl);

#--
dfrESsM<-dfrESs[dfrESs$x=="male",];
dfrESsM$yf<-factor(dfrESsM$y);
dfrESsM$zf<-factor(dfrESsM$z);
#----list for models
lstMods<-list();
#----gaussian with log link
nm<-"gaussian(log)";
lstMods[[nm]]<-gam(emp_sel ~ s(z,bs="tp",id=1) + s(z,bs="re",by=yf,id=1),data=dfrESsM,family=gaussian(link="log"));
#----
nm<-"gaussian(log) re 1";
lstMods[[nm]]<-gam(emp_sel~ s(z,bs="tp") + s(zf,yf,bs="re"), data=dfrESsM,family=gaussian(link="log"));

summary(lstMods[[nm]]);
plot(lstMods[[nm]],page=1);
gam.check(lstMods[[nm]]);

zs<-seq(27,182,by=5);
dfrPrd<-rbind(data.frame(z=zs,yf=0,zf=0));
prd<-predict(lstMods[[nm]],dfrPrd,type="response",se.fit=TRUE);
cis<-wtsUtilities::calcCIs(prd$fit,sdvs=prd$se.fit,pdfType="lognormal",ci=0.8);
dfrPrdRes<-cbind(dfrPrd,fit=prd$fit,se=prd$se.fit,lci=cis$lci,uci=cis$uci);
dfrPrdRes$yf<-factor(dfrPrdRes$yf);
pdf(file=paste0("EmpiricalSelectivity.",nm,".pdf"),width=8,height=6);
pl<-plotMDFR.XY(dfrPrdRes,x="z",value.var="fit",plotLines=TRUE,plotPoints=FALSE,showPlot=TRUE)
pl<-pl + ggplot2::geom_ribbon(mapping=aes_string(ymin="lci",ymax="uci"),alpha=0.4,fill="grey");
pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity");
pl<-pl + ggplot2::ggtitle(nm);
print(pl);
dev.off();

save(dfrStats,lstMods,file="EmpiricalSelectivtyFromBootstrapping.stats.RData");
