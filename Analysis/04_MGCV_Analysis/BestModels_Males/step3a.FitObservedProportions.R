#--fit various models for ln(r) using mgcv to fit GAMs for a specific sex
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
source(file.path(dirThs,"../r_estSelRatio.R"));

#--get trimmed data
dfrDat =  wtsUtilities::getObj(file.path(dirThs,"../rda_dfrTrimmedData.RData"));

#--extract characteristics from size and environmental data
#----size
grd_z = seq(5,180,5)+2.5; med_z = 100.0; #--for males
#--depth
med_d = median(dfrDat$d,na.rm=TRUE); rng_d = range(dfrDat$d,na.rm=TRUE);
grd_d = seq(from=rng_d[1],rng_d[2],length.out=50);
#--temperature
med_t = median(dfrDat$t,na.rm=TRUE); rng_t = range(dfrDat$t,na.rm=TRUE);
grd_t = seq(from=rng_t[1],rng_t[2],length.out=50);
#--phi
med_f = median(dfrDat$f,na.rm=TRUE); rng_f = range(dfrDat$f,na.rm=TRUE);
grd_f = seq(from=rng_f[1],rng_f[2],length.out=50);
#--sorting
med_s = median(dfrDat$s,na.rm=TRUE); rng_s = range(dfrDat$s,na.rm=TRUE);
grd_s = seq(from=rng_s[1],rng_s[2],length.out=50);

grids = list(z=grd_z,d=grd_d,t=grd_t,f=grd_f,s=grd_s);
lbls  = list(z="size (mm CW)",d="depth (m)",t="temperature (deg C)",f="-log2(phi)",s="sorting coeff.")

#--the logit-scale observed proportions "o" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_a is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

pltctr = 1;

#--set sex
x = "MALE";
dfrDat   = dfrDat |> dplyr::filter(x=="MALE") |> dplyr::mutate(obsR=exp(lnR));
dfrDatRE = dfrDat |> dplyr::mutate(h=as.factor(h));
ggplot(dfrDat,aes(x=z,y=obsR,colour=as.factor(y))) + geom_point(alpha=0.5);

#--set up BINOMIAL models
famBM = stats::binomial(link="logit");
mdlsBM = list();
rsltBM = list();

#--set up NEGATIVE BINOMIAL models
famNB = mgcv::nb(link="log");
dfrDatNB = dfrDat |> dplyr::filter(is.finite(log(obsR)));
ggplot(dfrDatNB,aes(x=z,y=obsR,size=n)) + geom_point(position="jitter") + 
  scale_y_log10() + scale_size_area()
ggplot(dfrDatNB,aes(x=z,y=log(obsR),size=n)) + geom_point(position="jitter") + 
  scale_size_area()
mdlsNB = list();
rsltNB = list();
#---------mdlNB_1DZa: ONLY SIZE TERM; k=20----------------------------------
#--ln(r) = ti(z)
k = 20;
frmla = obsR~s(z,bs="ts",k=k);
mdlNB_1DZa = estSelRatio(dfrDatNB,formula=frmla,family=famNB);
mdl = mdlNB_1DZa$model;

summary(mdl);
mdl$family$getTheta();
gam.check(mdl)
plot(mdl,residuals=TRUE)
gratia::draw(mdl) + scale_y_log10()
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
  plot(simResids)
  p1 = DHARMa::plotQQunif(simResids);
  p2 = DHARMa::plotResiduals(simResids);
  DHARMa::testUniformity(simResids);
  DHARMa::testDispersion(simResids,alternative="less");
  DHARMa::testZeroInflation(simResids);
  DHARMa::plotResiduals(simResids,form=dfrDat$z,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$d,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$t,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$f,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$s,quantreg=TRUE);


#---------mdlB_1DZa: ONLY SIZE TERM; k=20----------------------------------
#--ln(r) = ti(z)
  k = 20;
  frmla = p~s(z,bs="ts",k=k);
  mdlB_1DZa = estSelRatio(dfrDat,formula=frmla,family=famBM);
  frmla = p~ti(z,bs="ts",k=20) + ti(z,h,bs="fs");
  mdlRE = estSelRatio(dfrDatRE |> dplyr::filter(h %in% unique(dfrDat$h)[1:20]),formula=frmla,family=famBM);
#--results
  mdl = mdlB_1DZa;
  lst_sumry   = summary(mdl$model);
  txt_sumry   = paste0(capture.output(print(lst_sumry)),collapse="\n");
  txt_gamchk  = wtsMGCV::gam.check(mdl$model);
  ps_gamchk   = wtsMGCV::gam.check.plots(mdl$model);
  lst_prdtrms = wtsMGCV::predSmoothTerms(mdl$model,grids);
  ps_prdtrms  = wtsMGCV::plotSmoothTerms(lst_prdtrms,labs=lbls,dfrDat=dfrDat,ori="V",ci=0.80);
  df_prdsel_z = calcPrdSel(mdl$model,wtsMGCV::createGridTbl(z=grd_z,d=med_d,t=med_t,f=med_f,s=med_s));
  ps_prdsel_z = plotPrdSel1D(df_prdsel_z,x_=z,xlab="size (mm CW)",xints=c(25,180),dfrDat=dfrDat,ci=0.80);
#--DHARMa diagnostics
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl$model, plot = F, n=1000);
  plot(simResids)
  p1 = DHARMa::plotQQunif(simResids);
  p2 = DHARMa::plotResiduals(simResids);
  DHARMa::testUniformity(simResids);
  DHARMa::testDispersion(simResids,alternative="less");
  DHARMa::testZeroInflation(simResids);
  DHARMa::plotResiduals(simResids,form=dfrDat$z,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$d,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$t,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$f,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$s,quantreg=TRUE);

#---------mdlB_1DZ: ONLY SIZE TERM; k = 20----------------------------------
#--ln(r) = ti(z)
  k = 20;
  frmla = p~ti(z,bs="ts",k=k);
  mdlB_1DZb = estSelRatio(dfrDat,formula=frmla,family=fam);
#--results
  mdl = mdlB_1DZb;
  lst_sumry   = summary(mdl$model);
  txt_sumry   = paste0(capture.output(print(lst_sumry)),collapse="\n");
  txt_gamchk  = wtsMGCV::gam.check(mdl$model);
  ps_gamchk   = wtsMGCV::gam.check.plots(mdl$model);
  lst_prdtrms = wtsMGCV::predSmoothTerms(mdl$model,grids);
  ps_prdtrms  = wtsMGCV::plotSmoothTerms(lst_prdtrms,labs=lbls,dfrDat=dfrDat,ori="V",ci=0.80);
  df_prdsel_z = calcPrdSel(mdl$model,wtsMGCV::createGridTbl(z=grd_z,d=med_d,t=med_t,f=med_f,s=med_s));
  ps_prdsel_z = plotPrdSel1D(df_prdsel_z,x_=z,xlab="size (mm CW)",xints=c(25,180),dfrDat=dfrDat,ci=0.80);
#--DHARMa diagnostics
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl$model, plot = F, n=1000);
  plot(simResids)
  p1 = DHARMa::plotQQunif(simResids);
  p2 = DHARMa::plotResiduals(simResids);
  DHARMa::testOutliers(simResids,type="bootstrap",nBoot=1000)
  DHARMa::testUniformity(simResids);
  DHARMa::testDispersion(simResids,alternative="less");
  DHARMa::testZeroInflation(simResids);
  DHARMa::plotResiduals(simResids,form=dfrDat$z,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$d,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$t,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$f,quantreg=TRUE);
  DHARMa::plotResiduals(simResids,form=dfrDat$s,quantreg=TRUE);

  #----SIMPLE RE(h)----------------------------
k = 20;
frmla = p~ti(z,bs="ts",k=k) + s(h,bs="re");
mdlB_1DREZb = estSelRatio(dfrDatRE,frmla,family=fam,select=TRUE);




#---------ALL 1-WAY INTERACTIONS----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s)
k = 10; 
frmla  = p~ti(z,bs="ts",k=k)+ti(d,bs="ts",k=k)+ti(t,bs="ts",k=k)+ti(f,bs="ts",k=k)+ti(s,bs="ts",k=k);
mdlB_All1Da = estSelRatio(dfrDat,frmla,family=fam,select=TRUE);
#---------ALL 1-WAY INTERACTIONS----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s)
k = 20; 
frmla  = p~ti(z,bs="ts",k=k)+ti(d,bs="ts",k=k)+ti(t,bs="ts",k=k)+ti(f,bs="ts",k=k)+ti(s,bs="ts",k=k);
mdlB_All1Db = estSelRatio(dfrDat,frmla,family=fam,select=TRUE);
#---------ALL 1-WAY INTERACTIONS----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s)
k = 30; #--PROBLEMS WITH THIS!!
frmla  = p~ti(z,bs="ts",k=k)+ti(d,bs="ts",k=k)+ti(t,bs="ts",k=k)+ti(f,bs="ts",k=k)+ti(s,bs="ts",k=k);
mdlB_All1Dc = estSelRatio(dfrDat,frmla,family=fam);
#---------ALL 1-WAY INTERACTIONS + RE-----------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s)
k = 20; 
frmla  = p~ti(z,bs="ts",k=k)+ti(d,bs="ts",k=k)+ti(t,bs="ts",k=k)+ti(f,bs="ts",k=k)+ti(s,bs="ts",k=k) + 
           s(h,bs="re");
mdlB_All1DREb = estSelRatio(dfrDatRE,frmla,family=fam,select=TRUE);

#--results
mdl = mdlB_All1Db;
lst_sumry   = summary(mdl$model);
txt_sumry   = paste0(capture.output(print(lst_sumry)),collapse="\n");
txt_gamchk  = wtsMGCV::gam.check(mdl$model);
ps_gamchk   = wtsMGCV::gam.check.plots(mdl$model);
lst_prdtrms = wtsMGCV::predSmoothTerms(mdl$model,grids);
ps_prdtrms  = wtsMGCV::plotSmoothTerms(lst_prdtrms,labs=lbls,dfrDat=dfrDat,ori="V",ci=0.80);
df_prdsel_z = calcPrdSel(mdl$model,wtsMGCV::createGridTbl(z=grd_z,d=med_d,t=med_t,f=med_f,s=med_s));
ps_prdsel_z = plotPrdSel1D(df_prdsel_z,x_=z,xlab="size (mm CW)",xints=c(25,180),dfrDat=dfrDat,ci=0.80);
#--DHARMa diagnostics
simResids <- DHARMa::simulateResiduals(fittedModel = mdl$model, plot = F, n=1000);
plot(simResids)
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);
DHARMa::testOutliers(simResids,type="bootstrap",nBoot=1000);
DHARMa::testDispersion(simResids,alternative="greater");
DHARMa::testZeroInflation(simResids);
DHARMa::plotResiduals(simResids,form=dfrDat$z,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$d,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$t,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$f,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$s,quantreg=TRUE);

#---------ALL 2-WAY INTERACTIONS----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s)
k = 20; k1 = 10; k2 = 6;
frmla  = p~ti(z,bs="ts",k=k)+ti(d,bs="ts",k=k1)+ti(t,bs="ts",k=k1)+ti(f,bs="ts",k=k1)+ti(s,bs="ts",k=k1)+
           ti(z,d,bs="ts",k=k2)+ti(z,t,bs="ts",k=k2)+ti(z,f,bs="ts",k=k2)+ti(z,s,bs="ts",k=k2)+
           ti(d,t,bs="ts",k=k2)+ti(d,f,bs="ts",k=k2)+ti(d,s,bs="ts",k=k2)+
           ti(t,f,bs="ts",k=k2)+ti(t,s,bs="ts",k=k2)+
           ti(f,s,bs="ts",k=k2);
mdlB_All2Db = estSelRatio(dfrDat,frmla,family=fam,select=TRUE);
#--results
mdl = mdlB_All2Db;
lst_sumry   = summary(mdl$model);
txt_sumry   = paste0(capture.output(print(lst_sumry)),collapse="\n");
txt_gamchk  = wtsMGCV::gam.check(mdl$model);
ps_gamchk   = wtsMGCV::gam.check.plots(mdl$model);
lst_prdtrms = wtsMGCV::predSmoothTerms(mdl$model,grids);
ps_prdtrms  = wtsMGCV::plotSmoothTerms(lst_prdtrms,labs=lbls,dfrDat=dfrDat,ori="V",ci=0.80);
df_prdsel_z = calcPrdSel(mdl$model,wtsMGCV::createGridTbl(z=grd_z,d=med_d,t=med_t,f=med_f,s=med_s));
ps_prdsel_z = plotPrdSel1D(df_prdsel_z,x_=z,xlab="size (mm CW)",xints=c(25,180),dfrDat=dfrDat,ci=0.80);
#--DHARMa diagnostics
simResids <- DHARMa::simulateResiduals(fittedModel = mdl$model, plot = F, n=1000);
plot(simResids)
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);
DHARMa::testOutliers(simResids,type="bootstrap",nBoot=1000);
DHARMa::testDispersion(simResids,alternative="greater");
DHARMa::testZeroInflation(simResids);
DHARMa::plotResiduals(simResids,form=dfrDat$z,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$d,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$t,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$f,quantreg=TRUE);
DHARMa::plotResiduals(simResids,form=dfrDat$s,quantreg=TRUE);

#---------SIGNIFICANT 2-WAY INTERACTIONS----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + MISSING[ti(f))] + ti(s) + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + ti(t,f) + ti(t,s) + ti(f,s)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+ti(s,bs="ts")+
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
           ti(d,t,bs="ts")+ti(d,f,bs="ts")+ti(d,s,bs="ts")+
           ti(t,f,bs="ts")+ti(t,s,bs="ts")+
           ti(f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=NULL);
mdlsB[["signif 2-way interactions"]] = mdl;
summary(mdl$model);
#plot(mdl$model);
rm(mdl);

#--------ALL 3-WAY INTERACTIONS--------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s) +
#--        ti(z,d,t) + (z,d,f) + (z,d,s) + (z,t,f) + (z,t,s) + (z,f,s) +
#--        ti(d,t,f) + ti(d,t,s) +
#--        ti(t,f,s)
k = 20; k1 = 10; k2 = 6;
frmla  = p~ti(z,bs="ts",k=k)   + 
           ti(d,bs="ts",k=k1)   + ti(t,bs="ts",k=k1)   + ti(f,bs="ts",k=k1)   + ti(s,bs="ts",k=k1) +
           ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
           ti(d,t,bs="ts",k=k2) + ti(d,f,bs="ts",k=k2) + ti(d,s,bs="ts",k=k2) +
           ti(t,f,bs="ts",k=k2) + ti(t,s,bs="ts",k=k2)+
           ti(f,s,bs="ts",k=k2) +
           ti(z,d,t,bs="ts",k=c(k,k2,k2)) + ti(z,d,f,bs="ts",k=c(k,k2,k2)) + ti(z,d,s,bs="ts",k=c(k,k2,k2)) +
           ti(z,t,f,bs="ts",k=c(k,k2,k2)) + ti(z,t,s,bs="ts",k=c(k,k2,k2)) +
           ti(z,f,s,bs="ts"k=c(k,k2,k2)) +
           ti(d,t,f,bs="ts",k=k2) + ti(d,t,s,bs="ts",k=k2) + ti(d,f,s,bs="ts",k=k2) + ti(t,f,s,bs="ts",k=k2);
mdlB_All3D  = estSelRatio(dfrDat,frmla,family=fam);

k = 20;
frmla = p~s(z,bs="ts",k=k)   +  s(z,h,bs="fs");
mdlB_1DHZa = estSelRatio(dfrDatRE,frmla,family=fam);

#--------SIGNIFICANT 3-WAY INTERACTIONS--------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) + MISSING[ti(z,s)] + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + MISSING[ti(t,s)] + 
#--        MISSING[ti(f,s)]
#--        ti(z,d,t) + (z,d,f) + (z,d,s) + (z,t,f) + (z,t,s) + (z,f,s) +
#--        ti(d,t,f) + ti(d,t,s) +
#--        ti(t,f,s)
frmla  = p~ti(z,bs="ts")   + ti(d,bs="ts")   + ti(f,bs="ts")   + ti(s,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + # ti(z,s,bs="ts") +
           ti(d,t,bs="ts") + ti(d,f,bs="ts") + ti(d,s,bs="ts") +
           ti(t,f,bs="ts") + # ti(t,s,bs="ts") +
           #ti(f,s,bs="ts") +
           ti(z,d,t,bs="ts") + ti(z,d,f,bs="ts") + ti(z,d,s,bs="ts") +
           ti(z,t,f,bs="ts") + ti(z,t,s,bs="ts") +
           ti(z,f,s,bs="ts") +
           ti(d,t,f,bs="ts") + ti(d,t,s,bs="ts") + ti(d,f,s,bs="ts") + ti(t,f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["signif 3-way interactions"]] = mdl;
summary(mdl$model);
#plot(mdl$model);
rm(mdl);

#--------2nd PASS SIGNIFICANT 3-WAY INTERACTIONS--------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + MISSING[ti(s)] +
#--        ti(z,d) + ti(z,t) + ti(z,f) + MISSING[ti(z,s)] + 
#--        ti(d,t) + MISSING[ti(d,f)] + MISSING[ti(d,s)] + 
#--        ti(t,f) + MISSING[ti(t,s)] + 
#--        MISSING[ti(f,s)]
frmla  = p~ti(z,bs="ts")   + ti(d,bs="ts")   + ti(f,bs="ts")   + # ti(s,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + # ti(z,s,bs="ts") +
           ti(d,t,bs="ts") + # ti(d,f,bs="ts") + ti(d,s,bs="ts") +
           ti(t,f,bs="ts") + # ti(t,s,bs="ts") +
           #ti(f,s,bs="ts") +
           ti(z,d,t,bs="ts") + ti(z,d,f,bs="ts") + ti(z,d,s,bs="ts") +
           ti(z,t,f,bs="ts") + ti(z,t,s,bs="ts") +
           ti(z,f,s,bs="ts") +
           ti(d,t,f,bs="ts") + ti(d,t,s,bs="ts") + ti(d,f,s,bs="ts") + ti(t,f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["signif 3-way interactions: 2nd pass"]] = mdl;
summary(mdl$model);
#plot(mdl$model);
rm(mdl);

#---------ALL 2-WAY INTERACTIONS + RE(h) ----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f)) + ti(s) + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + ti(t,f) + ti(t,s) + ti(f,s) +
#--        RE(h)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+ti(f,bs="ts")+ti(s,bs="ts")+
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
           ti(d,t,bs="ts")+ti(d,f,bs="ts")+ti(d,s,bs="ts")+
           ti(t,f,bs="ts")+ti(t,s,bs="ts")+
           ti(f,s,bs="ts") +
           s(h,bs="re");
#mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=NULL);
mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="ML");
summary(mdl);
#plot(mdl);
mdlsB[["signif 2-way interactions + RE"]] = mdl;
rm(mdl);

#---------1st PASS SIGNIFICANT 2-WAY INTERACTIONS + RE(h) ----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f)) + MISSING[ti(s)] + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        MISSING[ti(d,t)] + MISSING[ti(d,f)] + ti(d,s) + 
##--       MISSING[ti(t,f)] + MISSING[ti(t,s)] + MISSING[ti(f,s)] +
#--        RE(h)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+ti(f,bs="ts")+ #ti(s,bs="ts")+
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
        #   ti(d,t,bs="ts")+ti(d,f,bs="ts")+
           ti(d,s,bs="ts")+
        #   ti(t,f,bs="ts")+ti(t,s,bs="ts")+
        #   ti(f,s,bs="ts") +
           s(h,bs="re");
#mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=NULL);
mdlREML.RE1 = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="REML");
summary(mdl);
#plot(mdl);
mdlsB[["signif 2-way interactions + RE"]] = mdl;
rm(mdl);

#---------2nd PASS SIGNIFICANT 2-WAY INTERACTIONS + RE(h) ----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + MISSING[ti(f))] + MISSING[ti(s)] + 
#--        ti(z,d) + MISSING[ti(z,t)] + MISSING[ti(z,f)] + ti(z,s) + 
#--        MISSING[ti(d,t)] + MISSING[ti(d,f)] + ti(d,s) + 
##--       MISSING[ti(t,f)] + MISSING[ti(t,s)] + MISSING[ti(f,s)] +
#--        RE(h)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+#ti(f,bs="ts")+ti(s,bs="ts")+
           ti(z,d,bs="ts")+#ti(z,t,bs="ts")+ti(z,f,bs="ts")+
           ti(z,s,bs="ts")+
        #   ti(d,t,bs="ts")+ti(d,f,bs="ts")+
           ti(d,s,bs="ts")+
        #   ti(t,f,bs="ts")+ti(t,s,bs="ts")+
        #   ti(f,s,bs="ts") +
           s(h,bs="re");
#mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=NULL);
mdlREML.RE2 = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="REML");
summary(mdlREML.RE2);

stats::AIC(mdlsB[[1]]$model,mdlsB[[2]]$model,mdlsB[[3]]$model,mdlsB[[4]]$model,mdlsB[[5]]$model,mdlsB[[6]]$model);
stats::BIC(mdlsB[[1]]$model,mdlsB[[2]]$model,mdlsB[[3]]$model,mdlsB[[4]]$model,mdlsB[[5]]$model,mdlsB[[6]]$model);

wtsUtilities::saveObj(mdlsB,     paste0("selModels.Binom.All.",tolower(x),"s.RData"));
#--wtsUtilities::saveObj(mdlsB[[3]],paste0("selModels.Binom.Best.",tolower(x),"s.RData"));

# ################################################################################
# #--set up NEGATIVE BINOMIAL model
# fam = mgcv::nb(link="log");
# dfrDatNB = dfrDat |> dplyr::mutate(lgtP = log(p/(1-p)));
# mdlsNB = list();
# 
# #--ln(r) = ti(z)
# frmla = p~ti(z,bs="ts");
# mdl = estSelRatio(dfrDatNB,formula=frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
# mdlsNB[["z"]] = mdl;
# summary(mdl$model);
# plot(mdl$model);
# rm(mdl);


