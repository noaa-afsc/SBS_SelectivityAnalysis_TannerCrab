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
ggplot(dfrDat,aes(x=z,y=obsR,colour=as.factor(y))) + geom_point(alpha=0.5);

#--remove zeros, infs, questionable observed Rs
dfrDatp   = dfrDat |> dplyr::filter(obsR<10, is.finite(lnR));
dfrDatpRE = dfrDatp |> dplyr::mutate(h=as.factor(h));
ggplot(dfrDatp,aes(x=z,y=obsR,colour=as.factor(y))) + geom_point(alpha=0.5);
ggplot(dfrDatp,aes(x=z,y=obsR,colour=as.factor(y),size=n)) + geom_point(alpha=0.5) + 
  scale_y_log10() + scale_size_area();

#--set up NEGATIVE BINOMIAL models for lnR
famNB = mgcv::nb(link="log");
mdlsNB = list();
rsltNB = list();
#---------mdlNB_1DZa: ONLY SIZE TERM; k=10----------------------------------
#--ln(r) = ti(z)
k = 10;
frmla = obsR~s(z,bs="ts",k=k);
mdlNB_1DZa = estSelRatio(dfrDatp,formula=frmla,family=famNB);
mdl = mdlNB_1DZa$model;
ggplot(mdlNB_1DZa$results,aes(x=z)) + 
  geom_point(aes(y=obsR),position="jitter") + 
  geom_point(aes(y=fits),colour="green")
ggplot(mdlNB_1DZa$results,aes(x=z)) + 
  geom_point(aes(y=lnR),position="jitter") + 
  geom_point(aes(y=log(fits)),colour="green")
dfrPrd = calcPrdSel(mdl,newDFR(z=grd_z),type="link")
ggplot(mdlNB_1DZa$results,aes(x=z)) + 
  geom_point(aes(y=lnR),position="jitter") + 
  geom_line(aes(y=prd),data=dfrPrd,colour="green") + 
  geom_ribbon(aes(ymin=prd-se,ymax=prd+se),data=dfrPrd,colour="green",alpha=0.5)

gam.check(mdl)
plot(mdl);
plot(mdl,residuals=TRUE)
gratia::draw(mdl,residuals=TRUE)
simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);
dfrPrd = calcPrdSel(mdl,newDFR(z=grd_z),type="link")
ggplot() + 
  geom_point(aes(x=z,y=lnR), data=dfrDatp,position="jitter") + 
  geom_line(aes(x=z,y=prdLnR),data=dfrPrd) + 
  geom_ribbon(aes(x=z,ymin=prdLnR-seLnR,ymax=prdLnR+seLnR),data=dfrPrd,alpha=0.5)
ggplot() + 
  geom_point(aes(x=z,y=obsR),data=dfrDatp,position="jitter") + 
  geom_line(aes(x=z,y=prdR), data=dfrPrd) + 
  geom_ribbon(aes(x=z,ymin=exp(prdLnR-seLnR),ymax=exp(prdLnR+seLnR)),data=dfrPrd,alpha=0.5) + 
  coord_cartesian(ylim=c(0,2))



