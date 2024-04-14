#--fit various models for ln(r) using mgcv to fit GAMs for a specific sex
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
source(file.path(dirThs,"../r_estSelRatio.R"));

#--get trimmed data
dfrDat =  wtsUtilities::getObj(file.path(dirThs,"../rda_dfrTrimmedData.RData")) |> 
            dplyr::mutate(obsR=exp(lnR),
                          nN=round(p*n), #--number caught in NMFS gear
                          nB=n-nN);      #--number caught in BSFRF gear

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

#--the logit-scale observed proportions "o=nN/(nN+nB)" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_n is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

pltctr = 1;

#--set sex
x = "MALE";
dfrDat   = dfrDat |> dplyr::filter(x=="MALE") |> dplyr::mutate(obsR=exp(lnR));
ggplot(dfrDat,aes(x=z,y=p,colour=as.factor(y))) + geom_point(alpha=0.5) + 
  geom_smooth() + 
  labs(x="size (mm CW)",y="observed proportion NMFS",
       colour="survey\nyear",size="total\nnumber") + 
  wtsPlots::getStdTheme();
ggplot(dfrDat,aes(x=z,y=p)) + geom_point(alpha=0.5) + 
  geom_smooth() + 
  labs(x="size (mm CW)",y="observed proportion NMFS",
       colour="survey\nyear",size="total\nnumber") + 
  wtsPlots::getStdTheme();

#--remove zeros, infs, questionable observed Rs
dfrDatp   = dfrDat |> dplyr::filter(obsR<10, is.finite(lnR));
dfrDatpRE = dfrDatp |> dplyr::mutate(h=as.factor(h));
ggplot(dfrDatp,aes(x=z,y=p,colour=as.factor(y))) + geom_point(alpha=0.5,position="jitter") + 
  geom_smooth(aes(colour=NA)) + 
  labs(x="size (mm CW)",y="observed proportion NMFS",
       colour="survey\nyear",size="total\nnumber") + 
  wtsPlots::getStdTheme();
ggplot(dfrDatp,aes(x=z,y=lgtp)) + geom_point(alpha=0.5,position="jitter") + 
  geom_smooth() + 
  labs(x="size (mm CW)",y="logit(p)",
       colour="survey\nyear",size="total\nnumber") + 
  wtsPlots::getStdTheme();
ggplot(dfrDatp,aes(x=z,y=lnR)) + geom_point(alpha=0.5,position="jitter") + 
  geom_smooth() + 
  labs(x="size (mm CW)",y="observed selectivity ratio (R)",
       colour="survey\nyear",size="total\nnumber") + 
  wtsPlots::getStdTheme();

#--set up BETA regression models for lnR
famBR = mgcv::betar(link="logit");
mdlsBR = list();
rsltBR = list();
#---------ONLY SIZE TERM; k=10----------------------------------
if (FALSE){
  #--ln(r) = ti(z)
  k = 10;
  frmla = p~s(z,bs="ts",k=k);
  mdlBR_1DZa = estSelRatio(dfrDatp,formula=frmla,family=famBR);
  mdl = mgcv::gam(data=dfrDatp,family=family,formula=formula,weights=NULL,
                method=method,select=select,scale=scale,offset=lnq);

  wtsUtilities::saveObj(mdlBR_1DZa,file.path(dirThs,"rda_mdlBR_1DZa.RData"))
}
if(!exists(mdlBR_1DZa)) mdlBR_1DZa = wtsUtilities::getObj(file.path(dirThs,"rda_mdlBR_1DZa.RData"))
mdl = mdlBR_1DZa$model;
ggplot(mdlBR_1DZa$results,aes(x=z)) + 
  geom_point(aes(y=p),position="jitter") + 
  geom_point(aes(y=fits),colour="green")
ggplot(mdlBR_1DZa$results,aes(x=z)) + 
  geom_point(aes(y=obsR),position="jitter") + 
  geom_point(aes(y=exp(fits-lnq)),colour="green")
dfrPrd = calcPrdSel(mdl,newDFR(z=grd_z),type="link")
ggplot(mdlBR_1DZa$results,aes(x=z)) + 
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
  geom_line(aes(x=z,y=prd),data=dfrPrd) + 
  geom_ribbon(aes(x=z,ymin=prd-se,ymax=prd+se),data=dfrPrd,alpha=0.5)
ggplot() + 
  geom_point(aes(x=z,y=obsR),data=dfrDatp,position="jitter") + 
  geom_line(aes(x=z,y=exp(prd)), data=dfrPrd) + 
  geom_ribbon(aes(x=z,ymin=exp(prd-se),ymax=exp(prd+se)),data=dfrPrd,alpha=0.5) + 
  coord_cartesian(ylim=c(0,2))


#--------ALL Z/E 2-WAY INTERACTIONS--------------------------
IF (FALSE){
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  k = 10; k1 = 6; 
  frmla  = p~ti(z,bs="ts",k=k) + 
             ti(d,bs="ts",k=k1) + ti(t,bs="ts",k=k1) + ti(f,bs="ts",k=k1) + ti(s,bs="ts",k=k1) +
             ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1));
  mdlBR_ZE2D  = estSelRatio(dfrDatp,frmla,family=famBR);
  wtsUtilities::saveObj(mdlBR_ZE2D,file.path(dirThs,"rda_mdlBR_ZE2D.RData"))
}
if(!exists(mdlBR_ZE2D)) mdlBR_ZE2D = wtsUtilities::getObj(file.path(dirThs,"rda_mdlBR_ZE2D.RData"))
mdl = mdlBR_ZE2D$model;
gam.check(mdl)
plot(mdl);
plot(mdl,residuals=TRUE)
gratia::draw(mdl,residuals=TRUE)
simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);

AIC(mdlBR_1DZa$model,mdlBR_ZE2D$model)

#--------Significant Z/E 2-WAY INTERACTIONS--------------------------
IF (FALSE){
  #--ln(r) = ti(z) + 
 #--         ------- + ti(t)  +  ti(f) + ---- +
  #--        ti(z,d) + ti(z,t) + ----- + ti(z,s)
  k = 10; k1 = 6; 
  frmla  = p~ti(z,bs="ts",k=k) + 
             ti(t,bs="ts",k=k1) + 
             ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1));
  mdlBR_ZE2Dsig  = estSelRatio(dfrDatp,frmla,family=famBR);
  wtsUtilities::saveObj(mdlBR_ZE2Dsig,file.path(dirThs,"rda_mdlBR_ZE2Dsig.RData"))
}
if(!exists(mdlBR_ZE2Dsig)) mdlBR_ZE2Dsig = wtsUtilities::getObj(file.path(dirThs,"rda_mdlBR_ZE2Dsig.RData"))
mdl = mdlBR_ZE2Dsig$model;
summary(mdl)
gam.check(mdl)
plot(mdl);
plot(mdl,residuals=TRUE)
gratia::draw(mdl,residuals=TRUE)
simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);

AIC(mdlBR_1DZa$model,mdlBR_ZE2D$model,mdlBR_ZE2Dsig$model)

#--------ALL 2-WAY INTERACTIONS--------------------------
IF (FALSE){
  #--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
  #--        ti(d,t) + ti(d,f) + ti(d,s) + 
  #--        ti(t,f) + ti(t,s) + 
  #--        ti(f,s)
  k = 10; k1 = 6; k2 = 6;
  frmla  = p~ti(z,bs="ts",k=k)   + 
             ti(d,bs="ts",k=k1)   + ti(t,bs="ts",k=k1)   + ti(f,bs="ts",k=k1)   + ti(s,bs="ts",k=k1) +
             ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
             ti(d,t,bs="ts",k=k2) + ti(d,f,bs="ts",k=k2) + ti(d,s,bs="ts",k=k2) +
             ti(t,f,bs="ts",k=k2) + ti(t,s,bs="ts",k=k2)+
             ti(f,s,bs="ts",k=k2);
  mdlBR_All2D  = estSelRatio(dfrDatp,frmla,family=famBR);
  wtsUtilities::saveObj(mdlBR_All2D,file.path(dirThs,"rda_mdlBR_All2D.RData"))
}
if(!exists(mdlBR_All2D)) mdlBR_All2D = wtsUtilities::getObj(file.path(dirThs,"rda_mdlBR_All2D.RData"))
mdl = mdlBR_All2D$model;
summary(mdl)
gam.check(mdl)
plot(mdl);
plot(mdl,residuals=TRUE)
gratia::draw(mdl,residuals=TRUE)
simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);
AIC(mdlBR_1DZa$model,mdlBR_ZE2D$model,mdlBR_ZE2Dsig$model,mdlBR_All2D$model)

#--------ALL Significant 2-WAY INTERACTIONS--------------------------
IF (FALSE){
  #--ln(r) = ti(z) + --(d) + ti(t) + --(f) + --(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
  #--        --(d,t) + ti(d,f) + --(d,s) + 
  #--        ti(t,f) + --(t,s) + 
  #--        ti(f,s)
  k = 10; k1 = 6; k2 = 6;
  frmla  = p~ti(z,bs="ts",k=k)   + 
             ti(t,bs="ts",k=k1)   + 
             ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
             ti(d,f,bs="ts",k=k2) + 
             ti(t,f,bs="ts",k=k2) + 
             ti(f,s,bs="ts",k=k2);
  mdlBR_All2Dsig  = estSelRatio(dfrDatp,frmla,family=famBR);
  wtsUtilities::saveObj(mdlBR_All2Dsig,file.path(dirThs,"rda_mdlBR_All2Dsig.RData"))
}
if(!exists(mdlBR_All2Dsig)) mdlBR_All2Dsig = wtsUtilities::getObj(file.path(dirThs,"rda_mdlBR_All2Dsig.RData"))
mdl = mdlBR_All2Dsig$model;
summary(mdl);
gam.check(mdl)
plot(mdl);
plot(mdl,residuals=TRUE)
gratia::draw(mdl,residuals=TRUE)
simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
DHARMa::plotQQunif(simResids);
DHARMa::plotResiduals(simResids);
AIC(mdlBR_1DZa$model,
    mdlBR_ZE2D$model,mdlBR_ZE2Dsig$model,
    mdlBR_All2D$model,mdlBR_All2Dsig$model)



#--------ALL 3-WAY INTERACTIONS--------------------------
#--NOTE: this model did not converge after ~18 hours and was subsequently terminated
IF (FALSE){
  #--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
  #--        ti(d,t) + ti(d,f) + ti(d,s) + 
  #--        ti(t,f) + ti(t,s) + 
  #--        ti(f,s) +
  #--        ti(z,d,t) + (z,d,f) + (z,d,s) + (z,t,f) + (z,t,s) + (z,f,s) +
  #--        ti(d,t,f) + ti(d,t,s) +
  #--        ti(t,f,s)
  k = 10; k1 = 10; k2 = 6;
  frmla  = p~ti(z,bs="ts",k=k)   + 
             ti(d,bs="ts",k=k1)   + ti(t,bs="ts",k=k1)   + ti(f,bs="ts",k=k1)   + ti(s,bs="ts",k=k1) +
             ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
             ti(d,t,bs="ts",k=k2) + ti(d,f,bs="ts",k=k2) + ti(d,s,bs="ts",k=k2) +
             ti(t,f,bs="ts",k=k2) + ti(t,s,bs="ts",k=k2)+
             ti(f,s,bs="ts",k=k2) +
             ti(z,d,t,bs="ts",k=c(k,k2,k2)) + ti(z,d,f,bs="ts",k=c(k,k2,k2)) + ti(z,d,s,bs="ts",k=c(k,k2,k2)) +
             ti(z,t,f,bs="ts",k=c(k,k2,k2)) + ti(z,t,s,bs="ts",k=c(k,k2,k2)) +
             ti(z,f,s,bs="ts",k=c(k,k2,k2)) +
             ti(d,t,f,bs="ts",k=k2) + ti(d,t,s,bs="ts",k=k2) + ti(d,f,s,bs="ts",k=k2) + ti(t,f,s,bs="ts",k=k2);
  mdlBR_All3D  = estSelRatio(dfrDatp,frmla,family=famBR);
}

