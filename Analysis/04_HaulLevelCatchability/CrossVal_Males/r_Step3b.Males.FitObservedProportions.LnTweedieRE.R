#--fit various models for ln(r) using mgcv to fit GAMs for males using the TWEEDIE distribution----
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

#--get censored data and prediction grids----
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Males.RData"));

#--remove zeros, infs, questionable observed Rs----
dfrDatp   = lst$dfrDat |> dplyr::filter(obsR<10, between(z,15,150));
lvls = c("any",unique(dfrDatp$h));
dfrDatpp = dfrDatp |> dplyr::mutate(h=factor(h,levels=lvls));

#--TWEEDIE (using mgcv::tw) regression  models for lnR----
famTW = mgcv::tw(link="log");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  # ks=c(20,10);
  # k1 = ks[1]; k2 = ks[2];
  # frmla  = obsR~ti(z,bs="ts",k=k1)   +
  #                ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
  #                ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = obsR~s(z,bs="ts",k=k1) + ti(z,h,bs="fs",k=k1);
if (FALSE){
  mdl_ZE2D  = mgcv::gam(frmla,family=famTW,data=dfrDatpp,select=FALSE,method="ML",fit=TRUE,
                        drop.unused.levels=FALSE);
  wtsUtilities::saveObj(mdl_ZE2D,file.path(dirThs,"rda_Step3b1.LnTweedieModels_RE.RData"));
}
  
#--make plots, etc.
if (FALSE){
  source(file.path(dirThs,"..","r_PredictionsAndPlots.R"));
  mdl = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b1.LnTweedieModels_RE.RData"));
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);

  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  dfrPrd = prdMod(mdl,trms=c("all"),type="response",lst=grdPrd,p=0.10);
  plotMod(dfrPrd)
  
  #wtsUtilities::saveObj(??,file.path(dirThs,"rda_Step3b3.LnTweedieModels_RE.RData"));
}

