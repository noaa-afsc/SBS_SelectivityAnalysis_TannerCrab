#--fit random effects model for logit-scale proportions p using mgcv to fit GAMs for males using the BINOMIAL distribution----
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

#--get censored data and prediction grids----
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04_HaulLevelCatchability/CrossVal_Females");
setwd(dirThs);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Females.RData"));

#--remove zeros, infs, questionable observed Rs----
dfrDatp   = lst$dfrDat |> dplyr::filter(between(z,15,130));
lvls = c("any",unique(dfrDatp$h));
dfrDatpp = dfrDatp |> dplyr::mutate(h=factor(h,levels=lvls));

#--BINOMIAL regression  models for logit-scale proportions p----
famB = stats::binomial(link="logit");
#--------MAIN SMOOTH(Z) + factor smooth REs by haul--------------------------
  #--lgtp = s(z,bs="ts",k=k1) + ti(z,h,bs="fs",k=k1)
  # ks=c(10,8);
  # k1 = ks[1]; k2 = ks[2];
  # frmla  = p = s(z,bs="ts",k=k1) + ti(z,h,bs="fs",k=k1)
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = p~s(z,bs="ts",k=k1) + ti(z,h,bs="fs",k=k1);

if (FALSE){
  mdl_ZE2D  = mgcv::gam(frmla,family=famB,data=dfrDatpp,select=FALSE,method="ML",fit=TRUE,
                        offset=lnq,weights=n,drop.unused.levels=FALSE);
  wtsUtilities::saveObj(mdl_ZE2D,file.path(dirThs,"rda_Step3b1.BinomialModels_RE.RData"));
}

#--make plots, etc.
if (FALSE){
  source(file.path(dirThs,"..","r_PredictionsAndPlots.R"));
  mdl = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b1.BinomialModels_RE.RData"));
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);

  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  dfrPrd = prdMod(mdl,trms=c("all"),type="response",lst=grdPrd,p=0.10);
  plotMod(dfrPrd)
  
  #wtsUtilities::saveObj(??,file.path(dirThs,"rda_Step3b3.BionomialModels_RE.RData"));
}

