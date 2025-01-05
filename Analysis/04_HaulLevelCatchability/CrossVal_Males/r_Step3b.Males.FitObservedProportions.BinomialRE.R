#--fit various models for ln(r) using mgcv to fit GAMs for males using the BINOMIAL distribution----
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

#--get censored data and prediction grids----
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Males.RData"));

#--remove zeros, infs, questionable observed Rs----
#dfrDatp   = lst$dfrDat |> dplyr::filter(obsR<10, is.finite(lnR),between(z,15,150));
dfrDatp   = lst$dfrDat |> dplyr::filter(between(z,15,150));
lvls = c("any",unique(dfrDatp$h));
dfrDatpp = dfrDatp |> dplyr::mutate(h=factor(h,levels=lvls));

#--BINOMIAL regression  models for lnR----
famB = stats::binomial(link="logit");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  # ks=c(20,10);
  # k1 = ks[1]; k2 = ks[2];
  # frmla  = p~ti(z,bs="ts",k=k1)   +
  #            ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
  #            ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
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

