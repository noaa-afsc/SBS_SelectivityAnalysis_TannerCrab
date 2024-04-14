#--fit various models for ln(r) using mgcv to fit GAMs for a specific sex
#----NORMAL distribution: FEMALES
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
#source(file.path(dirThs,"../r_estSelRatio.R"));

#--get trimmed data
dfrDat =  wtsUtilities::getObj(file.path(dirThs,"../rda_dfrTrimmedData.RData"));

#--extract characteristics from size and environmental data
#----size
grd_z = seq(5,140,5)+2.5; med_z = 70.0; #--for females
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
# where q = expF_b/expF_n is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

#--set sex
x = "FEMALE";
dfrDat   = dfrDat |> dplyr::filter(x=="FEMALE") |> dplyr::mutate(obsR=exp(lnR));
ggplot(dfrDat,aes(x=z,y=lnR,size=n,group=z)) + geom_boxplot() + geom_point()

#--remove zeros, infs, questionable observed Rs
dfrDatp   = dfrDat |> dplyr::filter(obsR<10, is.finite(lnR),between(z,15,130));
ggplot(dfrDatp,aes(x=z,y=lnR,size=n,group=z)) + geom_boxplot() + geom_point()
#--NORMAL models for lnR----------------------------------------
famNrml = stats::gaussian(link="identity");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
  #          ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) 
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = lnR~ti(z,bs="ts",k=k1)   +
             ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
             ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  mdlNrml_ZE2D  = mgcv::gam(frmla,family=famNrml,data=dfrDatp,select=TRUE,method="REML",fit=FALSE);

if (FALSE){
#--run cross validation-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
  set.seed(1111111);
  mdl = mdlNrml_ZE2D;
  dfrCrsVal = runCrossValidation(
                 mdl,
                 ks,
                 dfrData=dfrDatp,
                 numFolds=10,
                 concrv_opt=2,
                 doParallel=TRUE,
                 ncores=parallel::detectCores()-1,
                 max_ncmbs=NULL,
                 icmbs=NULL,
                 log=TRUE,
                 debug=FALSE);
  wtsUtilities::saveObj(dfrCrsVal,file.path(dirThs,"rda_GaussianModels_CrsVal.RData"));
}
if (FALSE){
#--evaluate best model-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
  source(file.path(dirThs,"../r_PlotStats_BestModels.R"));
  mdl = mdlNrml_ZE2D;
  if (!exists("dfrCrsVal")) dfrCrsVal = wtsUtilities::getObj(file.path(dirThs,"rda_GaussianModels_CrsVal.RData"));
  #--extract base model results
  dfrCrsVal1 = dfrCrsVal |> 
                 dplyr::filter(i==1) |> #--extract base model results
                 dplyr::select(fold,
                               base_rsqr=rsqr,base_aic=aic,base_rsqr_prd=rsqr_prd,
                               base_mspe_prd=mspe_prd,base_mase_prd=mase_prd);
  #--calculate stats differences by fold relative to the base model stats
  dfrCrsVald = dfrCrsVal |> 
                 dplyr::left_join(dfrCrsVal1,by=c("fold")) |> 
                 dplyr::mutate(impr_rsqr=rsqr-base_rsqr,
                               impr_aic=base_aic-aic,
                               impr_rsqr_prd=rsqr_prd-base_rsqr_prd,
                               impr_mspe_prd=(base_mspe_prd-mspe_prd)/base_mspe_prd,
                               impr_mase_prd=(base_mase_prd-mase_prd)/base_mase_prd);
  #--summarize stats differences (mean, median) over folds
  dfrCrsValdp = dfrCrsVald |> 
                 dplyr::group_by(i,frmla,smths) |> 
                 dplyr::summarize(scrSignifp=sum(as.numeric(signifp)),
                                  scrSignifs=sum(as.numeric(signifs)),
                                  scrConcrv_tst=sum(as.numeric(concrv_tst)),
                                  mn_impr_rsqr=mean(impr_rsqr),
                                  mn_impr_aic=mean(impr_aic),
                                  mn_impr_rsqr_prd=mean(impr_rsqr_prd),
                                  mn_impr_mspe_prd=mean(impr_mspe_prd),
                                  mn_impr_mase_prd=mean(impr_mase_prd),
                                  md_impr_rsqr=median(impr_rsqr),
                                  md_impr_aic=median(impr_aic),
                                  md_impr_rsqr_prd=median(impr_rsqr_prd),
                                  md_impr_mspe_prd=median(impr_mspe_prd),
                                  md_impr_mase_prd=median(impr_mase_prd)
                 ) |> 
                 dplyr::ungroup() |> 
                 dplyr::arrange(i); 
  #--find best model based on concurvity, significant smooths, and mn_impr_mspe_prd 
  #----by arranging all models
  dfrCrsValdp1 = dfrCrsValdp |> 
                   dplyr::arrange(desc(scrConcrv_tst),
                                  desc(scrSignifs),
                                  desc(mn_impr_mspe_prd));
  #--plot stats from models with better summary stats than the base model
  plotStats_BestModels(dfrCrsVald,dfrCrsValdp1);
  
  #--evaluate best model using full dataset
  best_idx = dfrCrsValdp1$i[1];#--index of best model in evaluated combinations
  best_mdl = evalBestModel(mdl,ks,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
  wtsUtilities::saveObj(dfrCrsValdp1,file.path(dirThs,"rda_GaussianModels_OrderedModels.RData"));
  wtsUtilities::saveObj(best_mdl,    file.path(dirThs,"rda_GaussianModels_BestModel.RData"));
}

#--------ALL 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
  #--        ti(d,t) + ti(d,f) + ti(d,s) + 
  #--        ti(t,f) + ti(t,s) + 
  #--        ti(f,s)
  # k = 10; k1 = 6; k2 = 6;
  # frmla  = lnR~ti(z,bs="ts",k=k)   + 
  #            ti(d,bs="ts",k=k1)   + ti(t,bs="ts",k=k1)   + ti(f,bs="ts",k=k1)   + ti(s,bs="ts",k=k1) +
  #            ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
  #            ti(d,t,bs="ts",k=k2) + ti(d,f,bs="ts",k=k2) + ti(d,s,bs="ts",k=k2) +
  #            ti(t,f,bs="ts",k=k2) + ti(t,s,bs="ts",k=k2)+
  #            ti(f,s,bs="ts",k=k2);
  # mdlNrml_All2D  = mgcv::gam(frmla,family=famTW,data=dfrDatp,select=TRUE,method="REML",fit=FALSE);
