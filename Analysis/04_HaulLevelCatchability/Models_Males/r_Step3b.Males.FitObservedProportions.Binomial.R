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

#--BINOMIAL regression  models for lnR----
famB = stats::binomial(link="logit");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  ks=c(20,10);
  k1 = ks[1]; k2 = ks[2];
  frmla  = p~ti(z,bs="ts",k=k1)   +
             ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
             ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  mdlB_ZE2D  = mgcv::gam(frmla,family=famB,data=dfrDatp,select=TRUE,method="REML",fit=FALSE,
                         offset=lnq,weights=n);

#--run cross-validation using concurvity and other criteria to rank models
if (FALSE){
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
#--run cross validation-------------------------------------------------
  set.seed(1111111);
  mdl = mdlB_ZE2D;
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
  wtsUtilities::saveObj(dfrCrsVal,file.path(dirThs,"rda_Step3b1.BinomialModels_CrsVal.RData"));
}
if (FALSE){
#--evaluate best model-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
  source(file.path(dirThs,"../r_PlotStats_BestModels.R"));
  mdl = mdlB_ZE2D;
  if (!exists("dfrCrsVal")) dfrCrsVal = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b1.BinomialModels_CrsVal.RData"));
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
  # #--find best model based on concurvity, significant smooths, and mn_impr_mspe_prd 
  # #----by arranging all models
  # dfrCrsValdp1 = dfrCrsValdp |> 
  #                  dplyr::arrange(desc(scrConcrv_tst),
  #                                 desc(scrSignifs),
  #                                 desc(mn_impr_mspe_prd));
  #--find best model based on concurvity, significant smooths, and mn_impr_mspe_prd 
  #----by arranging all models
  dfrCrsValdp1 = dfrCrsValdp |> 
                   dplyr::arrange(desc(md_impr_mspe_prd));
  dfrCrsValdp1 = dfrCrsValdp1[1:min(10,nrow(dfrCrsValdp1)),];
  #--plot stats from models with better summary stats than the base model
  plotStats_BestModels(dfrCrsVald,dfrCrsValdp1);
  
  best_idx = dfrCrsValdp1$i[2];#--index of best model in evaluated combinations
  best_mdl = evalBestModel(mdl,ks,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
  wtsUtilities::saveObj(dfrCrsValdp1,file.path(dirThs,"rda_Step3b1.BinomialModels_OrderedModels.RData"));
  wtsUtilities::saveObj(best_mdl,    file.path(dirThs,"rda_Step3b1.BinomialModels_BestModel.RData"));
}

