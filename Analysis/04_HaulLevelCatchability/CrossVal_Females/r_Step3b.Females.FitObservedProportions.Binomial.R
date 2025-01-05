#--fit various models for ln(r) using mgcv to fit GAMs for FEMALES using the BINOMIAL distribution----
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

#--get censored data and prediction grids----
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Females.RData"));

#--remove zeros, infs, questionable observed Rs----
#dfrDatp   = lst$dfrDat |> dplyr::filter(obsR<10, is.finite(lnR),between(z,15,150));
dfrDatp   = lst$dfrDat |> dplyr::filter(between(z,15,130));

#--BINOMIAL regression  models for lnR----
famB = stats::binomial(link="logit");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = p~ti(z,bs="ts",k=k1)   +
             ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
             ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  # frmla  = p~ti(z,bs="ts",k=k1)   +
  #            ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2) ;
  mdlB_ZE2D  = mgcv::gam(frmla,family=famB,data=dfrDatp,select=FALSE,method="ML",fit=FALSE,
                         offset=lnq,weights=n);

#--run cross-validation using concurvity and other criteria to rank models
if (FALSE){
#--run cross validation-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_EvalAllModelsFunctions.R"));
  source(file.path(dirThs,"../r_RunCrossValidationFunctions.R"));
  source(file.path(dirThs,"../r_LogLikeFunctions.R"));
  setwd(dirThs);
  set.seed(1111111);
  mdl = mdlB_ZE2D;
  dfrCrsVal = runCrossValidation(
                 mdl,
                 ks,
                 dfrData=dfrDatp,
                 col_link="lnR",
                 numFolds=20,
                 selectBy="by_obs",
                 concrv_opt=2,
                 doParallel=TRUE,
                 ncores=parallel::detectCores()-1,
                 max_ncmbs=NULL,
                 icmbs=NULL,
                 log=TRUE,
                 debug=TRUE);
  wtsUtilities::saveObj(dfrCrsVal,file.path(dirThs,"rda_Step3b1.BinomialModels_CrsVal.RData"));
}
  
if (FALSE){
  if (!exists("dfrCrsVal")) dfrCrsVal = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b1.BinomialModels_CrsVal.RData"));
  #--keep scores, drop nested lists
  dfrScrs = dfrCrsVal |> dplyr::select(!c(lstMdl,lstTrn,lstTst)); rm(dfrCrsVal);
  wtsUtilities::saveObj(dfrScrs,file.path(dirThs,"rda_Step3b2.BinomialModels_Scrs.RData"));
}
  
if (FALSE){
#--evaluate best model-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_EvalAllModelsFunctions.R"));
  source(file.path(dirThs,"../r_RunCrossValidationFunctions.R"));
  source(file.path(dirThs,"../r_LogLikeFunctions.R"));
  source(file.path(dirThs,"../r_PlotStats_BestModels.R"));
  mdl = mdlB_ZE2D;
  if (!exists("dfrScrs")) dfrScrs = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b2.BinomialModels_Scrs.RData"));

  #--calculate mean scores
  dfrMnScrs = dfrScrs |> dplyr::group_by(i,smths) |> 
                dplyr::summarize(n=n(),
                                mnScrTst=mean(scrTst),
                                mnAIC=mean(aic),
                                mnRsrPrd=mean(rsqr_prd),
                                numConcrvTst=sum(concrv_tst)) |> 
                dplyr::arrange(dplyr::desc(numConcrvTst),dplyr::desc(mnScrTst)) |> 
                dplyr::filter(n>18,numConcrvTst>10) |> 
                dplyr::arrange(dplyr::desc(mnScrTst));
  sel_mdls = unique(c(dfrMnScrs[1:5,"smths"]$smths,"ti(z)"));
  p1 = ggplot(dfrScrs |> dplyr::filter(smths %in% sel_mdls) |> 
                dplyr::mutate(smths=factor(smths,levels=sel_mdls)),
              aes(x=smths,y=scrTst)) + geom_boxplot() + geom_point() + 
         geom_point(data=dfrMnScrs |> dplyr::filter(smths %in% sel_mdls) |> 
                      dplyr::mutate(smths=factor(smths,levels=sel_mdls)),
                    mapping=aes(x=smths,y=mnScrTst),shape=23,size=6) +
         geom_hline(data=dfrMnScrs |> dplyr::filter(smths %in% sel_mdls[1]),
                    aes(yintercept=mnScrTst),linetype=3) + 
         labs(y="GCV score",subtitle="binomial models") + 
         wtsPlots::getStdTheme() + 
         theme(axis.text.x=element_text(size=12,angle=345,hjust=0),
               axis.title.x=element_blank());
  ggsave("pltBestModels.Females.Binomial.pdf",plot=p1,width=6.5,height=4)
  
  #--evaluate best model
  best_smth = "ti(z)+ti(t)+ti(f)";#--user must determine this based on results above
  best_idx  = (dfrMnScrs |> dplyr::filter(smths==best_smth))$i;
  best_mdl = evalBestModel(mdl,ks,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
  
  wtsUtilities::saveObj(best_mdl,file.path(dirThs,"rda_Step3b3.BinomialModels_BestModel.RData"));
}

