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
#dfrDatp   = lst$dfrDat |> dplyr::filter(between(z,15,150));

#--TWEEDIE (using mgcv::tw) regression  models for lnR----
famTW = mgcv::tw(link="log");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
 #--         ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  ks=c(20,10);
  k1 = ks[1]; k2 = ks[2];
  frmla  = obsR~ti(z,bs="ts",k=k1)   +
                 ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
                 ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  # frmla  = obsR~ti(z,bs="ts",k=k1)   +
  #            ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2) ;
  mdl_ZE2D  = mgcv::gam(frmla,family=famTW,data=dfrDatp,select=FALSE,method="ML",fit=FALSE);

#--run cross-validation using concurvity and other criteria to rank models
if (FALSE){
#--run cross validation-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_EvalAllModelsFunctions.R"));
  source(file.path(dirThs,"../r_RunCrossValidationFunctions.R"));
  source(file.path(dirThs,"../r_LogLikeFunctions.R"));
  setwd(dirThs);
  set.seed(1111111);
  mdl = mdl_ZE2D;
  runCrossValidation(
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
    debug=TRUE
  );
}
  
if (FALSE){
  #--process cross-validation folds to extract scores
  dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
  numFolds = 20;
  lstCrsVal = list();
  for (f in 1:numFolds){
    dfrFld = wtsUtilities::getObj(file.path(dirThs,"ModelResults_LnTweedieModels",
                                            paste0("fold_",f,".RData"))) |> 
               dplyr::select(!c(lstMdl,lstTrn,lstTst));
    lstCrsVal[[f]] = dfrFld;
    rm(dfrFld)
  }
  dfrScrs = dplyr::bind_rows(lstCrsVal);
  rm(lstCrsVal);
  wtsUtilities::saveObj(dfrScrs,file.path(dirThs,"rda_Step3b2.LnTweedieModels_Scrs.RData"));
}

  
if (FALSE){
#--evaluate mean scores and choose best model-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_EvalAllModelsFunctions.R"));
  source(file.path(dirThs,"../r_EvalBestModel.R"));
  mdl = mdl_ZE2D;
  if (!exists("dfrScrs")) dfrScrs = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b2.LnTweedieModels_Scrs.RData"));
  
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
  wtsUtilities::saveObj(dfrMnScrs,file.path(dirThs,"rda_Step3b3a.LnTweedieModels_MnScrs.RData"));
  
  #--compare top 5 models + base by scores
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
  ggsave("pltBestModels.Males.LnTweedie.pdf",plot=p1,width=6.5,height=4);
  
  #--evaluate best model
  best_smth = "ti(z)+ti(f)";#--user must determine this based on results above
  best_idx  = (dfrMnScrs |> dplyr::filter(smths==best_smth))$i;
  best_mdl = evalBestModel(mdl,ks,best_idx);
  wtsUtilities::saveObj(best_mdl,file.path(dirThs,"rda_Step3b3b.LnTweedieModels_BestModel.RData"));
}

if (FALSE){
  source(file.path(dirThs,"../r_PlotModelSmooths.R"));
  source(file.path(dirThs,"../r_PredictionsAndPlots.R"));
  best_mdl = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b3b.LnTweedieModels_BestModel.RData"));
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts1 = plotModelSmooths(best_mdl);
  
  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  dfrPrd = prdMod(best_mdl,trms=c("(Intercept)","ti(z)"),type="link",lst=grdPrd,p=0.10) |> 
            dplyr::mutate(emp_sel=exp(emp_sel),
                          lci=exp(lci),
                          uci=exp(uci));
  plts2 = plotMod(dfrPrd);
  dfrPrd = prdMod(best_mdl,trms=c("all"),type="link",lst=grdPrd,p=0.10) |> 
            dplyr::mutate(emp_sel=exp(emp_sel),
                          lci=exp(lci),
                          uci=exp(uci));
  plts3 = plotMod(dfrPrd);
}
  
if (FALSE){
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_EvalAllModelsFunctions.R"));
  source(file.path(dirThs,"../r_EvalBestModel.R"));
  source(file.path(dirThs,"../r_PredictionsAndPlots.R"));
  mdl = mdl_ZE2D;
  dfrMnScrs = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b3a.LnTweedieModels_MnScrs.RData"));
  sel_mdls = unique(c(dfrMnScrs[1:5,"smths"]$smths,"ti(z)"));
  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  lstMdls = list();
  for (sel_mdl in sel_mdls[!is.na(sel_mdls)]){
    #--sel_mdl = sel_mdls[!is.na(sel_mdls)][1];
    idx  = (dfrMnScrs |> dplyr::filter(smths==sel_mdl))$i;
    lstMdls[[sel_mdl]] = evalBestModel(mdl,ks,idx);
  }

  lstPrd = list();
  for (sel_mdl in sel_mdls[!is.na(sel_mdls)]){
    lstPrd[[sel_mdl]]  = prdMod(lstMdls[[sel_mdl]],trms=c("(Intercept)","ti(z)"),type="link",lst=grdPrd,p=0.10) |> 
                           dplyr::mutate(trms=sel_mdl);
  }
  dfrPrd = dplyr::bind_rows(lstPrd) |> 
            dplyr::mutate(emp_sel=exp(emp_sel),
                          lci=exp(lci),
                          uci=exp(uci));
  ggplot(dfrPrd,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=trms,fill=trms)) + 
    geom_ribbon(alpha=0.3) + geom_line() + 
    geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
    scale_y_continuous(limits=c(0,1.5),oob=scales::squish) + 
    labs(x="size (mm CW)",y="base selectivity",colour="model\nterms",fill="model\nterms") + 
    wtsPlots::getStdTheme() + 
    theme(legend.position="inside",
          legend.position.inside=c(0.01,0.99),
          legend.justification.inside=c(0,1));
}

