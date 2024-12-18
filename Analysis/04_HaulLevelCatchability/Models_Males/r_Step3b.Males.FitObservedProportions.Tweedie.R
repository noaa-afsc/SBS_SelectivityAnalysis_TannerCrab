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
dfrDatp   = lst$dfrDat |> dplyr::filter(obsR<10, is.finite(lnR),between(z,15,150));

#--TWEEDIE regression models for lnR---------------------
famTW = mgcv::tw(link="log");
#--------ALL Z 2-WAY INTERACTIONS--------------------------
  #--ln(r) = ti(z) + 
  #          ti(d) + ti(t) + ti(f) + ti(s) +
  #--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s)
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = obsR~ti(z,bs="ts",k=k1)   +
             ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2)   + ti(f,bs="ts",k=k2)   + ti(s,bs="ts",k=k2) +
             ti(z,d,bs="ts",k=c(k1,k2)) + ti(z,t,bs="ts",k=c(k1,k2)) + ti(z,f,bs="ts",k=c(k1,k2)) + ti(z,s,bs="ts",k=c(k1,k2));
  mdlTW_ZE2D = mgcv::gam(frmla,family=famTW,data=dfrDatp,select=TRUE,method="REML",fit=FALSE);

if (FALSE){
#--run cross validation-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
  set.seed(1111111);
  mdl = mdlTW_ZE2D;
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
  wtsUtilities::saveObj(dfrCrsVal,file.path(dirThs,"rda_Step3b.TweedieModels_CrsVal.RData"));
}
if (FALSE){
#--evaluate best model-------------------------------------------------
  source(file.path(dirThs,"../r_gam.prefit.functions.R"));
  source(file.path(dirThs,"../r_SelectModelByConcurvityFunctions.R"));
  source(file.path(dirThs,"../r_PlotStats_BestModels.R"));
  mdl = mdlTW_ZE2D;
  if (!exists("dfrCrsVal")) dfrCrsVal = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b.TweedieModels_CrsVal.RData"));
  #--extract base model results
  dfrCrsVal1 = dfrCrsVal |> 
                 dplyr::filter(i==1) |> #--extract base model results
                 dplyr::select(fold,
                               base_rsqr=rsqr,base_aic=aic,base_rsqr_prd=rsqr_prd,
                               base_mspe_prd=mspe_prd,base_mase_prd=mase_prd);
  #--calculate stats differences by fold relative to the base model stats
  dfrCrsVald = dfrCrsVal |> 
                 dplyr::left_join(dfrCrsVal1,by=c("fold")) |> 
                 dplyr::mutate(impr_rsqr=rsqr-base_rsqr,            #--improve > 0
                               impr_aic=base_aic-aic,               #--improve > 0
                               impr_rsqr_prd=rsqr_prd-base_rsqr_prd,#--improve > 0
                               impr_mspe_prd=(base_mspe_prd-mspe_prd)/base_mspe_prd,  #--improve > 0
                               impr_mase_prd=(base_mase_prd-mase_prd)/base_mase_prd); #--improve > 0
  #--summarize stats differences (mean, median) over folds
  dfrCrsValdp = dfrCrsVald |> 
                 dplyr::group_by(i,frmla,smths) |> 
                 dplyr::summarize(scrSignifp=sum(as.numeric(signifp)),
                                  scrSignifs=sum(as.numeric(signifs)),
                                  scrConcrv_tst=sum(as.numeric(concrv_tst)),
                                  mn_impr_rsqr=mean(impr_rsqr),                #--improve > 0
                                  mn_impr_aic=mean(impr_aic),                  #--improve > 0
                                  mn_impr_rsqr_prd=mean(impr_rsqr_prd),        #--improve > 0
                                  mn_impr_mspe_prd=mean(impr_mspe_prd),        #--improve > 0
                                  mn_impr_mase_prd=mean(impr_mase_prd),        #--improve > 0
                                  md_impr_rsqr=median(impr_rsqr),              #--improve > 0
                                  md_impr_aic=median(impr_aic),                #--improve > 0
                                  md_impr_rsqr_prd=median(impr_rsqr_prd),      #--improve > 0
                                  md_impr_mspe_prd=median(impr_mspe_prd),      #--improve > 0
                                  md_impr_mase_prd=median(impr_mase_prd)       #--improve > 0
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
  
  best_idx = dfrCrsValdp1$i[1];#--index of best model in evaluated combinations
  best_mdl = evalBestModel(mdl,ks,best_idx);
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = best_mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);
  plts = getModelPlots(best_mdl);
  wtsUtilities::saveObj(dfrCrsValdp1,file.path(dirThs,"rda_Step3b.TweedieModels_OrderedModels.RData"));
  wtsUtilities::saveObj(best_mdl,    file.path(dirThs,"rda_Step3b.TweedieModels_BestModel.RData"));
}

if (FALSE){
#--best model + haul-level random effects
  lvls = c("any",unique(dfrDatp$h));
  dfrDatpp = dfrDatp |> dplyr::mutate(h=factor(h,levels=lvls));
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = obsR~ti(z,bs="ts",k=k1) + ti(t,bs="ts",k=k2) + ti(f,bs="ts",k=k2) + 
                  s(s,h,bs="fs",k=k2);
  mdl_bestRE = mgcv::gam(frmla,family=famTW,data=dfrDatpp,select=TRUE,method="REML",
                         drop.unused.levels=FALSE);#--use for RE with "any" factor level
  wtsUtilities::saveObj(mdl_bestRE, file.path(dirThs,"rda_Step3b.TweedieModels_BestModelRE.RData"));
#----function to predict values based on a model
  prdMod<-function(mdl,trms,lst,type="response",keep=NULL,p=0.05){
    dfr = wtsMGCV::createGridTbl(lst);
    if (any(trms=="all")){
      #--add intercept and all smooth terms
      trmsp = "(Intercept)";
      trms = wtsMGCV::getSmoothTerms(mdl);
      for (trm in trms) trmsp = c(trmsp,trm);
      trms = trmsp;
    }
    prd = dplyr::bind_cols(
              dfr,
              tibble::as_tibble(
                mgcv::predict.gam(mdl,dfr,type=type,terms=trms,se.fit=TRUE),
              ) |> 
              dplyr::mutate(type="fit",
                            lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                            uci=qnorm(p,fit,se.fit,lower.tail=FALSE),
                            terms=paste(trms,collapse=" + "))
          ) |> 
            dplyr::rename(emp_sel=fit);
    if (!is.null(keep)){
      prd = prd |> 
              dplyr::distinct(pick(tidyselect::any_of(keep),emp_sel,se.fit,lci,uci,terms));
      drop = names(lst)[!(names(lst) %in% keep)];
      for (drp in drop) prd[[drp]] = NA;
    }
    return(prd);
  }
  plotMod<-function(tmp,ylims=c(0,1.5)){
    if (all(is.na(tmp$y))) tmp$y = "all";
    p = ggplot(tmp,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y));
    if ("n" %in% names(tmp)){
      p = p + geom_point(aes(size=n)) + scale_size_area() + 
              geom_line();
    }
    p = p + 
           geom_ribbon(alpha=0.3) + 
           geom_line() + 
           geom_hline(yintercept=0.5,linetype=3) + 
           scale_y_continuous(limits=ylims,oob=scales::squish) + 
           labs(x="size (mm CW)",y="empirical\nselectivity",
                colour="study\nyear",fill="study\nyear",size="crab\nsampled") + 
           theme(legend.position="inside",
                 legend.position.inside=c(0.01,0.99),
                 legend.justification.inside=c(0,1),
                 legend.byrow=TRUE,
                 legend.box="horizontal");
    return(p);
  }
  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"))
  dfrPrd = prdMod(mdl_bestRE,trms=c("all"),type="response",lst=grdPrd);
  plotMod(dfrPrd)
}
  
#--------ALL 2-WAY INTERACTIONS--------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s)
# k = 10; k1 = 6; k2 = 6;
# frmla  = obsR~ti(z,bs="ts",k=k)   + 
#            ti(d,bs="ts",k=k1)   + ti(t,bs="ts",k=k1)   + ti(f,bs="ts",k=k1)   + ti(s,bs="ts",k=k1) +
#            ti(z,d,bs="ts",k=c(k,k1)) + ti(z,t,bs="ts",k=c(k,k1)) + ti(z,f,bs="ts",k=c(k,k1)) + ti(z,s,bs="ts",k=c(k,k1)) +
#            ti(d,t,bs="ts",k=k2) + ti(d,f,bs="ts",k=k2) + ti(d,s,bs="ts",k=k2) +
#            ti(t,f,bs="ts",k=k2) + ti(t,s,bs="ts",k=k2)+
#            ti(f,s,bs="ts",k=k2);
