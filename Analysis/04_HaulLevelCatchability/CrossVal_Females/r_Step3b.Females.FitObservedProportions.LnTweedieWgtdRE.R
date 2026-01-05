#--fit various models for ln(r) using mgcv to fit GAMs for females using the TWEEDIE distribution----
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);

#--get censored data and prediction grids----
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Females.RData"));

#--remove questionable observed Rs----
dfrDatp   = lst$dfrDat |> dplyr::filter(nB>2, between(z,15,130));
#--change hauls to a factor with additional level "any"
lvls = c("any",unique(dfrDatp$h));
dfrDatpp = dfrDatp |> dplyr::mutate(h=factor(h,levels=lvls));

#--TWEEDIE (using mgcv::tw) regression  models for lnR----
famTW = mgcv::tw(link="log");
#--------RE Model--------------------------
  ks=c(10,8);
  k1 = ks[1]; k2 = ks[2];
  frmla  = obsR~s(z,bs="ts",k=k1) + ti(z,h,bs="fs",k=k1);
  frmla  = obsR~s(z,bs="ts",k=k1) + s(h,bs="re",k=k1);
if (FALSE){
  mdl  = mgcv::gam(frmla,family=famTW,data=dfrDatpp,select=FALSE,method="ML",fit=TRUE,
                        drop.unused.levels=FALSE,weights=nB);
  wtsUtilities::saveObj(mdl,file.path(dirThs,"rda_Step3b1.LnTweedieModelsWgtd_RE.RData"));
}

#--make plots, etc.
if (FALSE){
  source(file.path(dirThs,"..","r_PredictionsAndPlots.R"));
  mdl = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b1.LnTweedieModelsWgtd_RE.RData"));
  #--diagnostic plots
  simResids <- DHARMa::simulateResiduals(fittedModel = mdl, plot = F, n=1000);
  DHARMa::plotQQunif(simResids);
  DHARMa::plotResiduals(simResids);

  grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=lvls)
  dfrPrd = prdMod(mdl,trms=c("all"),type="link",lst=grdPrd,p=0.10) |> 
            dplyr::mutate(emp_sel=exp(emp_sel),
                          lci=exp(lci),
                          uci=exp(uci));
  plotMod(dfrPrd) 
  ggplot(dfrPrd,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,group=h)) + 
    geom_line(linetype=3) + 
    geom_line(data=dfrPrd |> dplyr::filter(h=="any"),colour="blue") + 
    geom_ribbon(data=dfrPrd |> dplyr::filter(h=="any"),alpha=0.2,fill="blue") + 
    geom_point(aes(x=z,y=log(obsR)),data=dfrDatpp |> dplyr::filter(is.finite(lnR)),inherit.aes=FALSE) + 
    scale_x_continuous(limits=c(0,130)) + 
    scale_y_continuous(limits=c(-6,3))

    ggplot(dfrPrd |> dplyr::filter(h=="any"),aes(x=z,y=emp_sel,ymin=lci,ymax=uci,group=h)) + 
    geom_line() + geom_ribbon(alpha=0.2) + 
    geom_point(aes(x=z,y=log(obsR)),data=dfrDatpp,inherit.aes=FALSE) + 
    scale_x_continuous(limits=c(25,130)) + 
  
  
  wtsUtilities::saveObj(dfrPrd,file.path(dirThs,"rda_Step3b2p.LnTweedieModelsWgtd_RE.RData"));
  ggplot(dfrDatp,aes(x=z,y=obsR,size=nB)) + geom_point(alpha=0.2)
}

