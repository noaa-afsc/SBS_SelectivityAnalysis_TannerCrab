#--compare male selectivity curves
require(DHARMa);
require(dplyr);
require(ggplot2);
require(gratia)
require(mgcv);
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
source(file.path(dirThs,"..","r_PredictionsAndPlots.R"));

#--get prediction grids----
lst    = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Males.RData"));
grdPrd = list(z=lst$grids$z,d=lst$meds$d,t=lst$meds$t,f=lst$meds$f,s=lst$meds$s,h=factor("any"));

#--get models
paths = list(Binomial="rda_Step3b3b.BinomialModels_BestModel.RData",
             Tweedie ="rda_Step3b3b.LnTweedieModels_BestModel.RData",
             BinomialRE="rda_Step3b1.BinomialModels_RE.RData",
             TweedieRE ="rda_Step3b1.LnTweedieModels_RE.RData");
lstMdls = list();
for (pth in names(paths)) lstMdls[[pth]]=wtsUtilities::getObj(file.path(dirThs,paths[[pth]]));

#--get base selectivity estimates excluding covariate effects
lstPrd = list();
for (mdl in names(lstMdls)){
  dfrPrd = prdMod(lstMdls[[mdl]],trms=c("(Intercept)","ti(z)","s(z)"),type="link",lst=grdPrd,p=0.10) |> 
             dplyr::mutate(model_name=mdl);
  lstPrd[[mdl]] = dfrPrd;
  rm(dfrPrd);
}
dfrSelFcns = dplyr::bind_rows(lstPrd) |> 
                dplyr::mutate(emp_sel=exp(emp_sel),
                              lci=exp(lci),
                              uci=exp(uci));
rm(lstPrd);

#--compare estimated selectivity curves without covariates or REs
ggplot(dfrSelFcns,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=model_name,fill=model_name)) + 
  geom_ribbon(alpha=0.3) + geom_line() + 
  geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
  scale_y_continuous(limits=c(0,1.2),oob=scales::squish) + 
  labs(x="size (mm CW)",y="estimated selectivity",
      colour="model",fill="model") + 
  wtsPlots::getStdTheme() + 
  theme(legend.position="inside",
       legend.position.inside=c(0.01,0.99),
       legend.justification.inside=c(0,1),
       legend.byrow=TRUE,
       legend.box="horizontal");

#--get base selectivity estimates at median covariate values and "all" REs
lstPrd = list();
for (mdl in names(lstMdls)){
  dfrPrd = prdMod(lstMdls[[mdl]],trms=c("all"),type="link",lst=grdPrd,p=0.10) |> 
             dplyr::mutate(model_name=mdl);
  lstPrd[[mdl]] = dfrPrd;
  rm(dfrPrd);
}
dfrSelFcns = dplyr::bind_rows(lstPrd) |> 
                dplyr::mutate(emp_sel=exp(emp_sel),
                              lci=exp(lci),
                              uci=exp(uci));
rm(lstPrd);

#--compare estimated selectivity curves
ggplot(dfrSelFcns,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=model_name,fill=model_name)) + 
  geom_ribbon(alpha=0.3) + geom_line() + 
  geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
  scale_y_continuous(limits=c(0,1.2),oob=scales::squish) + 
  labs(x="size (mm CW)",y="estimated selectivity",
      colour="model",fill="model") + 
  wtsPlots::getStdTheme() + 
  theme(legend.position="inside",
       legend.position.inside=c(0.01,0.99),
       legend.justification.inside=c(0,1),
       legend.byrow=TRUE,
       legend.box="horizontal");

