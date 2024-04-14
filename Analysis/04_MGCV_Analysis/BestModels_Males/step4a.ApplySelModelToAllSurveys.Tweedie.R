#--apply selectivity model to NMFS survey data for single sex
#----TWEEDIE distribution
require(dplyr);
require(ggplot2);
require(Hmisc);
require(mgcv);
require(tcsamSurveyData);

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
source(file.path(dirThs,"../r_PlotEstimatedSelectivityFunctions.R"))

#--define sex
x = "male";
#--get best model
mdl = wtsUtilities::getObj(file.path(dirThs,"rda_TweedieModels_BestModel.RData"));


#--define grid for sizes
grd_z = seq(17.5,152.5,5);
dfrZ = tibble::tibble(z=grd_z);

#--get survey haul data with spatial covariates depth, temperature, phi, and sorting
#--(year,hauljoin,depth,temp,phi,sorting)
#--expand table with size grid info
dfrHDwSCs = wtsUtilities::getObj(file.path(dirThs,
                                           "../../03_Sediment_Analyses/rda_dfrHD_All_NMFS_WithInterpolatedSedValues.RData")) |> 
              sf::st_drop_geometry() |>
              dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting) |>
              dplyr::filter(y>=1982) |> 
              dplyr::mutate(yd=factor(10*floor(y/10))) |> 
              tidyr::expand_grid(dfrZ);

#--predict haul-specific selectivity ratio (year,hauljoin,size,R)
prd = predict(mdl,newdata=dfrHDwSCs,se.fit=TRUE,type="link");#--return predictions on lnR scale 
seR = sqrt(exp(prd$se.fit^2)-1)*exp(prd$fit+(prd$se.fit^2)/2);

#--predicted selectivity by haul--
seFactor = 2;
dfrRbyH = dfrHDwSCs |> 
            dplyr::mutate(prdLnR=prd$fit,  #--predicted ln-scale selectivity ratio
                          seLnR=prd$se.fit,#--standard error
                          prdR=exp(prdLnR),#--arithmetic-scale predicted selectivity ratio
                          seR=seR,         #--arithmetic-scale standard error
                          lower=exp(prdLnR-seFactor*seLnR),
                          upper=exp(prdLnR+seFactor*seLnR));

#--get weighting by NMFS haul catch-at-size
source(file.path(dirThs,"../r_GetNumbersCaught_NMFS.R"))
dfrZCsByH = getNumbersCaught_NMFS(c(1982,2023),
                                  sex=toupper(x),
                                  sz_rng=c(15,154),
                                  cutpts = seq(-0.5,204.5,5));

#--add in haul-level catch-at-size numbers
dfrRbyHp = dfrRbyH |> 
               dplyr::inner_join(dfrZCsByH |> dplyr::select(HAULJOIN,bin,numIndivs), 
                                 by=c("h"="HAULJOIN","z"="bin")) |> 
               dplyr::mutate(invwgt1=seR,
                             invwgt2=1/numIndivs,
                             invwgt3=seR/numIndivs);
#--number and std. dev-weighted mean selectivity across all hauls--
dfrMnRp = dfrRbyHp |> 
              dplyr::group_by(z) |> 
              dplyr::summarize(mnR1=Hmisc::wtd.mean(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE),
                               lower1=Hmisc::wtd.quantile(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper1=Hmisc::wtd.quantile(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE,probs=0.90),
                               mnR2=Hmisc::wtd.mean(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE),
                               lower2=Hmisc::wtd.quantile(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper2=Hmisc::wtd.quantile(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE,probs=0.90),
                               mnR3=Hmisc::wtd.mean(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE),
                               lower3=Hmisc::wtd.quantile(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper3=Hmisc::wtd.quantile(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE,probs=0.90)
                              ) |> 
              dplyr::ungroup();
p = plotWgtdSels(dfrRbyHp,dfrMnRp,1,"1/se(R)"); print(p);
p = plotWgtdSels(dfrRbyHp,dfrMnRp,2,"N");       print(p);
p = plotWgtdSels(dfrRbyHp,dfrMnRp,3,"N/se(R)"); print(p);

#--number and std. dev-weighted mean selectivity across all hauls, by year
dfrMnRbyYp = dfrRbyHp |> 
              dplyr::group_by(yd,y,z) |> 
              dplyr::summarize(mnR1=Hmisc::wtd.mean(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE),
                               lower1=Hmisc::wtd.quantile(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper1=Hmisc::wtd.quantile(prdR,weights=1/invwgt1,normwt=FALSE,na.rm=TRUE,probs=0.90),
                               mnR2=Hmisc::wtd.mean(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE),
                               lower2=Hmisc::wtd.quantile(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper2=Hmisc::wtd.quantile(prdR,weights=2/invwgt2,normwt=FALSE,na.rm=TRUE,probs=0.90),
                               mnR3=Hmisc::wtd.mean(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE),
                               lower3=Hmisc::wtd.quantile(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE,probs=0.10),
                               upper3=Hmisc::wtd.quantile(prdR,weights=3/invwgt3,normwt=FALSE,na.rm=TRUE,probs=0.90)
                              ) |> 
              dplyr::ungroup() |> 
              dplyr::mutate(y=as.factor(y));
p = plotWgtdSelsByY(dfrRbyHp,dfrMnRbyYp,1,"1/se(R)"); print(p);
p = plotWgtdSelsByY(dfrRbyHp,dfrMnRbyYp,1,"N");       print(p);
p = plotWgtdSelsByY(dfrRbyHp,dfrMnRbyYp,1,"N/se(R)"); print(p);



