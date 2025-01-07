#--calculate predicted annual NMFS selectivity functions from haul-specific versions----
#----for best ln-scale Tweedie model
require(dplyr);
require(ggplot2);
require(mgcv);

#--read project setup info----
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

if (FALSE){  
  #--load sediment data (phi, sorting) interpolated to haul locations----
  s  = wtsUtilities::getObj(file.path(dirPrj,"rda_ProjectSetup.RData"));
  dfrHD_sed = wtsUtilities::getObj(file.path(s$dirs$SedAnls,
                                             "rda_dfrHD_All_NMFS_WithInterpolatedSedValues.RData")) |> 
                sf::st_drop_geometry() |>
                dplyr::filter(YEAR>=1982);
  
  #--load interpolation grids----
  lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Males.RData"));
  
  #--get best lnTweedie model
  best_mdl = wtsUtilities::getObj(file.path(dirThs,"rda_Step3b3b.LnTweedieModels_BestModel.RData"));
  
  #--get prediction "grid"----
  dfrPrdGrd = dfrHD_sed |> 
                dplyr::select(y=YEAR,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting,HAULJOIN) |> 
                dplyr::cross_join(tibble::tibble(z=lst$grids$z) |> dplyr::filter(z>=25)) |>
                dplyr::mutate(h="any");
  
  #--make predictions----
  source(file.path(dirThs,"../r_PredictionsAndPlots.R"));
  dfrPrd = prdMod(best_mdl,trms=c("all"),type="link",lst=dfrPrdGrd,p=0.10) |> 
              dplyr::mutate(emp_sel=exp(emp_sel),
                            lci=exp(lci),
                            uci=exp(uci));
  wtsUtilities::saveObj(dfrPrd,file.path(dirThs,"rda_Step4a.Males.LnTweedieModel_Best.AllPredicted.SelFcns.RData"));
}

#--calculate unweighted averages across hauls by year
dfrPrd = wtsUtilities::getObj(file.path(dirThs,"rda_Step4a.Males.LnTweedieModel_Best.AllPredicted.SelFcns.RData"));
dfrPrdAvg = dfrPrd |> 
              dplyr::group_by(y,z) |> 
              dplyr::summarize(emp_sel=mean(emp_sel,na.rm=TRUE)) |> 
              dplyr::ungroup();

yrs = c(1982,2023)
ggplot(dfrPrd |> dplyr::filter(y %in% yrs),aes(x=z,y=emp_sel,colour=as.character(y),group=HAULJOIN)) + 
  geom_line(linetype=3) + 
  geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
  geom_line(data=dfrPrdAvg |> dplyr::filter(y %in% yrs),aes(group=NA),linewidth=2) + 
  labs(x="size (mm CW)",y="predicted selectivity",colour="year") + 
  wtsPlots::getStdTheme()+
  theme(legend.position="none");
yrs = 1982:2023
ggplot(dfrPrdAvg |> dplyr::filter(y %in% yrs),aes(x=z,y=emp_sel,colour=as.character(y))) + 
  geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
  geom_line(data=,linewidth=1) + 
  labs(x="size (mm CW)",y="predicted selectivity",colour="year") + 
  wtsPlots::getStdTheme()+
  theme(legend.position="none");

#--calculate CPUE-weighted averages across hauls by year
dfrCPUE = wtsUtilities::getObj(file.path(dirThs,"../rda_Step4.dfrCPUE_NMFS.RData")) |> 
            dplyr::filter(SEX=="MALE") |> 
            dplyr::mutate(N=numCPUE*AREA_SWEPT_VARIABLE) |> 
            dplyr::select(y=YEAR,HAULJOIN,z=SIZE,N,numCPUE);
dfrPrd1 = dfrPrd |> 
            dplyr::inner_join(dfrCPUE,by=c("y","HAULJOIN","z"));
dfrPrd1Avg = dfrPrd1 |> 
               dplyr::group_by(y,z) |> 
               dplyr::summarise(emp_selC=sum(numCPUE*emp_sel,na.rm=TRUE)/sum(numCPUE),
                                emp_selN=sum(N*emp_sel,na.rm=TRUE)/sum(N)) |> 
               dplyr::ungroup();
yrs = 1982:2023
ggplot(dfrPrd1Avg |> dplyr::filter(y %in% yrs),aes(x=z,y=emp_selC,colour=as.character(y))) + 
  geom_hline(yintercept=c(0.0,0.5,1.0),linetype=3) + 
  geom_line(data=,linewidth=1) + 
  labs(x="size (mm CW)",y="predicted selectivity",colour="year") + 
  wtsPlots::getStdTheme()+
  theme(legend.position="none");
