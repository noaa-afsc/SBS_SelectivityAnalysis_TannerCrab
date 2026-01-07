#--use Thygesen et al approach to estimate relative selectivity ratio
##--using log-Gaussian Cox processes
require(TMB);
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04a_ThygesenIntercalibration");
dirCPP = file.path(dirThs,"Intercalibration/gearcalib/src");

#--compile dynamic library
setwd(dirCPP);
compile("gearcalib.cpp",CXX="clang++")
if(is.loaded("EvalDoubleFunObject")) dyn.unload("gearcalib.so") #<-may error out
dyn.load("gearcalib.so")
setwd(dirThs);

#--get data
##--dfrCPUE: fleet, y, gis_station, sampling_factor, area_swept_variable, x, m, s, z, n, val, type
dfrCPUE = wtsUtilities::getObj(file.path(dirPrj,"Analysis/01_SBS_Data/rda_Step3_SBS_CrabAbundance.RData"))$dfrCPUE;
dfrCPUE = dfrCPUE |> dplyr::mutate(gear=factor(fleet,levels=c("BSFRF","NMFS")),
                                   group=paste0(y,"_",gis_station),
                                   ASp=area_swept_variable/sampling_factor);
dfrGr = dfrCPUE |> dplyr::select(group,x,n,ASp) |> 
          dplyr::group_by(group,x) |> 
          dplyr::summarize(totN=sum(n,na.rm=FALSE),
                           totA=sum(ASp,na.rm=FALSE)) |> 
          dplyr::ungroup() |>
          dplyr::arrange(group,x) |> 
          dplyr::filter(totN>0,is.finite(totA)) |> dplyr::select(group,x);
dfrM.N = dfrGr |> dplyr::filter(x=="male") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,179)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,n) |> 
          tidyr::pivot_wider(names_from=z,values_from=n) |> dplyr::arrange(group,gear);
dfrM.A = dfrGr |> dplyr::filter(x=="male") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,179)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,ASp) |> 
          tidyr::pivot_wider(names_from=z,values_from=ASp) |> dplyr::arrange(group,gear);
dfrF.N = dfrGr |> dplyr::filter(x=="female") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,119)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,n) |> 
          tidyr::pivot_wider(names_from=z,values_from=n) |> dplyr::arrange(group,gear);
dfrF.A = dfrGr |> dplyr::filter(x=="female") |> 
          dplyr::inner_join(dfrCPUE |> dplyr::filter(dplyr::between(z,24,119)),by=c("group","x")) |> 
          dplyr::select(group,gear,z,ASp) |> 
          tidyr::pivot_wider(names_from=z,values_from=ASp) |> dplyr::arrange(group,gear);

#--source the R code
source(file.path(dirThs,"Intercalibration/gearcalib/R/gearcalib.R"));

#--convert to NEW gearcalib "d" format:
##--d$N: matrix with N in haul x size
##--d$SweptArea: matrix with area_swept/sampling_factor in haul x size
##--d$group: factor encoding hauls (size = nHauls)
##--d$Gear:  factor encoding gear  (size = nHauls)
###--males
zM = as.double(colnames(dfrM.N)[3:ncol(dfrM.N)]) + 3;
dM = list();
dM$L         = zM;
dM$N         = as.matrix(dfrM.N[,3:ncol(dfrM.N)]);
dM$SweptArea = as.matrix(dfrM.A[,3:ncol(dfrM.A)]);
dM$group     = factor(dfrM.N$group);
dM$Gear      = dfrM.N$gear; 
###--females
zF = as.double(colnames(dfrF.N)[3:ncol(dfrF.N)]) + 3;
dF = list();
dF$L         = zF;
dF$N         = as.matrix(dfrF.N[,3:ncol(dfrF.N)]);
dF$SweptArea = as.matrix(dfrF.A[,3:ncol(dfrF.A)]);
dF$group     = factor(dfrF.N$group);
dF$Gear      = dfrF.N$gear; 

#--fit the model for males
fitM = gearcalibFit(dM,fit0=TRUE);
botM = boot(dM);
pltM = plot.gearcalibFit(fitM,boot=botM,Lvec=zM,ymax=NULL)

#--fit the model for females
fitF = gearcalibFit(dF,fit0=TRUE);
botF = boot(dF);
pltF = plot.gearcalibFit(fitF,boot=botF,Lvec=zF,ymax=1)

res = list(m=list(fit=fitM,boot=botM,plots=pltM),
           f=list(fit=fitF,boot=botF,plots=pltF));
wtsUtilities::saveObj(res,file.path(dirThs,"results.RData"));

