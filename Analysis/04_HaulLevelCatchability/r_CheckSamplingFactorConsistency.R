##--compare results for dfrPropsAllYears from r_Step1b_DoPreliminaryCalcs.R
##--for  changes to how SAMPLING_FACTOR is treated as of 2026-02-02:
###--Previously, SAMPLING_FACTOR was assumed to be the same for all crab in a size bin, sex, haul combination
####--and treated as a grouping variable when calculating the total number of crab in a size bin, sex, haul combination
###--However, this was found to be false in a relatively few (but still probably significant) number of cases
####--To remedy this, SAMPLING_FACTOR was summed at the same time as numIndivs to give 
####--the total number of sampled crab (i.e., sum(numIndivs)) 
####--and the total SAMPLING_FACTOR (i.e. sum(SAMPLING_FACTOR)) for each size bin, sex, haul combination
###--NOTE: the effective SAMPLING_FACTOR (as in a multiplier on the total number of crab sampled in a size bin, sex, haul combination 
####--to get the total number of crab actually caught in the size bin, sex, haul combination) is 
####--effSF = sum(SAMPLING_FACTOR)/sum(numIndivs) --i.e., just the average sampling factor
require(ggplot2);
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04_HaulLevelCatchability")
dfrPropsAllYears = wtsUtilities::getObj(file.path(dirThs,"rda_Step1b_dfrPropsAllSBSYears.RData"));          #--new
obj              = wtsUtilities::getObj(file.path(dirThs,"rda_Step1b_dfrPropsAllSBSYears.20260111.RData")); #--old

dfrCmp = dplyr::bind_rows(dfrPropsAllYears |> dplyr::mutate(case="new"),
                          obj |> dplyr::mutate(case="old")) |> 
           dplyr::select(case,y=YEAR,h=HAULJOIN,x=SEX,z=SIZE,nN=numNMFS,sfN=sfNMFS,nB=numBSFRF,sfB=sfBSFRF);
dfrCmpW = dfrCmp |> tidyr::pivot_wider(names_from="case",values_from=c("nN","sfN","nB","sfB"));
dfrCmpWnN  = dfrCmpW |> dplyr::filter(nN_new !=nN_old);  #--204 rows
dfrCmpWsfN = dfrCmpW |> dplyr::filter(sfN_new!=sfN_old); #--204 rows
dfrCmpWnB  = dfrCmpW |> dplyr::filter(nB_new !=nB_old);  #--466 rows
dfrCmpWsfB = dfrCmpW |> dplyr::filter(sfB_new!=sfB_old); #--466 rows

dfrNewMnPropsAllYears = wtsUtilities::getObj(file.path(dirThs,"rda_Step1c_dfrMnPropsAllSBSYears.RData"));
dfrOldMnPropsAllYears = wtsUtilities::getObj(file.path(dirThs,"rda_Step1c_dfrMnPropsAllSBSYears.20260111.RData"));
dfrCmpW = dplyr::bind_rows(dfrNewMnPropsAllYears |> dplyr::mutate(case="new"),
                          dfrOldMnPropsAllYears |> dplyr::mutate(case="old")) |> 
            dplyr::select(case,y=YEAR,x=SEX,z=SIZE,numTot,mnPropNMFS,mnCPUE_NMFS,mnCPUE_BSFRF,mnLnR,mnCPUE_Prop,mnCPUE_lnR) |> 
          tidyr::pivot_wider(names_from="case",values_from=c("numTot","mnPropNMFS","mnCPUE_NMFS","mnCPUE_BSFRF","mnLnR","mnCPUE_Prop","mnCPUE_lnR"));
dfrCmpWnumTot      = dfrCmpW |> dplyr::filter(numTot_new      !=numTot_old);      #-- 153 rows
dfrCmpWmnPropNMFS  = dfrCmpW |> dplyr::filter(mnPropNMFS_new  !=mnPropNMFS_old);  #-- 116 rows
dfrCmpWmLnR        = dfrCmpW |> dplyr::filter(mnLnR_new       !=mnLnR_old);       #-- 116 rows
dfrCmpWmnCPUE_Prop = dfrCmpW |> dplyr::filter(mnCPUE_Prop_new !=mnCPUE_Prop_old); #-- 267 rows
dfrCmpWmnCPUE_lnR  = dfrCmpW |> dplyr::filter(mnCPUE_lnR_new  !=mnCPUE_lnR_old);  #-- 264 rows
ggplot(dfrCmpWmLnR,aes(x=z,shape=as.factor(y)))       + geom_point(aes(y=mnLnR_new),colour="green")      + geom_point(aes(y=mnLnR_old),colour="blue");
ggplot(dfrCmpWmnCPUE_Prop,aes(x=z,shape=as.factor(y))) + geom_point(aes(y=mnCPUE_Prop_new),colour="green") + geom_point(aes(y=mnCPUE_Prop_old),colour="blue");
ggplot(dfrCmpWmnCPUE_lnR,aes(x=z,shape=as.factor(y))) + geom_point(aes(y=mnCPUE_lnR_new),colour="green") + geom_point(aes(y=mnCPUE_lnR_old),colour="blue");

