#--set working directory to file location

#--checking r_Step1
dirThs = "."; #<-make sure to setwd to r_Step1a file folder first
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step1a_SBS_DataWithSedData.RData"));
View(lst$dfrHD_NMFS_SBS)
dfrChk = lst$dfrHD_NMFS_SBS |> dplyr::group_by(YEAR,HAULJOIN) |> dplyr::summarize(nc=dplyr::n()) |> dplyr::ungroup() |> dplyr::filter(nc>1)
nrow(dfrChk)  #--0
View(lst$dfrHD_BSFRF_SBS)
dfrChk = lst$dfrHD_BSFRF_SBS |> dplyr::group_by(YEAR,HAULJOIN) |> dplyr::summarize(nc=dplyr::n()) |> dplyr::ungroup() |> dplyr::filter(nc>1)
nrow(dfrChk)  #--0

dirThs = "."; #<-make sure to setwd to r_Step1b file folder first
dfrPropsAllYears   = wtsUtilities::getObj(file.path(dirThs,"rda_Step1b_dfrPropsAllSBSYears.RData"));
View(dfrPropsAllYears)
nrow(dfrPropsAllYears)
dfrChk = dfrPropsAllYears |> dplyr::group_by(YEAR,HAULJOIN,SEX,SIZE) |> dplyr::summarize(nc=dplyr::n()) |> dplyr::ungroup() |> dplyr::filter(nc>1)
nrow(dfrChk) #--now (2026-01-08) 0 after change to scaling numIndiv by sampling factor before aggregation 
             #--was 496 for BBRKC!  134 for Tanner crab! <==problem is multiple sampling factors involved (worse for BBRKC because of aggregation over sex)
dfrChkp = dfrChk |> dplyr::inner_join(dfrPropsAllYears,by=c("YEAR","HAULJOIN","SEX","SIZE"))

#--checking r_Step2...
dirThs = "."; #<-make sure to setwd to r_Step2 file folder first
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step2_TrimmedDataList.RData"));
View(lst$lstTrimmedFinal$dfrDat)
dfrChk = lst$lstTrimmedFinal$dfrDat |> group_by(y,h,z) |> dplyr::summarize(nc=dplyr::n()) |> dplyr::ungroup() |> dplyr::filter(nc>1)
nrow(dfrChk) #--now (2026-01-08) 0

#--checking r_Step3a
dirThs = "."; #<-make sure to setwd to r_Step3a file folder first
lst = wtsUtilities::getObj(file.path(dirThs,"../rda_Step2_TrimmedDataList.RData"));
dfrDat =  lst$lstTrimmedFinal$dfrDat |> 
            dplyr::mutate(obsR=exp(lnR),
                          nN=round(p*n), #--number caught in NMFS gear
                          nB=n-nN);      #--number caught in BSFRF gear
dfrChk = dfrDat |> group_by(y,h,z) |> dplyr::summarize(nc=dplyr::n()) |> dplyr::ungroup() |> dplyr::filter(nc>1)
nrow(dfrChk)  #--now (2026-01-08) 0
