#--calculate empirical selectivity from data
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
#--define output filename
ofn = file.path(dirThs,"rda_Step1_EmpiricalSelectivityFromData.RData");

require(tcsamSurveyData);
source(file.path(dirThs,"r_CalcEmpiricalSelectivity.R"));
source(file.path(dirThs,"r_ExtractHaulsAndIndivs.R"));

dfn = file.path(dirPrj,"Analysis/00a_SBS_Data/rda_Step1_SBS_RawData.RData");#--SBS data list object
lst = wtsUtilities::getObj(dfn);
#--set up processing parameters
aggBySex           =FALSE;
aggByMaturity      =TRUE;
aggByShellCondition=TRUE;
cutpts             =seq(from=25,to=185,by=5);
truncate.low       =TRUE;
truncate.high      =FALSE;
verbosity          =0;

#--determine unique years
uYs<-unique(lst$dfrSD_SBS$YEAR);

#--calculate size compositions
dfrZCs<-NULL;
dfrESs<-NULL;
nBs<-0;
for (iY in uYs){
  #--get actual SBS stations, hauls, and indiv data for year iY
  #--iY = uYs[1];
  dfrYSD_SBS = lst$dfrSD_SBS |> dplyr::filter(YEAR==iY);
  dfrYHD_BSFRF_SBS = lst$dfrHD_BSFRF_SBS |> dplyr::filter(YEAR==iY);
  dfrYID_BSFRF_SBS = lst$dfrID_BSFRF_SBS |> dplyr::filter(HAULJOIN %in% dfrYHD_BSFRF_SBS$HAULJOIN);
  dfrYHD_NMFS_SBS  = lst$dfrHD_NMFS_SBS |> dplyr::filter(YEAR==iY);
  dfrYID_NMFS_SBS  = lst$dfrID_NMFS_SBS |> dplyr::filter(HAULJOIN %in% dfrYHD_NMFS_SBS$HAULJOIN);

  #--get number of stations/hauls in EBS for year iY
  nS_SBS <- nrow(dfrYSD_SBS);
  for (iB in 0){
    cat("Processing",iB,"of",nBs,"for",iY,"\n");
    #--DON'T resample from SBS stations
    # idS_SBS <- ceiling(stats::runif(n=nS_SBS,min=0,max=nS_SBS));#random index to row of dfrSDy_SBS
    # dfrRSD_SBS<-dfrYSD_SBS[idS_SBS,];#--resampled SBS stations
    dfrRSD_SBS<-dfrYSD_SBS;#--use original SBS stations (will be resampled stations in step3a)
    #----extract BSFRF hauls and individuals (don't really need to do this, but consistent with step3a)
    lst1<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_BSFRF_SBS,dfrYID_BSFRF_SBS,resampleIndivs=FALSE);
    dfrRSD_BSFRF_SBS<-lst1$dfrSDr;#--rows have been reordered in ascending order
    dfrRHD_BSFRF_SBS<-lst1$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
    dfrRID_BSFRF_SBS<-lst1$dfrIDr;#--newHAULJOINs match those in dfrRHD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
    rm(lst1);
    #----DON'T resample NMFS hauls and individuals
    lst1<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_NMFS_SBS,dfrYID_NMFS_SBS,resampleIndivs=FALSE);
    dfrRSD_NMFS_SBS<-lst1$dfrSDr;#--rows have been reordered in ascending order
    dfrRHD_NMFS_SBS<-lst1$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
    dfrRID_NMFS_SBS<-lst1$dfrIDr;#--newHAULJOINs match those in dfrRHD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
    rm(lst1);

    #--dfrRSD_BSFRF_SBS and dfrRSD_NMFS_SBS are identical
    #--one needs to replace dfrRSD_SBS so re-ordering is correct
    dfrRSD_SBS<-dfrRSD_BSFRF_SBS;
    rm(dfrRSD_BSFRF_SBS,dfrRSD_NMFS_SBS);

    #--now replace HAULJOINS with newHAULJOINS
    dfrRHD_BSFRF_SBS$HAULJOIN<-dfrRHD_BSFRF_SBS$newHAULJOIN; dfrRHD_BSFRF_SBS<-wtsUtilities::deleteCol(dfrRHD_BSFRF_SBS,"newHAULJOINS");
    dfrRID_BSFRF_SBS$HAULJOIN<-dfrRID_BSFRF_SBS$newHAULJOIN; dfrRID_BSFRF_SBS<-wtsUtilities::deleteCol(dfrRID_BSFRF_SBS,"newHAULJOINS");
    dfrRHD_NMFS_SBS$HAULJOIN <-dfrRHD_NMFS_SBS$newHAULJOIN;  dfrRHD_NMFS_SBS <-wtsUtilities::deleteCol(dfrRHD_NMFS_SBS, "newHAULJOINS");
    dfrRID_NMFS_SBS$HAULJOIN <-dfrRID_NMFS_SBS$newHAULJOIN;  dfrRID_NMFS_SBS <-wtsUtilities::deleteCol(dfrRID_NMFS_SBS, "newHAULJOINS");
    #--rename GIS_STATIONS to avoid duplicates
    dfrRSD_SBS$GIS_STATION      <-as.character(1:nS_SBS);
    dfrRHD_BSFRF_SBS$GIS_STATION<-as.character(1:nS_SBS);
    dfrRHD_NMFS_SBS$GIS_STATION <-as.character(1:nS_SBS);

    #--create single stratum for SBS stations
    dfrRSD_SBS$STRATUM<-"1";
    dfrRSD_SBS$STRATUM_CODE<-"1";
    dfrRSD_SBS$STRATUM_AREA          <-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
    dfrRSD_SBS$STRATUM_AREA_BYSTATION<-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
    dfrRHD_BSFRF_SBS$STRATUM<-"1";
    dfrRHD_NMFS_SBS$STRATUM <-"1";

    #--calculate empirical selectivity
    lst1<-calcEmpiricalSelectivity(dfrRSD_SBS,
                                  dfrRHD_BSFRF_SBS,
                                  dfrRID_BSFRF_SBS,
                                  dfrRHD_NMFS_SBS,
                                  dfrRID_NMFS_SBS,
                                  aggBySex=aggBySex,
                                  aggByMaturity=aggByMaturity,
                                  aggByShellCondition=aggByShellCondition,
                                  cutpts=cutpts,
                                  truncate.low=truncate.low,
                                  truncate.high=truncate.high,
                                  showPlot=FALSE,
                                  verbosity=0);
    lst1$dfrZCs[["iB"]]<-iB;
    lst1$dfrESs[["iB"]]<-iB;
    dfrZCs<-rbind(dfrZCs,lst1$dfrZCs);
    dfrESs<-rbind(dfrESs,lst1$dfrESs);
    rm(lst1);
  }#--iB

}#--iY

wtsUtilities::saveObj(list(dfrZCs=dfrZCs,dfrESs=dfrESs),fn=ofn);


