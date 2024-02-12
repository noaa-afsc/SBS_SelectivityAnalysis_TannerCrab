#--do bootstrapping for empirical selectivity
require(tcsamSurveyData);

dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

source(file.path(dirThs,"r_extractHaulsAndIndivs.R"));
source(file.path(dirThs,"r_calcEmpiricalSelectivity.R"));

dfn = file.path(dirPrj,"Analysis/00a_SBS_Data/rda_Step1_SBS_RawData.RData");#--SBS data list object
lst = wtsUtilities::getObj(dfn);
dfrSD_SBS       = lst$dfrSD_SBS;
dfrHD_BSFRF_SBS = lst$dfrHD_BSFRF_SBS;
dfrID_BSFRF_SBS = lst$dfrID_BSFRF_SBS;
dfrHD_NMFS_SBS  = lst$dfrHD_NMFS_SBS;
dfrID_NMFS_SBS  = lst$dfrID_NMFS_SBS;

#--set up processing parameters
aggBySex           =FALSE;
aggByMaturity      =TRUE;
aggByShellCondition=TRUE;
cutpts             =seq(from=25,to=185,by=5);
truncate.low       =TRUE;
truncate.high      =FALSE;
verbosity          =0;

#--determine unique years
uYs<-unique(dfrSD_SBS$YEAR);

#--do bootstrapping
dfrZCs<-NULL;
dfrESs<-NULL;
nBs<-1000;
for (iY in uYs){
  #--get actual SBS stations, hauls, and indiv data for year iY
  dfrYSD_SBS<-dfrSD_SBS[dfrSD_SBS$YEAR==iY,];
  dfrYHD_BSFRF_SBS<-dfrHD_BSFRF_SBS[dfrHD_BSFRF_SBS$YEAR==iY,];
  dfrYID_BSFRF_SBS<-dfrID_BSFRF_SBS[dfrID_BSFRF_SBS$HAULJOIN %in% dfrYHD_BSFRF_SBS$HAULJOIN,];
  dfrYHD_NMFS_SBS <-dfrHD_NMFS_SBS[dfrHD_NMFS_SBS$YEAR==iY,];
  dfrYID_NMFS_SBS <-dfrID_NMFS_SBS[dfrID_NMFS_SBS$HAULJOIN %in% dfrYHD_NMFS_SBS$HAULJOIN,];

  #--get number of stations/hauls in EBS for year iY
  nS_SBS <- nrow(dfrYSD_SBS);
  for (iB in 1:nBs){
    cat("Processing",iB,"of",nBs,"for",iY,"\n");
    #--resample from SBS stations
    idS_SBS <- ceiling(stats::runif(n=nS_SBS,min=0,max=nS_SBS));#random index to row of dfrSDy_SBS
    dfrRSD_SBS<-dfrYSD_SBS[idS_SBS,];#--resampled SBS stations

    #----resample BSFRF hauls and individuals
    lst<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_BSFRF_SBS,dfrYID_BSFRF_SBS,resampleIndivs=TRUE);
    dfrRSD_BSFRF_SBS<-lst$dfrSDr;#--rows have been reordered in ascending order
    dfrRHD_BSFRF_SBS<-lst$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
    dfrRID_BSFRF_SBS<-lst$dfrIDr;#--newHAULJOINs match those in dfrRHD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
    rm(lst);
    #----resample NMFS hauls and individuals
    lst<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_NMFS_SBS,dfrYID_NMFS_SBS,resampleIndivs=TRUE);
    dfrRSD_NMFS_SBS<-lst$dfrSDr;#--rows have been reordered in ascending order
    dfrRHD_NMFS_SBS<-lst$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
    dfrRID_NMFS_SBS<-lst$dfrIDr;#--newHAULJOINs match those in dfrRHD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
    rm(lst);

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
    lst<-calcEmpiricalSelectivity(dfrRSD_SBS,
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
    lst$dfrZCs[["iB"]]<-iB;
    lst$dfrESs[["iB"]]<-iB;
    dfrZCs<-rbind(dfrZCs,lst$dfrZCs);
    dfrESs<-rbind(dfrESs,lst$dfrESs);
    rm(lst);
  }#--iB

}#--iY

save(dfrZCs,dfrESs,
     file=file.path(dirThs,"rda_Step2a_EmpiricalSelectivityFromBootstrapping.RData"));


