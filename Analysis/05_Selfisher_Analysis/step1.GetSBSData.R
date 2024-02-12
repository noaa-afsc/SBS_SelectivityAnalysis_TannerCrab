#--This script saves SBS data objects in the following file
output_rdataFile<-"data.SBS_Data.RData";#file name to save SBS data objects to

#--load necessary libraries
require(tcsamSurveyData);

minYr<-2013;
maxYr<-2017;
sex <- "ALL";
mat <- "ALL";
s_c <- "ALL";
minZ <- 0;
maxZ <- Inf;
calcMaleMat<- FALSE;
verbosity<-0;

#--NMFS data
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
#--BSFRF data
dirData_BSFRF   <-"~/Work/StockAssessments-Crab/Data/Survey.BSFRF/2013-2017";
fnHaulData_BSFRF<-file.path(dirData_BSFRF,"Stockhausen_BSFRF_Data_2013-2017.MS_reviewed.AllHauls.csv");

#--read NMFS strata definitions
dfrStrata<-read.csv(fnStrata,check.names=FALSE,stringsAsFactors=FALSE)[,1:8];
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="BTC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity);
rm(dfrStrata);

#--read in BSFRF SBS data
tmp<-read.csv(file=fnHaulData_BSFRF,check.names=FALSE,stringsAsFactors=FALSE);
dfrHaulData_BSFRF <- convertFormat.BSFRF2NMFS(tmp,types="SBS",verbosity=verbosity);
#------correct some known errors in 2013-2017 data
dfrHaulData_BSFRF<-dfrHaulData_BSFRF[!(dfrHaulData_BSFRF$GIS_STATION %in% c("AZ-0405","AZ0405")),];
dfrHaulData_BSFRF$GIS_STATION[dfrHaulData_BSFRF$GIS_STATION=="GF-1918"]<-"GF1918";
dfrHaulData_BSFRF$GIS_STATION[dfrHaulData_BSFRF$GIS_STATION=="HG-1918"]<-"HG1918";
dfrHaulData_BSFRF$GIS_STATION[dfrHaulData_BSFRF$GIS_STATION=="IH-1918"]<-"IH1918";
dfrHaulData_BSFRF$GIS_STATION[dfrHaulData_BSFRF$GIS_STATION=="JI-1918"]<-"JI1918";
dfrHaulData_BSFRF$GIS_STATION[dfrHaulData_BSFRF$GIS_STATION=="H21"]    <-"H-21";
#----select haul data (selects 467 hauls)
dfrHD_BSFRF<-selectHauls.TrawlSurvey(dfrSD,
                                     tbl=dfrHaulData_BSFRF,
                                     YearRange=c(minYr,maxYr),
                                     export=FALSE,
                                     verbosity=verbosity);
rm(tmp);
#------remove hauls with unknown area swept
idx<-is.na(dfrHD_BSFRF$AREA_SWEPT_VARIABLE);
dfrHD_BSFRF<-dfrHD_BSFRF[!idx,];#removes 145 hauls, keeps 322 hauls
rm(idx);

#--get unique SBS stations by year and start date
dfrUniqBSFRFStns<-unique(dfrHD_BSFRF[,c("YEAR","STRATUM","GIS_STATION")]);               #320
dfrUniqStationsSBS<-unique(dfrHD_BSFRF[,c("YEAR","STRATUM","GIS_STATION","START_DATE")]);#321 rows

#--read in NMFS haul data corresponding to stations at which BSFRF SBS hauls were conducted
dfrHaulData_NMFS<-read.csv(file=fnHaulData_NMFS,check.names=FALSE,stringsAsFactors=FALSE);
#----Select haul data
#------Need to select only NMFS hauls corresponding to unique BSFRF stations (320)
#------Only matches on YEAR and GIS_STATION (not START_DATE)
dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrUniqBSFRFStns,
                                    tbl=dfrHaulData_NMFS,
                                    YearRange=c(minYr,maxYr),
                                    export=FALSE,
                                    verbosity=verbosity);#320 rows

#--Make sure stations in BSFRF and NMFS haul data match
#--in terms of YEAR, GIS_STATION, *AND* START_DATE
#----Re-select BSFRF stations to match subset of NMFS stations
qry<-"select b.YEAR,b.STRATUM,b.GIS_STATION,b.START_DATE
      from dfrUniqStationsSBS as b, dfrHD_NMFS as n
      where b.YEAR        = n.YEAR        and
            b.GIS_STATION = n.GIS_STATION and
            b.START_DATE  = n.START_DATE;";
dfrUniqStationsSBS<-sqldf::sqldf(qry); #316 rows

#----re-select BSFRF haul data to match selected stations
qry<-"select b.YEAR,b.STRATUM,b.GIS_STATION,b.HAULJOIN,b.HAUL_TYPE,
             b.START_DATE,b.MID_LATITUDE,b.MID_LONGITUDE,
             b.BOTTOM_DEPTH,b.GEAR_TEMPERATURE,b.AREA_SWEPT_VARIABLE
      from dfrHD_BSFRF as b, dfrUniqStationsSBS as u
      where b.YEAR        = u.YEAR        and
            b.GIS_STATION = u.GIS_STATION and
            b.START_DATE  = u.START_DATE;";
dfrHD_BSFRF_SBS<-sqldf::sqldf(qry);#316 rows

#--check BSFRF stations with multiple hauls
#----NOTE: these hauls occurred on same date as a NMFS haul at each station
qry<-"select YEAR, STRATUM, GIS_STATION, count(*) as HAUL_COUNT
      from dfrHD_BSFRF_SBS
      group by YEAR, STRATUM, GIS_STATION
      order by YEAR, STRATUM, GIS_STATION;"
dfrHaulCnts_SBS <- sqldf::sqldf(qry);
numDuplicateHauls_SBS <- sum(dfrHaulCnts_SBS[dfrHaulCnts_SBS$HAUL_COUNT>1,"HAUL_COUNT"]-1);#now 0

#----re-select NMFS haul data to match selected stations
dfrHD_NMFS_SBS<-selectHauls.TrawlSurvey(dfrUniqStationsSBS,
                                        tbl=dfrHaulData_NMFS,
                                        YearRange=c(minYr,maxYr),
                                        export=FALSE,
                                        verbosity=verbosity); #316 now

#--Create SBS strata based on jointly-sampled stations
#----get the survey grid layers (grid and stations)
surveyGridLayers <- tcsamSurveyData::gisGetSurveyGridLayers();
#----merge the BSFRF haul data with the station grid, keeping
#----only stations occupied by the BSFRF SBS study, to pick up
#----information on areal coverage for each station
polysUniqStns <- wtsGIS::mergeDataframeWithLayer(dfrUniqStationsSBS,
                                                 surveyGridLayers$grid,
                                                 dataID="GIS_STATION",
                                                 geomsID="STATION_ID",
                                                 sfJoinType="left join");
#----calculate the area of each stratum by summing over the area associated with each station
#-----NOTE: STRATUM_AREA will be in square nautical miles
tmp1<-polysUniqStns[,c("YEAR","STRATUM","GIS_STATION","AREA"),drop=TRUE];#keep some columns, drop geometry
qry<-"select
         YEAR,STRATUM,sum(AREA)/(1852*1852) as STRATUM_AREA
      from tmp1
      group by YEAR, STRATUM
      order by YEAR,STRATUM;";
tmp2 <- sqldf::sqldf(qry);
#----select strata info that corresponds to SBS stations
qry<-"select
        s.YEAR,s.STRATUM,s.STRATUM_CODE,s.STRATUM_AREA,
        s.GIS_STATION,s.STATION_LONGITUDE,s.STATION_LATITUDE,
        s.STATION_AREA
      from dfrSD as s, dfrUniqStationsSBS as u
      where s.YEAR=u.YEAR and s.GIS_STATION=u.GIS_STATION;";
tmp3<-sqldf::sqldf(qry);
#----replace STRATUM_AREA from full NMFS survey with that from SBS survey
qry<-"select
        s.YEAR,s.STRATUM,s.STRATUM_CODE,n.STRATUM_AREA,
        s.GIS_STATION,s.STATION_LONGITUDE,s.STATION_LATITUDE,
        s.STATION_AREA
      from tmp3 as s, tmp2 as n
      where s.YEAR=n.YEAR and s.STRATUM=n.STRATUM;";
dfrSD_SBS<-sqldf::sqldf(qry); #316 rows
rm(tmp1,tmp2,tmp3,dfrSD);

#--Select individual data
#----BSFRF
dfrID_BSFRF_SBS<-selectIndivs.TrawlSurvey(dfrHD_BSFRF_SBS,
                                          tbl=dfrHaulData_BSFRF,
                                          sex=sex,
                                          maturity=mat,
                                          shell_condition=s_c,
                                          calcMaleMaturity=calcMaleMat,
                                          minSize=minZ,
                                          maxSize=maxZ,
                                          export=FALSE,
                                          verbosity=verbosity);
#------add individual weights using WatZ regressions
dfrID_BSFRF_SBS$CALCULATED_WEIGHT<-dfrID_BSFRF_SBS$numIndivs*
                                   tcsamFunctions::calc.WatZ(dfrID_BSFRF_SBS$SIZE,
                                                             dfrID_BSFRF_SBS$SEX,
                                                             dfrID_BSFRF_SBS$MATURITY);
#----NMFS
dfrID_NMFS_SBS<-selectIndivs.TrawlSurvey(dfrHD_NMFS_SBS,
                                         tbl=dfrHaulData_NMFS,
                                         sex=sex,
                                         maturity=mat,
                                         shell_condition=s_c,
                                         calcMaleMaturity=calcMaleMat,
                                         minSize=minZ,
                                         maxSize=maxZ,
                                         export=FALSE,
                                         verbosity=verbosity);

save(dfrID_BSFRF_SBS,dfrID_NMFS_SBS,
     dfrHD_BSFRF_SBS,dfrHD_NMFS_SBS,
     dfrSD_SBS,
     minYr,maxYr,file=output_rdataFile);

objs<-ls();
idx<-objs %in% c("dfrID_BSFRF_SBS","dfrID_NMFS_SBS","dfrHD_BSFRF_SBS","dfrHD_NMFS_SBS","dfrSD_SBS","minYr","maxYr");
rm(list=objs[!idx]);#--clean up
rm(idx,objs);

