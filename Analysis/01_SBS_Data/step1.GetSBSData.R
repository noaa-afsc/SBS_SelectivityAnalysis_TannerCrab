#--Extract SBS data for analysis
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

#--load necessary libraries
require(tcsamSurveyData);

#--processing parameters
minYr = 2013;
maxYr = 2018;
sex   = "ALL";
mat   =  "ALL";
s_c   =  "ALL";
minZ  =  0;
maxZ  =  Inf;
calcMaleMat = FALSE;
verbosity   = 0;

#--NMFS data
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
#--BSFRF data
dirData_BSFRF   <-"~/Work/StockAssessments-Crab/Data/Survey.BSFRF/AllSurveys";
fnHaulData_BSFRF<-file.path(dirData_BSFRF,"rda_BSFRF.SBS.TannerCrabHaulData.RData");

#--read NMFS strata definitions
dfrStrata<-readr::read_csv(fnStrata) |>
             dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,
                           LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM) |> 
             dplyr::filter(dplyr::between(SURVEY_YEAR,minYr,maxYr));
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="BTC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity) |> 
         dplyr::filter(dplyr::between(YEAR,minYr,maxYr));
rm(dfrStrata);

#--read in NMFS haul data corresponding to stations at which BSFRF SBS hauls were conducted
#----read in haul data file
dfrHaulData_NMFS<-readr::read_csv(file=fnHaulData_NMFS,
                                  guess_max=Inf,
                                  skip=5) |> 
                    dplyr::filter(dplyr::between(AKFIN_SURVEY_YEAR,minYr,maxYr));
#----Select NMFS haul data for SBS survey years
dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                    tbl=dfrHaulData_NMFS,
                                    YearRange=c(minYr,maxYr),
                                    export=FALSE,
                                    verbosity=verbosity);

#--read in BSFRF SBS haul data
dfrHaulData_BSFRF = wtsUtilities::getObj(fnHaulData_BSFRF);
dfrHD_chk = dplyr::distinct(dfrHaulData_BSFRF,AKFIN_SURVEY_YEAR,HAULJOIN,GIS_STATION);#--390 hauls
#----select unique BSFRF SBS hauls that match a NMFS station
#----can match multiple hauls at a single station
dfrHD_BSFRF = selectHauls.TrawlSurvey(dfrSD,
                                      tbl=dfrHaulData_BSFRF,
                                      YearRange=c(minYr,maxYr),
                                      export=FALSE,
                                      verbosity=verbosity);
#--selects 390 hauls that match YEAR, GIS_STATION in dfrSD

#--the following should work for all hauls to filter out duplicates at a station
dfrHD_BSFRF_SBS = dfrHD_BSFRF |> dplyr::inner_join(dfrHD_NMFS |> dplyr::select(YEAR,GIS_STATION,START_DATE,NMFS_HOUR=START_HOUR),
                                                   by=c("YEAR","GIS_STATION","START_DATE")) |>
                                 dplyr::filter(abs(as.numeric(START_HOUR)-as.numeric(NMFS_HOUR))==
                                               min(abs(as.numeric(START_HOUR)-as.numeric(NMFS_HOUR))),
                                               .by=c(YEAR,GIS_STATION,START_DATE));
#----the above results in 384 hauls, not 390 hauls as expected, so some BSFRF hauls that were labeled SBS must not 
#------occur on the same day as the NMFS haul at the same station
dfrDropped = dfrHD_BSFRF |> dplyr::anti_join(dfrHD_BSFRF_SBS,by="HAULJOIN")
dfrDropped |> dplyr::select(YEAR,GIS_STATION,START_DATE,START_HOUR,HAULJOIN);
#   YEAR GIS_STATION START_DATE START_HOUR                  HAULJOIN
# 1 2015        G-09   06062015       0640  SBS;OH3;06062015;0640;20
# 2 2016        C-04   06202016       1326 SBS;SSB;06202016;1326;105
# 3 2016        F-08   06062016       0712  SBS;SSB;06062016;0712;32
# 4 2016        J-03   06092016       0749  SBS;SSB;06092016;0749;42
# 5 2016        J-03   06092016       1016  SBS;SSB;06092016;1016;43
# 6 2017        D-02   07122017       0703  SBS;HMB;07122017;0703;41

#--now need to select NMFS strata and hauls that match the selected BSFRF hauls
dfrSD_SBS = dfrSD |> dplyr::inner_join(dfrHD_BSFRF_SBS |> dplyr::select(YEAR,GIS_STATION),
                                       by=c("YEAR","GIS_STATION"));
dfrHD_NMFS_SBS = dfrHD_NMFS |> dplyr::inner_join(dfrHD_BSFRF_SBS |> dplyr::select(YEAR,GIS_STATION),
                                                 by=c("YEAR","GIS_STATION"));
#----the above results in 384 hauls, as well

#--cleanup a bit
rm(dfrSD,dfrHD_BSFRF,dfrHD_NMFS,dfrHD_chk);

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
lstAll = list(minYr=minYr,
              maxYr=maxYr,
              sex=sex,
              maturity=mat,
              shel_condition=s_c,
              calcMaleMaturity=calcMaleMat,
              minZ=minZ,
              maxZ=maxZ,
              dfrSD_SBS=dfrSD_SBS,
              dfrHD_BSFRF_SBS=dfrHD_BSFRF_SBS,
              dfrID_BSFRF_SBS=dfrID_BSFRF_SBS,
              dfrHD_NMFS_SBS=dfrHD_NMFS_SBS,
              dfrID_NMFS_SBS=dfrID_NMFS_SBS);
wtsUtilities::saveObj(lstAll,file.path(dirThs,"rda_step1_SBS_RawData.RData"));

