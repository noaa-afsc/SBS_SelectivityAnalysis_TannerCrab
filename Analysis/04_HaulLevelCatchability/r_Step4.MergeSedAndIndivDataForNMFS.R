#--get indiv data with sediment characteristics for all 83-112 NMFS hauls----
  dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

  #--get NMFS individual data----
  require(tcsamSurveyData);
  dirData_NMFS    = "~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
  fnStrata        = file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
  fnHaulData_NMFS = file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
  verbosity       = 0;
  minYr = 1982;
  maxYr = 2024;
  
  ##--read NMFS strata definitions and create NMFS strata table for SBS years----
  dfrStrata<-readr::read_csv(fnStrata) |>
               dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,
                             LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM) |> 
               dplyr::filter(dplyr::between(SURVEY_YEAR,minYr,maxYr));
  dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                  species="BTC",
                                  strataType="2015",
                                  export=FALSE,
                                  verbosity=verbosity);
  rm(dfrStrata);
  ##--NMFS haul data----
  ##----read in NMFS haul data 
  dfrHaulData_NMFS<-readr::read_csv(file=fnHaulData_NMFS,
                                    guess_max=Inf,
                                    skip=5) |> 
                      dplyr::filter(dplyr::between(AKFIN_SURVEY_YEAR,minYr,maxYr));
  ##------extract NMFS haul-specific data tbl
  dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                      tbl=dfrHaulData_NMFS,
                                      YearRange=c(minYr,maxYr),
                                      export=FALSE,
                                      verbosity=verbosity);
  ##--select NMFS individual data at haul stations
  sex   = "ALL";
  mat   =  "ALL";
  s_c   =  "ALL";
  minZ  =  25;
  maxZ  =  Inf;
  calcMaleMat = FALSE;
  dfrID_NMFS<-selectIndivs.TrawlSurvey(dfrHD_NMFS,
                                       tbl=dfrHaulData_NMFS,
                                       sex=sex,
                                       maturity=mat,
                                       shell_condition=s_c,
                                       calcMaleMaturity=calcMaleMat,
                                       minSize=minZ,
                                       maxSize=maxZ,
                                       export=FALSE,
                                       verbosity=verbosity) |> 
                dplyr::filter(SEX %in% c("MALE","FEMALE"));
  
  #--calculate CPUE by size bin by haul
  dfrCPUE = calcCPUE.ByHaul(dfrHD_NMFS,
                            dfrID_NMFS,
                            bySex=TRUE,
                            byMaturity=FALSE,
                            byShellCondition=FALSE,
                            bySize=TRUE,
                            cutpts=seq(25,185,5)) |> 
              dplyr::mutate(SIZE=SIZE+2.5);#--shift to bin center
  wtsUtilities::saveObj(dfrCPUE,file.path(dirThs,"rda_Step4.dfrCPUE_NMFS.RData"));
