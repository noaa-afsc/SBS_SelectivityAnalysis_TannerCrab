getNumbersCaught_NMFS<-function(yr_rng=c(1982,2023),
                                sex="MALE",
                                sz_rng=c(15,154),
                                cutpts = seq(-0.5,204.5,5)){
  minYr = yr_rng[1];
  maxYr = yr_rng[2];
  #--NMFS data
  dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
  fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
  fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
  #--read NMFS strata definitions
  dfrStrata<-readr::read_csv(fnStrata) |>
               dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,
                             LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM) |> 
               dplyr::filter(dplyr::between(SURVEY_YEAR,minYr,maxYr));
  dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                  species="BTC",
                                  strataType="2015",
                                  export=FALSE,
                                  verbosity=0) |> 
           dplyr::filter(dplyr::between(YEAR,minYr,maxYr));
  rm(dfrStrata);
  #--read in NMFS haul data corresponding to stations at which BSFRF SBS hauls were conducted
  #----read in haul data file
  dfrHaulData_NMFS<-readr::read_csv(file=fnHaulData_NMFS,
                                    guess_max=Inf,
                                    skip=5) |> 
                      dplyr::filter(dplyr::between(AKFIN_SURVEY_YEAR,minYr,maxYr));
  #----Select NMFS haul data
  dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                      tbl=dfrHaulData_NMFS,
                                      YearRange=c(minYr,maxYr),
                                      export=FALSE,
                                      verbosity=0);
  dfrID_NMFS<-selectIndivs.TrawlSurvey(dfrHD_NMFS,
                                       tbl=dfrHaulData_NMFS,
                                       sex=sex,
                                       maturity="ALL",
                                       shell_condition="ALL",
                                       calcMaleMaturity=FALSE,
                                       minSize=sz_rng[1],
                                       maxSize=sz_rng[2],
                                       export=FALSE,
                                       verbosity=0);
  rm(dfrHaulData_NMFS);
  #--calculate numbers-at-size catch by haul
  dfrZCsByH = dfrHD_NMFS |> 
                dplyr::inner_join(dfrID_NMFS |> 
                                    dplyr::mutate(bin=cutpts[cut(SIZE,cutpts,labels=FALSE)]+3),
                                  by="HAULJOIN") |> 
               dplyr::group_by(YEAR,HAULJOIN,bin) |> 
               dplyr::summarize(numIndivs=sum(numIndivs,na.rm=TRUE)) |> 
               dplyr::ungroup();
  return(dfrZCsByH);
}
