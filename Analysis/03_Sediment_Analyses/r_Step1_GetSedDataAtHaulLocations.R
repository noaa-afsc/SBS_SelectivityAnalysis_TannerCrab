#--extract sediment characteristics at all NMFS haul locations
getSedData<-function(){
  out = list();
  
  #--get interpolation rasters
  top<-file.path("02_SedimentAnalyses","Rasters");
  fns<-c("CoKr_EditedGrainSizeRaster.tif",
         "CoKr_EditedSortingRaster.tif");
  rasters<-list();
  for (fn in fns){
    fnp<-stringr::str_remove(fn,(".tif"))
    rasters[[fnp]]<-stars::read_stars(file.path(top,fn));
  }
  crs = sf::st_crs(rasters[[1]]);
  
  #--create ggpolot2 basemap layers
  bmls = wtsGIS::gg_CreateBasemapLayers();
  out = c(out,list(bmls=bmls));

  #--get **all NMFS haul** data
  minYr = 1975;
  maxYr = 2023;
  verbosity   = 0;

  #--NMFS data
  dirData_NMFS    = "~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
  fnStrata        = file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
  fnHaulData_NMFS = file.path(dirData_NMFS,"TannerCrab_HaulData.csv");

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

  #--create sf dataset
  dfrHD_NMFS = wtsGIS::createSF_points(dfrHD_NMFS,xCol="MID_LONGITUDE",yCol="MID_LATITUDE");
  dfrHD_NMFS = sf::st_transform(dfrHD_NMFS,crs);#--transform to raster CRS (Alaska Albers)
  
  #--extract values from rasters
  dfrHD = dfrHD_NMFS;
  dfrHD = cbind(dfrHD,phi    =stars::st_extract(rasters[[1]],dfrHD)[[1]]);
  dfrHD = cbind(dfrHD,sorting=stars::st_extract(rasters[[2]],dfrHD)[[1]]);
  out = c(out,list(dfrHD=dfrHD));
  
  wtsUtilities::saveObj(out,"rda_Step1_SedDataAtAllNMFSHaulLocations.RData");
  rm(dfrHaulData_NMFS, dfrHD_NMFS, dfrSD, dirData_NMFS, fnHaulData_NMFS, fnStrata);
  
  return(out);
}
out = getSedData();
rm(out);

