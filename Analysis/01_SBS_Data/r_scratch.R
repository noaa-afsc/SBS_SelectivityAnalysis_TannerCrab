dfrStrata |> dplyr::distinct(SURVEY_YEAR,DISTRICT,TOTAL_AREA_SQ_NM) |> 
             tidyr::pivot_wider(names_from="SURVEY_YEAR",values_from="TOTAL_AREA_SQ_NM")

dfrSD |> dplyr::group_by(YEAR,STRATUM,STRATUM_AREA_BYSTATION) |> 
         dplyr::summarize(tot_area=wtsUtilities::Sum(STATION_AREA)) |> 
         dplyr::ungroup() |> 
         dplyr::mutate(tot_area=1.852*1.852*tot_area)

grids$grid$STATION_AREA = grids$grid |> sf::st_area() |> units::set_units(nmile*nmile);
dfrSD |> dplyr::inner_join(grids$grid |> sf::st_drop_geometry(),
                           by=c("GIS_STATION"="STATION_ID")) |>

#-------------------------------------------------------------------------------                             
  #----get the survey grid layers (grid and stations)
  grid = tcsamSurveyData::gisGetSurveyGridLayers()$grid |>
           dplyr::transmute(AREA=TOTAL_AREA,STATION_ID=STATION_ID);
  dfrUniqStns<-unique(dfrSD[,c("YEAR","STRATUM","GIS_STATION")]);
  polysUniqStns <- wtsGIS::mergeDataframeWithLayer(dfrUniqStns,
                                                   grid,
                                                   dataID="GIS_STATION",
                                                   geomsID="STATION_ID");

  #----calculate the area of each stratum by summing over the area associated with each station
  #-----NOTE: STATION_AREA, STRATUM_AREA_BYSTATION will be in square nautical miles
  if ("TOTAL_AREA" %in% names(polysUniqStns)) polysUniqStns[["AREA"]] = polysUniqStns[["TOTAL_AREA"]];
  if ("TOTAL_AREA" %in% names(polysUniqStns)) polysUniqStns[["AREA"]] = polysUniqStns[["TOTAL_AREA"]];
  tmp1<-polysUniqStns[,c("YEAR","STRATUM","GIS_STATION","AREA"),drop=TRUE];#keep some columns, drop geometry
  tmp1$STATION_AREA <- tmp1$AREA/(1852*1852);#convert to sq. nm.
  qry<-"select YEAR,STRATUM,
        sum(STATION_AREA) as STRATUM_AREA_BYSTATION
        from tmp1
        group by YEAR,STRATUM
        order by YEAR,STRATUM;";
  tmp2<-sqldf::sqldf(qry);
  qry<-"select t1.YEAR,
               t1.STRATUM,
               t1.GIS_STATION,
               t1.STATION_AREA,
               t2.STRATUM_AREA_BYSTATION
        from tmp1 as t1, tmp2 as t2
        where t1.YEAR=t2.YEAR and t1.STRATUM=t2.STRATUM;"
  tmp3<-sqldf::sqldf(qry);

  #----add station area and stratum area based on sum over station areas
  qry<-"select
          s.YEAR,s.STRATUM,s.STRATUM_CODE,s.STRATUM_AREA,
          s.GIS_STATION,s.STATION_LONGITUDE,s.STATION_LATITUDE,
          t.STATION_AREA as STATION_AREA,t.STRATUM_AREA_BYSTATION as STRATUM_AREA_BYSTATION
        from dfrSD as s, tmp3 as t
        where s.YEAR=t.YEAR and s.GIS_STATION=t.GIS_STATION;";
  dfr<-sqldf::sqldf(qry);
                            