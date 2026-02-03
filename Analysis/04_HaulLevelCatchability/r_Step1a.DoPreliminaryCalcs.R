#--do preliminary calculations for haul-level catchability analysis
##--interpolate sediment data to SBS haul data locations----

#--read project setup info----
dirPrj = rstudioapi::getActiveProject();
dirThs = file.path(dirPrj,"Analysis/04_HaulLevelCatchability")
fn = file.path(dirPrj,"rda_ProjectSetup.RData");
s  = wtsUtilities::getObj(fn);

#--load "raw" SBS haul data
lst = wtsUtilities::getObj(file.path(s$dirs$SBS_Data,"rda_Step1_SBS_RawData.RData"));

#--load sediment data (phi, sorting) interpolated to haul locations
dfrHD_sed = wtsUtilities::getObj(file.path(s$dirs$SedAnls,
                                           "rda_dfrHD_All_NMFS_WithInterpolatedSedValues.RData"));

#--join sediment data to NMFS SBS haul data by YEAR and GIS_STATION (equivalent to using HAULJOIN)
#----note: geometry has to be dropped because dfrHD_NMFS_SBS is a dataframe, not a tibble or sf object
dfrHD_NMFS_SBS  = lst$dfrHD_NMFS_SBS |>  
                    dplyr::inner_join((dfrHD_sed |> sf::st_drop_geometry() |>
                                          dplyr::select(YEAR,GIS_STATION,phi,sorting)),
                                       by=c("YEAR","GIS_STATION"));
lst$dfrHD_NMFS_SBS = dfrHD_NMFS_SBS;

#--join sediment data to BSFRF SBS haul data by YEAR and GIS_STATION (can't use HAULJOIN--not equivalent)
#----note: geometry has to be dropped because dfrHD_BSFRF_SBS is a dataframe, not a tibble or sf object
dfrHD_BSFRF_SBS = lst$dfrHD_BSFRF_SBS |>  
                    dplyr::inner_join((dfrHD_sed |> sf::st_drop_geometry() |>
                                          dplyr::select(YEAR,GIS_STATION,phi,sorting)),
                                       by=c("YEAR","GIS_STATION"));

lst$dfrHD_BSFRF_SBS = dfrHD_BSFRF_SBS;
##--save SBS data with interpolated sediment info----
wtsUtilities::saveObj(lst,file.path(dirThs,"rda_Step1a_SBS_DataWithSedData.RData"));
rm(dfrHD_sed,dfrHD_BSFRF_SBS,dfrHD_NMFS_SBS);

