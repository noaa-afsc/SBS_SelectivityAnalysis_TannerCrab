#--Extract sediment values at NMFS *SBS* haul locations
require(tcsamSurveyData);

dirThs = dirname((rstudioapi::getActiveDocumentContext())$path);

#--get interpolation rasters
top<-"./Rasters";
fns<-c("CoKr_EditedGrainSizeRaster.tif",
       "CoKr_EditedSortingRaster.tif");
rasters<-list();
for (fn in fns){
  fnp<-stringr::str_remove(fn,(".tif"))
  rasters[[fnp]]<-stars::read_stars(file.path(top,fn));
  plot(rasters[[fnp]]);
}
crs = sf::st_crs(rasters[[1]]);

#--get NMFS SBS data
lst = wtsUtilities::getObj(file.path(dirThs,"../../00a_SBS_Data/rda_Step1_SBS_RawData.RData"));
dfrHD_NMFS_SBS = lst$dfrHD_NMFS_SBS;

#--create sf dataset
dfrHD_NMFS_SBS = wtsGIS::createSF_points(dfrHD_NMFS_SBS,xCol="MID_LONGITUDE",yCol="MID_LATITUDE");
dfrHD_NMFS_SBS = sf::st_transform(dfrHD_NMFS_SBS,crs);#--transform to raster CRS (Alaska Albers)

#--extract values from rasters
library(ggplot2);
dfrHD = dfrHD_NMFS_SBS;
dfrHD = cbind(dfrHD,phi    =stars::st_extract(rasters[[1]],dfrHD)[[1]]);
dfrHD = cbind(dfrHD,sorting=stars::st_extract(rasters[[2]],dfrHD)[[1]]);
wtsUtilities::saveObj(dfrHD,
                      file.path(dirThs,"rda_dfrHD_NMFS_SBS_WithInterpolatedSedValues.RData"));

#--plot normalized temp, phi, and sorting as a function of depth
dfrHDp = dfrHD %>% sf::st_drop_geometry() %>%
                   dplyr::select(BOTTOM_DEPTH,GEAR_TEMPERATURE,phi,sorting) %>%
                   dplyr::rename(depth=BOTTOM_DEPTH,tmp=GEAR_TEMPERATURE,srt=sorting) %>%
                   dplyr::mutate(dtmp=(tmp-min(tmp,na.rm=TRUE))/(max(tmp,na.rm=TRUE)-min(tmp,na.rm=TRUE)),
                                 dphi=(phi-min(phi,na.rm=TRUE))/(max(phi,na.rm=TRUE)-min(phi,na.rm=TRUE)),
                                 dsrt=(srt-min(srt,na.rm=TRUE))/(max(srt,na.rm=TRUE)-min(srt,na.rm=TRUE))) %>%
                   dplyr::select(depth,dtmp,dphi,dsrt) %>%
                   dplyr::rename(`bottom temperature`=dtmp,phi=dphi,sorting=dsrt);
dfrHDpp = dfrHDp %>%
             tidyr::pivot_longer(c(`bottom temperature`,phi,sorting),
                                 names_to="variable",
                                 values_to="values");
p_vsd = ggplot2::ggplot(dfrHDpp,mapping=aes(x=depth,y=values,colour=variable,fill=variable)) +
          ggplot2::geom_point(alpha=0.2) +
          ggplot2::geom_smooth(alpha=0.4) + 
          ggplot2::labs(x="bottom depth (m)",y="normalized value",colour="quantity",fill="quantity") +
          ggplot2::theme(legend.position="bottom");
print(p_vsd);
ggsave(file.path(dirThs,"Figures/fig03_ColinearityCheck.png"),p_vsd,
       device="png",width=7.75,height=6.5,units="in");

#--make pairs plots 
p_col = GGally::ggpairs(dfrHDp)
ggsave(file.path(dirThs,"Figures/fig03_ColinearityCheck.png"),p_vsd,
       device="png",width=7.75,height=6.5,units="in");
