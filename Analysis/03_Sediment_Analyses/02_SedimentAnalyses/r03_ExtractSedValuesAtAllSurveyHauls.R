#--Extract sediment values at *all* NMFS survey haul locations
require(tcsamSurveyData);

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

#--get (all) NMFS data
verbosity=0;
minYr=1975;
maxYr=2023;
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
#----read NMFS strata definitions
dfrStrata<-readr::read_csv(fnStrata)[,1:8];
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="BTC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity);
#--read in NMFS haul data
dfrHaulData_NMFS<-readr::read_csv(file=fnHaulData_NMFS,skip=5);
dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                    tbl=dfrHaulData_NMFS,
                                    YearRange=c(minYr,maxYr),
                                    export=FALSE,
                                    verbosity=verbosity);
rm(dfrHaulData_NMFS,dfrStrata);
#--create sf dataset
dfrHD_NMFS = wtsGIS::createSF_points(dfrHD_NMFS,xCol="MID_LONGITUDE",yCol="MID_LATITUDE");
dfrHD_NMFS = sf::st_transform(dfrHD_NMFS,crs);#--transform to raster CRS (Alaska Albers)

#--extract values from rasters
library(ggplot2);
dfrHD = dfrHD_NMFS;
dfrHD = cbind(dfrHD,phi    =stars::st_extract(rasters[[1]],dfrHD)[[1]]);
dfrHD = cbind(dfrHD,sorting=stars::st_extract(rasters[[2]],dfrHD)[[1]]);
wtsUtilities::saveObj(dfrHD,"rda_dfrHD_NMFS_All_WithInterpolatedSedValues.RData")

#--plot normalized temp, phi, and sorting as a function of depth
dfrHDp = dfrHD %>% sf::st_drop_geometry() %>%
                   dplyr::select(BOTTOM_DEPTH,GEAR_TEMPERATURE,phi,sorting) %>%
                   dplyr::rename(depth=BOTTOM_DEPTH,tmp=GEAR_TEMPERATURE,srt=sorting) %>%
                   dplyr::mutate(dtmp=(tmp-min(tmp,na.rm=TRUE))/(max(tmp,na.rm=TRUE)-min(tmp,na.rm=TRUE)),
                                 dphi=(phi-min(phi,na.rm=TRUE))/(max(phi,na.rm=TRUE)-min(phi,na.rm=TRUE)),
                                 dsrt=(srt-min(srt,na.rm=TRUE))/(max(srt,na.rm=TRUE)-min(srt,na.rm=TRUE))) %>%
                   dplyr::select(depth,dtmp,dphi,dsrt) %>%
                   dplyr::rename(`bottom temperature`=dtmp,phi=dphi,sorting=dsrt) %>%
                   tidyr::pivot_longer(c(`bottom temperature`,phi,sorting),
                                       names_to="variable",
                                       values_to="values");
p_vsd = ggplot2::ggplot(dfrHDp,mapping=aes(x=depth,y=values,colour=variable,fill=variable)) +
          ggplot2::geom_point(alpha=0.2) +
          ggplot2::geom_smooth(alpha=0.4) + 
          ggplot2::labs(x="bottom depth (m)",y="normalized value",colour="quantity",fill="quantity") +
          ggplot2::theme(legend.position="bottom");
print(p_vsd);
ggsave("./Figures/fig03_ColinearityCheck.png",p_vsd,device="png",width=7.75,height=6.5,units="in");

#--make pairs plots 
dfrHDpp = dfrHD %>% sf::st_drop_geometry() %>%
                    dplyr::select(BOTTOM_DEPTH,GEAR_TEMPERATURE,phi,sorting) %>%
                    dplyr::rename(depth=BOTTOM_DEPTH,tmp=GEAR_TEMPERATURE,srt=sorting) %>%
                    dplyr::mutate(dtmp=(tmp-min(tmp,na.rm=TRUE))/(max(tmp,na.rm=TRUE)-min(tmp,na.rm=TRUE)),
                                  dphi=(phi-min(phi,na.rm=TRUE))/(max(phi,na.rm=TRUE)-min(phi,na.rm=TRUE)),
                                  dsrt=(srt-min(srt,na.rm=TRUE))/(max(srt,na.rm=TRUE)-min(srt,na.rm=TRUE))) %>%
                    dplyr::select(depth,dtmp,dphi,dsrt) %>%
                    dplyr::rename(`bottom temperature`=dtmp,phi=dphi,sorting=dsrt);
GGally::ggpairs(dfrHDpp)
