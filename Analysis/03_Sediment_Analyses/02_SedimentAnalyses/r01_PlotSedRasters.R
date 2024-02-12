#--plot sediment rasters
require(ggplot2);
require(sf)
require(stars);

dirThs = dirname((rstudioapi::getActiveDocumentContext())$path);

#--get layers from SedimentsAnalysis.gdb ArcGIS file geodatabase
dsn  = file.path(dirThs,"../SedimentsAnalysis.gdb");

#--get info on the avaialble layers
lyrs = sf::st_layers(dsn);#--a list (with lists in some elements)
tbl_lyrs = tibble::tibble(name=lyrs$name,
                          geomtype=unlist(lyrs$geomtype),
                          driver=lyrs$driver,
                          features=lyrs$features,
                          fields=lyrs$fields);
#--View(tbl_lyrs);

#--read combined sediment data feature classes
fc<-list();
fc[["MeanPhi"]] <- sf::st_read(dsn,layer="CombinedSedData_CrabSurveyArea_MeanPhi");
fc[["Sorting"]] <- sf::st_read(dsn,layer="CombinedSedData_CrabSurveyArea_Sorting");

#--get interpolation rasters
top<-file.path(dirThs,"Rasters");
fns<-c("CoKr_EditedGrainSizeRaster.tif",
       "CoKr_EditedSortingRaster.tif");
rasters<-list();
for (fn in fns){
  fnp<-stringr::str_remove(fn,(".tif"))
  rasters[[fnp]]<-stars::read_stars(file.path(top,fn));
  plot(rasters[[fnp]]);
}
crs = sf::st_crs(rasters[[1]]);
bbox = sf::st_bbox(rasters[[1]]);

bmls = wtsGIS::gg_CreateBasemapLayers(final.crs=crs,bbox=bbox);
basemap = ggplot()+bmls$land+bmls$theme+bmls$map_scale;
sgls = tcsamSurveyData::gisGetSurveyGridLayers();

#--map phi data
p_phi1 = ggplot()+
        geom_sf(data=fc[[1]],mapping=aes(fill=mean_phi),shape=21)+
        bmls$land +
        geom_sf(data=sgls$grid,fill=NA,color="black") +
        scale_fill_viridis_c(option="plasma",limits=c(-2,8),oob=scales::squish)+
        labs(fill="phi") +
        bmls$map_scale +
        bmls$theme;
#--map phi interpolation raster
p_phi2 = ggplot()+
        geom_stars(data=rasters[[1]])+
#        geom_sf(data=dfrHD,mapping=aes_string(fill="phi"),shape=21)+
        bmls$land +
        geom_sf(data=sgls$grid,fill=NA,color="black") +
        scale_fill_viridis_c(option="plasma",limits=c(-2,8),oob=scales::squish)+
        labs(fill="phi") +
        bmls$map_scale +
        bmls$theme;
pg = ggpubr::ggarrange(p_phi1,p_phi2,ncol=1,common.legend=TRUE,legend="right");
ggsave("./Figures/fig01_MapPhi.png",pg,device="png",width=6.5,height=6.5,units="in");

#--map sorting data
p_srt1 = ggplot()+
        geom_sf(data=fc[[2]],mapping=aes_string(fill="sorting"),shape=21)+
        bmls$land +
        geom_sf(data=sgls$grid,fill=NA,color="black") +
        scale_fill_viridis_c(option="plasma",limits=c(0,3),oob=scales::squish)+
        labs(fill="sorting") +
        bmls$map_scale +
        bmls$theme;
#--map sorting interpolation raster
p_srt2 = ggplot()+
        geom_stars(data=rasters[[2]])+
#        geom_sf(data=dfrHD,mapping=aes_string(fill="phi"),shape=21)+
        bmls$land +
        geom_sf(data=sgls$grid,fill=NA,color="black") +
        scale_fill_viridis_c(option="plasma",limits=c(0,3),oob=scales::squish)+
        labs(fill="sorting") +
        bmls$map_scale +
        bmls$theme;
pg = ggpubr::ggarrange(p_srt1,p_srt2,ncol=1,common.legend=TRUE,legend="right");
ggsave("./Figures/fig02_MapSorting.png",pg,device="png",width=6.5,height=6.5,units="in");
