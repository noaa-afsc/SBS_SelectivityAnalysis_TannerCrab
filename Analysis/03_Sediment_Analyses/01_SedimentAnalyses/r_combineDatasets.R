#--combine EBSSED-2 and bers_www datasets for sorting and mean phi
#----read in bers_WWD.csv and extract sorting, mean phi (grainsize)
#------note that missing values are coded "-99"
bers<-readr::read_csv("bers_WWD.csv");
names(bers)<-tolower(names(bers));
tbl_bers<-bers[,c("latitude","longitude","waterdepth","grainsze","sorting")];#--drop missing values
names(tbl_bers)<-c("latitude","longitude","depth","mean_phi","sorting");

#----extract sorting and mean phi datasets from EBSSED-2 ESRI file geodatabase
dsn <- "./EBS_Sed.gdb";             #--path to file geodatabase folder
lyrs  <- rgdal::ogrListLayers(dsn); #--layer (feature class) names in dsn
nlyrs <-length(lyrs);               #--number of layers
print(lyrs);
lyrs<-lyrs[6:7];  #--extract just sorting and mean phi layers
  for (lyr in lyrs){
    cat("Layer: '",lyr,"'\n",sep="")
    lyrInfo<-rgdal::ogrInfo(dsn=dsn,layer=lyr);
    print(lyrInfo);
  }
#------sorting
lyr_sorting <- sf::read_sf(dsn,lyrs[1])[,,drop=TRUE];#--read, but drop geometry
names(lyr_sorting)<-tolower(names(lyr_sorting));
tbl_sorting <- lyr_sorting[,c("latitude","longitude","sorting")];
#------mean phi
lyr_meanphi <- sf::read_sf(dsn,lyrs[2])[,,drop=TRUE];#--read, but drop geometry
names(lyr_meanphi)<-tolower(names(lyr_meanphi));
tbl_meanphi <- lyr_meanphi[,c("latitude","longitude","mean_phi")];

#----combine datasets
tbl_sorting<-rbind(tbl_sorting,tbl_bers[,c("latitude","longitude","sorting")]);
tbl_meanphi<-rbind(tbl_meanphi,tbl_bers[,c("latitude","longitude","mean_phi")]);

#--remove missing data (coded -99)
idx<-tbl_sorting$sorting  != -99; tbl_sorting<-tbl_sorting[idx,];
idx<-tbl_meanphi$mean_phi != -99; tbl_meanphi<-tbl_meanphi[idx,];

#----remove duplicates
tbl_sorting <- dplyr::distinct(tbl_sorting);
tbl_meanphi <- dplyr::distinct(tbl_meanphi);

#----export datasets
readr::write_csv(tbl_sorting,"combined_sorting.csv");
readr::write_csv(tbl_meanphi,"combined_meanphi.csv");

