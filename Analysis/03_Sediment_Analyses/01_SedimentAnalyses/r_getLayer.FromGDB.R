require(rgdal)

# The input file geodatabase
dsn <- "./SedimentsAnalysis.gdb"

#--check rgdal can handle gdb files
subset(ogrDrivers(), grepl("GDB", name))
#--list all feature classes in a file geodatabase
lyrs  <- ogrListLayers(dsn)
nlyrs <-length(lyrs);
print(lyrs);
if (nlyrs<0){ #--won't happen
  for (lyr in lyrs){
    cat("Layer: '",lyr,"'\n",sep="")
    lyrInfo<-rgdal::ogrInfo(dsn=dsn,layer=lyr);
    print(lyrInfo);
  }
}

require(sf)
fc_CVR_DWB_MeanPhi <- sf::st_read(dsn, layer = "CVR_DWB_MeanPhi");
fc_CVR_EBK_MeanPhi <- sf::st_read(dsn, layer = "CVR_EBK_MeanPhi");
fc_CVR_DWB_Sorting <- sf::st_read(dsn, layer = "CVR_DWB_Sorting");
fc_CVR_EBK_Sorting <- sf::st_read(dsn, layer = "CVR_EBK_Sorting");

