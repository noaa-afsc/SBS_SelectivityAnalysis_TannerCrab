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

#--read in "CVR" () feature classes and calculate rms values
require(sf)
fc<-list();
fc[["CVR_DWB_MeanPhi"]] <- sf::st_read(dsn, layer = "CVR_DWB_MeanPhi");
fc[["CVR_EBK_MeanPhi"]] <- sf::st_read(dsn, layer = "CVR_EBK_MeanPhi");
fc[["CVR_DWB_Sorting"]] <- sf::st_read(dsn, layer = "CVR_DWB_Sorting");
fc[["CVR_EBK_Sorting"]] <- sf::st_read(dsn, layer = "CVR_EBK_Sorting");

nfc<-length(fc);
rms<-vector("numeric",nfc);
for (ifc in 1:nfc){
  rms[ifc]<-sqrt(sum(fc[[ifc]]$Error^2)/(nrow(fc[[ifc]])-1));
}
dfr_stats<-data.frame(type=names(fc),rms=rms,stringsAsFactors = FALSE);
View(dfr_stats);

#--make qq plots
require(ggplot2)
for (ifc in 1:nfc){
  dfr<-fc[[ifc]];
  p <-ggplot(dfr,mapping=aes(sample=Error));
  p <- p + stat_qq() + stat_qq_line();
  p <- p + labs(subtitle = names(fc)[ifc]);
  print(p);
}


