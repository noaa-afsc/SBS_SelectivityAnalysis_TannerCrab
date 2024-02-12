#--get layers from SedimentsAnalysis.gdb ArcGIS file geodatabase
require(sf);

dirThs = dirname((rstudioapi::getActiveDocumentContext())$path);

# The input file geodatabase
dsn  = file.path(dirThs,"../SedimentsAnalysis.gdb");

#--get info on the avaialble layers
lyrs = sf::st_layers(dsn);#--a list (with lists in some elements)
tbl_lyrs = tibble::tibble(name=lyrs$name,
                          geomtype=unlist(lyrs$geomtype),
                          driver=lyrs$driver,
                          features=lyrs$features,
                          fields=lyrs$fields);
View(tbl_lyrs);

#--read in "CVR" () feature classes and calculate rms values
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


