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

#--read in required layers
fc_CVR_DWB_MeanPhi <- sf::st_read(dsn, layer = "CVR_DWB_MeanPhi");
fc_CVR_EBK_MeanPhi <- sf::st_read(dsn, layer = "CVR_EBK_MeanPhi");
fc_CVR_DWB_Sorting <- sf::st_read(dsn, layer = "CVR_DWB_Sorting");
fc_CVR_EBK_Sorting <- sf::st_read(dsn, layer = "CVR_EBK_Sorting");

