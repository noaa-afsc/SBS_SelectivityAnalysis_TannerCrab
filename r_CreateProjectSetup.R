#--create project setup object
dirPrj = rstudioapi::getActiveProject();
lstAll = list();
##--sub-folders
dirs = list();
adir = "Analysis";
dirs$SBS_Data = file.path(dirPrj,adir,"01_SBS_Data");
dirs$EmpSel   = file.path(dirPrj,adir,"02_Empirical_Selectivity");
dirs$SedAnls  = file.path(dirPrj,adir,"03_Sediment_Analyses");
dirs$HaulLvl  = file.path(dirPrj,adir,"04_HaulLevelCatchability");
lstAll$dirs = dirs;

wtsUtilities::saveObj(lstAll,file.path(dirPrj,"rda_ProjectSetup.RData"));
