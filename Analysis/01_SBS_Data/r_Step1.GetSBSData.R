#--Extract SBS data for analysis
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

#--SBS data
dirData_BSFRF   <-"~/Work/StockAssessments-Crab/Data/Survey.BSFRF/AllSurveys_TannerCrab";
fnHaulData_SBS<-file.path(dirData_BSFRF,"rda_SBS.TannerCrab.DataObjects.RData");
lstAll = wtsUtilities::getObj(fnHaulData_SBS);
wtsUtilities::saveObj(lstAll,file.path(dirThs,"rda_step1_SBS_RawData.RData"));

