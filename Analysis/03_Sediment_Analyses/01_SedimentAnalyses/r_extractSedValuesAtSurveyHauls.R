#--get interpolation rasters
top<-"./Rasters";
fns<-c("MeanPhi_DWB_CrabSurveyArea1.tif",
       "MeanPhi_EBK_CrabSurveyArea1.tif",
       "Sorting_DWB_CrabSurveyArea3.tif",
       "Sorting_EBK_CrabSurveyArea1.tif");
rasters<-list();
for (fn in fns){
  fnp<-stringr::str_remove(fn,(".tif"))
  rasters[[fnp]]<-stars::read_stars(file.path(top,fn));
  plot(rasters[[fnp]]);
}

#--get NMFS data
require(tcsamSurveyData);
minYr=2017;
maxYr=2017;
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
#----read NMFS strata definitions
dfrStrata<-read.csv(fnStrata,check.names=FALSE,stringsAsFactors=FALSE)[,1:8];
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="BTC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity);
#--read in NMFS haul data
dfrHaulData_NMFS<-read.csv(file=fnHaulData_NMFS,check.names=FALSE,stringsAsFactors=FALSE);
dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                    tbl=dfrHaulData_NMFS,
                                    YearRange=c(minYr,maxYr),
                                    export=FALSE,
                                    verbosity=verbosity);
rm(dfrHaulData_NMFS,dfrStrata);

