#--Extract SBS data for analysis
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

#--load necessary libraries
require(tcsamSurveyData);

#--processing parameters
minYr = 2013;
maxYr = 2018;
sex   = "ALL";
mat   =  "ALL";
s_c   =  "ALL";
minZ  =  0;
maxZ  =  Inf;
calcMaleMat = FALSE;
verbosity   = 0;

#--NMFS data
dirData_NMFS   <-"~/Work/StockAssessments-Crab/Data/Survey.NMFS.EBS/Current";
fnStrata       <-file.path(dirData_NMFS,"TannerCrab_SurveyStrata.csv");
fnHaulData_NMFS<-file.path(dirData_NMFS,"TannerCrab_HaulData.csv");
#--BSFRF data
dirData_BSFRF   <-"~/Work/StockAssessments-Crab/Data/Survey.BSFRF/AllSurveys";
fnHaulData_BSFRF<-file.path(dirData_BSFRF,"rda_BSFRF.SBS.TannerCrabHaulData.RData");

#--read NMFS strata definitions
dfrStrata<-readr::read_csv(fnStrata) |>
             dplyr::select(STATION_ID,DISTRICT,TOWS,SURVEY_YEAR,
                           LATITUDE,LONGITUDE,STRATUM,TOTAL_AREA_SQ_NM) |> 
             dplyr::filter(dplyr::between(SURVEY_YEAR,minYr,maxYr));
dfrSD<-selectStrata.TrawlSurvey(dfrStrata,
                                species="BTC",
                                strataType="2015",
                                export=FALSE,
                                verbosity=verbosity) |> 
         dplyr::filter(dplyr::between(YEAR,minYr,maxYr));
rm(dfrStrata);

#--read in NMFS haul data corresponding to stations at which BSFRF SBS hauls were conducted
#----read in haul data file
dfrHaulData_NMFS<-readr::read_csv(file=fnHaulData_NMFS,
                                  guess_max=Inf,
                                  skip=5) |> 
                    dplyr::filter(dplyr::between(AKFIN_SURVEY_YEAR,minYr,maxYr));
#----Select NMFS haul data for SBS survey years
dfrHD_NMFS<-selectHauls.TrawlSurvey(dfrSD,
                                    tbl=dfrHaulData_NMFS,
                                    YearRange=c(minYr,maxYr),
                                    export=FALSE,
                                    verbosity=verbosity);

#--read in BSFRF SBS haul data
dfrHaulData_BSFRF = wtsUtilities::getObj(fnHaulData_BSFRF);
dfrHD_chk = dplyr::distinct(dfrHaulData_BSFRF,AKFIN_SURVEY_YEAR,HAULJOIN,GIS_STATION);#--390 hauls
#----select unique BSFRF SBS hauls that match a NMFS station
#----can match multiple hauls at a single station
dfrHD_BSFRF = selectHauls.TrawlSurvey(dfrSD,
                                      tbl=dfrHaulData_BSFRF,
                                      YearRange=c(minYr,maxYr),
                                      export=FALSE,
                                      verbosity=verbosity);
#--selects 390 hauls that match YEAR, GIS_STATION in dfrSD

#--NMFS: check sampling factors for uniqueness with hauls
dfr = dfrHaulData_NMFS |> dplyr::select(year=AKFIN_SURVEY_YEAR,hj=HAULJOIN,
                                        stn=GIS_STATION,asw=AREA_SWEPT,sub=SUBSAMPLE,
                                        x=SEX,z=WIDTH,sf=SAMPLING_FACTOR);
#----extract distinct sampling factors by haul
dfr1 = dfr |> dplyr::distinct(sf)
ggplot(dfr1,aes(x=sf)) + geom_histogram()
#----extract hauls with sampling factors > 1
dfr1 = dfr |> dplyr::distinct(year,hj,x,sf) |> 
              dplyr::count(year,hj,x,name="n") |> 
              dplyr::filter(n>1);
dfr2 = dfr |> dplyr::inner_join(dfr1,by=c("year","hj","x")) |>
              dplyr::mutate(sfr=1/sf,
                            cnt=1) |> 
              dplyr::group_by(year,hj,x,sf) |> 
              dplyr::summarize(tsf=sum(sf,na.rm=TRUE),  #--expanded total number
                               tct=sum(cnt,na.rm=TRUE));#--total number measured

#----NMFS
dfrID_NMFS<-selectIndivs.TrawlSurvey(dfrHD_NMFS,
                                         tbl=dfrHaulData_NMFS,
                                         sex=sex,
                                         maturity=mat,
                                         shell_condition=s_c,
                                         calcMaleMaturity=calcMaleMat,
                                         minSize=minZ,
                                         maxSize=maxZ,
                                         export=FALSE,
                                         verbosity=verbosity);
#----check weights
#------check individual weights using WatZ regressions as implemented in tcsamFunctions::calc.WatZ
dfrID_NMFS$CALC_WGT1<-dfrID_NMFS$numIndivs*
                                   tcsamFunctions::calc.WatZ(dfrID_NMFS$SIZE,
                                                             dfrID_NMFS$SEX,
                                                             dfrID_NMFS$MATURITY);
chkWgts1   = sd(dfrID_NMFS$CALCULATED_WEIGHT - dfrID_NMFS$CALC_WGT1,na.rm=TRUE);
dfrID_chk1 = dfrID_NMFS |> dplyr::filter(abs(CALCULATED_WEIGHT - dfrID_NMFS$CALC_WGT1) > 0.000001)
#------check individual weights using WatZ regressions as implemented by NMFS survey
calcWatZ2<-function(z,sex,clutch_size,
                    male=list(all=list(a=0.00027,b=3.022134)),
                    female=list(`1`=list(a=0.000562,b=2.816928),
                                `2`=list(a=0.000441,b=2.898686))){
    idx.m<-toupper(sex)=='MALE';
    idx.f<-toupper(sex)=='FEMALE';
    idx.c<-clutch_size > 1;
    wgt<-(male$all$a     *z^male$all$b   )*(idx.m)+
         (female[[1]]$a  *z^female[[1]]$b)*(idx.f&(!idx.c))+
         (female[[2]]$a  *z^female[[2]]$b)*(idx.f&( idx.c));
    
    return(wgt);
}
dfrID_NMFS$CALC_WGT2<-dfrID_NMFS$numIndivs*
                                   calcWatZ2(dfrID_NMFS$SIZE,
                                             dfrID_NMFS$SEX,
                                             dfrID_NMFS$CLUTCH_SIZE);
chkWgts2   = sd(dfrID_NMFS$CALCULATED_WEIGHT - dfrID_NMFS$CALC_WGT2,na.rm=TRUE);
dfrID_chk2 = dfrID_NMFS |> dplyr::filter(abs(CALCULATED_WEIGHT - dfrID_NMFS$CALC_WGT2) > 0.000001)

#--estimate observation variance for sampling factors > 1
dfrHDp_NMFS = dfrHD_NMFS |> dplyr::filter(HAULJOIN %in% dfr1$hj);
nHs = nrow(dfrHDp_NMFS);
nBHs=5000;#--number of bootstrapped hauls
nBZs=1;   #--number of bootstrapped size comps for each bootstrapped haul
b = 1;
lstBHs = list();
while (b <= nBHs){
  #--testing b = 1;
  irw  = ceiling(runif(1,0,nHs));    #--index of selected haul
  dfrH = dfrHDp_NMFS[irw,];          #--selected haul
  hj_  = dfrH$HAULJOIN;
  dfrUXSs = dfrID_NMFS |> dplyr::filter(HAULJOIN==hj_) |> 
             dplyr::distinct(SEX,SAMPLING_FACTOR) |> 
             dplyr::filter(SAMPLING_FACTOR>1);
  lstXS = list();
  for (irwXS in 1:nrow(dfrUXSs)){
  #--testing: irwXS=1;
    dfrUXS = dfrUXSs[irwXS,];
    ux  = dfrUXS$SEX;
    usf = dfrUXS$SAMPLING_FACTOR;
  #--calculate total sex-specific weight and abundance for haul
    totWA = dfrID_NMFS |> dplyr::filter(HAULJOIN==hj_,SEX==ux) |> 
               dplyr::summarize(wgt=sum(SAMPLING_FACTOR*CALCULATED_WEIGHT,na.rm=TRUE),
                                abd=sum(SAMPLING_FACTOR,na.rm=TRUE)); 
    #--select only sub-sampled crab of sex ux, sampling factor usf
    dfrI = dfrID_NMFS |> dplyr::filter(HAULJOIN==hj_,SEX==ux,SAMPLING_FACTOR==usf);
    nIs = nrow(dfrI);
    if (nIs>0){
      #--calculate "true" sex-specific basket weight for portion of haul to be subsampled
      wgtTrSSBW_ = (dfrI |> 
                     dplyr::summarize(wgt=sum(SAMPLING_FACTOR*CALCULATED_WEIGHT,na.rm=TRUE)))$wgt;

      #--start bootstrap results with "true" case
      lstBZs = list();
      #--check actual sampling factor calculation
      dfrRS = dfrI;#--original measured individuals
      #--calculated weight of (bootstrapped) measured individuals in haul
      wgtBSS_  = (dfrRS |> 
                   dplyr::summarize(wgt=sum(CALCULATED_WEIGHT,na.rm=TRUE)))$wgt;
      lstBZs[[1]] = tibble::tibble(bh=b,                      #--haul sampling bootstrap counter
                                   bhxi=0,                    #--individual (within sex) sampling bootstrap counter
                                   hj=hj_,                    #--haul join
                                   x=ux,                      #--sex
                                   wgtT=totWA$wgt,            #--reported total sex-specific weight for haul
                                   abdT=totWA$abd,            #--reported total abundance
                                   n=nIs,                     #--number of  individuals actually measured in subsample
                                   wgtTrSSBW=wgtTrSSBW_,      #--"true" sex-specific basket weight for individuals to be subsampled
                                   sfTr=usf,                  #--"true" sampling factor
                                   wgtBSS=wgtBSS_,            #--weight for bootstrapped subsampled individuals
                                   sfB=wgtTrSSBW_/wgtBSS_,    #--sampling factor from (non) bootstrapped individuals (should match "true" sampling factor)
                                   abdTB=abdT-nIs*(sfTr-sfB));#--"new" estimate of total abundance (matches reported since no resampling)
      #--bootstrap by resampling measured individuals to get bootstrapped sampling factors
      for (iBZ in 1:nBZs){
        idx_rs = floor(runif(nIs,0,nIs))+1;  #--indices to resample individuals
        dfrRS  = dfrI[idx_rs,];              #--"new" observation of size composition
        #--calculate weight of measured individuals
        wgtBSS_  = (dfrRS |> 
                      dplyr::summarize(wgt=sum(CALCULATED_WEIGHT,na.rm=TRUE)))$wgt;
        lstBZs[[1+iBZ]] = tibble::tibble(bh=b,                      #--haul sampling bootstrap counter
                                         bhxi=iBZ,                  #--individual (within sex) sampling bootstrap counter
                                         hj=hj_,                    #--hauljoin
                                         x=ux,                      #--sex
                                         wgtT=totWA$wgt,             #--reported total sex-specific weight for haul
                                         abdT=totWA$abd,            #--reported total abundance
                                         n=nIs,                     #--number of  individuals actually measured in subsample
                                         wgtTrSSBW=wgtTrSSBW_,      #--"true" sex-specific basket weight for individuals to be subsampled
                                         sfTr=usf,                  #--"true" sampling factor
                                         wgtBSS=wgtBSS_,            #--weight for bootstrapped subsampled individuals
                                         sfB=wgtTrSSBW_/wgtBSS_,    #--sampling factor from bootstrapped individuals
                                         abdTB=abdT-nIs*(sfTr-sfB));#--estimated total abundance for bootstrap realization
      }#--iBZ
      lstXS[[paste(ux,usf)]] = dplyr::bind_rows(lstBZs); rm(lstBZs);
    }#--nIs>0
  }#--irwXS
  lstBHs[[b]] = dplyr::bind_rows(lstXS);  rm(lstXS);
  b = b+1;
}#--b
dfrBHs = dplyr::bind_rows(lstBHs);
if (FALSE) wtsUtilities::saveObj(dfrBHs,file.path(dirThs,"rda_BootstrappedSFs.RData"));

#--analyze & plot results
dfrStats = dfrBHs |> dplyr::filter(bhxi>0,n>5) |> 
                     dplyr::group_by(x) |> 
                     dplyr::summarize(relErrMnSF=mean((sfTr-sfB)/sfTr,na.rm=TRUE),
                                      relErrSdSF=  sd((sfTr-sfB)/sfTr,na.rm=TRUE),
                                      relErrMnTA=mean((abdT-abdTB)/abdT,na.rm=TRUE),
                                      relErrSdTA=  sd((abdT-abdTB)/abdT,na.rm=TRUE));

dfr = dfrBHs |> 
        dplyr::filter(bhxi>0,n>5) |> 
        dplyr::mutate(relErrSF=(sfTr-sfB)/sfTr,
                      relErrTA=(abdT-abdTB)/abdT);
ggplot(dfr,aes(x=relErrSF,y=after_stat(ndensity),colour=x,fill=x)) + 
  geom_histogram(binwidth=0.05,position="identity",alpha=0.5) + 
  labs(x="relative variation in sampling factor",
       y="normalized density",
       colour="sex",fill="sex") + 
  wtsPlots::getStdTheme();
ggplot(dfr,aes(x=relErrTA,y=after_stat(ndensity),colour=x,fill=x)) + 
  geom_histogram(binwidth=0.05,position="identity",alpha=0.5) + 
  labs(x="relative variation in total abundance",
       y="normalized density",
       colour="sex",fill="sex") + 
  wtsPlots::getStdTheme();

