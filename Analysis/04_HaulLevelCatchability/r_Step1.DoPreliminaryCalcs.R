#--do preliminary calculations for haul-level catchability analysis

#--read project setup info----
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
fn = file.path(dirPrj,"rda_ProjectSetup.RData");
s  = wtsUtilities::getObj(fn);

##--Function to extract  hauls and individuals from SBS+sediment data----
#'
#' @title Extract hauls and individuals from SBS+sediment data
#' 
#' @description Function to extract  hauls and individuals from SBS+sediment data.
#' 
#' @param dfrSDr - (possibly resampled) survey strata definitions
#' @param dfrHD - dataframe with haul data
#' @param dfrID - dataframe with individual data
#' @param resampledIndivs - flag (T/F) to resample individuals within hauls
#' 
#' @return dataframe
#' 
#' @details 
#' 
#' @export
#' 
extractHaulsAndIndivs<-function(dfrSDr,dfrHD,dfrID,resampleIndivs=FALSE){
  #--order resampled stations by GIS_STATION (could be duplicates)
  qry<- "select * from dfrSDr order by GIS_STATION;";
  dfrSDr<-sqldf::sqldf(qry);
  #--select resampled hauls
  hdvars<-paste0("h.",names(dfrHD),collapse=",");
  qry<-"select
          &&hdvars
        from dfrSDr as s, dfrHD as h
        where s.GIS_STATION=h.GIS_STATION
        order by h.GIS_STATION;"
  qry<-gsub("&&hdvars",hdvars,qry,fixed=TRUE);
  dfrHDr<-sqldf::sqldf(qry);
  #--rename GIS_STATIONs to avoid duplicates
  # stns<-as.character(1:nrow(dfrSDr));
  # dfrSDr$GIS_STATION<-stns;#--can't do this here: must do it outside function
  # dfrHDr$GIS_STATION<-stns;

  #--select resampled indivs
  dfrIDr<-NULL;
  nH<-nrow(dfrHDr);
  if (testing) cat("Resampling from",nH,"hauls\n");
  for (h in 1:nH){
    dfrIDp<-dfrID[dfrID$HAULJOIN==dfrHDr$HAULJOIN[h],];
    nI<-nrow(dfrIDp);
    if (nI>0){
      idI <- 1:nI;
      if (resampleIndivs) idI <- ceiling(stats::runif(n=nI,min=0,max=nI));#random index to row of dfrIDp
      dfrIDrp<-dfrIDp[idI,];    #extracted individuals
      dfrIDrp$newHAULJOIN<-h;   #append "new" HAULJOIN
      dfrIDr<-rbind(dfrIDr,dfrIDrp);
    }
  }
  dfrHDr$newHAULJOIN<-1:nrow(dfrHDr);#append "new" HAULJOIN
  return(list(dfrSDr=dfrSDr,dfrHDr=dfrHDr,dfrIDr=dfrIDr));
}

#--interpolate sediment data to SBS haul data locations----
#--load "raw" SBS haul data
lst = wtsUtilities::getObj(file.path(s$dirs$SBS_Data,"rda_Step1_SBS_RawData.RData"));

#--load sediment data (phi, sorting) interpolated to haul locations
dfrHD_sed = wtsUtilities::getObj(file.path(s$dirs$SedAnls,
                                           "rda_dfrHD_All_NMFS_WithInterpolatedSedValues.RData"));

#--join sediment data to NMFS SBS haul data by YEAR and GIS_STATION (equivalent to using HAULJOIN)
#----note: geometry has to be dropped because dfrHD_NMFS_SBS is a dataframe, not a tibble or sf object
dfrHD_NMFS_SBS  = lst$dfrHD_NMFS_SBS |>  
                    dplyr::inner_join((dfrHD_sed |> sf::st_drop_geometry() |>
                                          dplyr::select(YEAR,GIS_STATION,phi,sorting)),
                                       by=c("YEAR","GIS_STATION"));
lst$dfrHD_NMFS_SBS = dfrHD_NMFS_SBS;

#--join sediment data to BSFRF SBS haul data by YEAR and GIS_STATION (can't use HAULJOIN--not equivalent)
#----note: geometry has to be dropped because dfrHD_BSFRF_SBS is a dataframe, not a tibble or sf object
dfrHD_BSFRF_SBS = lst$dfrHD_BSFRF_SBS |>  
                    dplyr::inner_join((dfrHD_sed |> sf::st_drop_geometry() |>
                                          dplyr::select(YEAR,GIS_STATION,phi,sorting)),
                                       by=c("YEAR","GIS_STATION"));

lst$dfrHD_BSFRF_SBS = dfrHD_BSFRF_SBS;
##--save SBS data with interpolated sediment info----
wtsUtilities::saveObj(lst,file.path(dirThs,"rda_Step1a_SBS_DataWithSedData.RData"));
rm(dfrHD_sed,dfrHD_BSFRF_SBS,dfrHD_NMFS_SBS);

#--calculate haul-level size-specific empirical proportions from SBS + sed data----
dfrSD_SBS = lst$dfrSD_SBS;
dfrHD_BSFRF_SBS = lst$dfrHD_BSFRF_SBS;
dfrID_BSFRF_SBS = lst$dfrID_BSFRF_SBS;
dfrHD_NMFS_SBS  = lst$dfrHD_NMFS_SBS;
dfrID_NMFS_SBS  = lst$dfrID_NMFS_SBS;
rm(lst);

##--determine unique years----
uYs<-sort(unique(dfrSD_SBS$YEAR));

##--set up processing parameters----
aggBySex           =FALSE;
aggByMaturity      =TRUE;
aggByShellCondition=TRUE;
cutpts             =seq(from=5,to=185,by=5);
truncate.low       =TRUE;
truncate.high      =FALSE;
verbosity          =0;

nc<-length(cutpts);
delta<-cutpts[2]-cutpts[1];
sizebins<-(cutpts[2:nc]+cutpts[1:(nc-1)])/2;
minSize<-c(min(sizebins),min(sizebins)); names(minSize)<-c("MALE","FEMALE");
maxSize<-c(max(sizebins),132.5);         names(maxSize)<-c("MALE","FEMALE");
center<-as.numeric(c((maxSize["MALE"]+minSize["MALE"])/2,
                     (maxSize["FEMALE"]+minSize["FEMALE"])/2)); names(center)<-c("MALE","FEMALE");
width <-as.numeric(c((maxSize["MALE"]-minSize["MALE"]),  
                     (maxSize["FEMALE"]-minSize["FEMALE"])));   names(width) <-c("MALE","FEMALE");

##--drop "missing" sex data----
nmiss_BSFRF = dfrID_BSFRF_SBS |> dplyr::filter(SEX=="MISSING") |> nrow();
nmiss_NMFS  = dfrID_NMFS_SBS  |> dplyr::filter(SEX=="MISSING") |> nrow();
dfrID_BSFRF_SBS = dfrID_BSFRF_SBS |> dplyr::filter(SEX!="MISSING");
dfrID_NMFS_SBS  = dfrID_NMFS_SBS  |> dplyr::filter(SEX!="MISSING");

##--apply cut points to sizes----
dfrID_BSFRF_SBS$SIZE<-cutpts[cut(dfrID_BSFRF_SBS$SIZE,breaks=cutpts,include.lowest=TRUE)]+delta/2;
dfrID_NMFS_SBS$SIZE <-cutpts[cut(dfrID_NMFS_SBS$SIZE, breaks=cutpts,include.lowest=TRUE)]+delta/2;

##--calculate proportions----
testing = TRUE;
`%<>%` <- magrittr::`%<>%`;
lstPropsAllYears<-list();
for (iY in uYs){
  if (testing) cat("Processing",iY,"\n");
  #--get actual SBS stations, hauls, and indiv data for year iY
  dfrYSD_SBS<-dfrSD_SBS[dfrSD_SBS$YEAR==iY,];
  dfrYHD_BSFRF_SBS<-dfrHD_BSFRF_SBS[dfrHD_BSFRF_SBS$YEAR==iY,];
  dfrYID_BSFRF_SBS<-dfrID_BSFRF_SBS[dfrID_BSFRF_SBS$HAULJOIN %in% dfrYHD_BSFRF_SBS$HAULJOIN,];
  dfrYHD_NMFS_SBS <-dfrHD_NMFS_SBS[dfrHD_NMFS_SBS$YEAR==iY,];
  dfrYID_NMFS_SBS <-dfrID_NMFS_SBS[dfrID_NMFS_SBS$HAULJOIN %in% dfrYHD_NMFS_SBS$HAULJOIN,];

  #--get number of stations/hauls in EBS for year iY
  nS_SBS <- nrow(dfrYSD_SBS);
  #--DON'T resample from SBS stations
  # idS_SBS <- ceiling(stats::runif(n=nS_SBS,min=0,max=nS_SBS));#random index to row of dfrSDy_SBS
  # dfrRSD_SBS<-dfrYSD_SBS[idS_SBS,];#--resampled SBS stations
  dfrRSD_SBS<-dfrYSD_SBS;#--use original SBS stations (will be resampled stations in step3a)
  #----extract BSFRF hauls and individuals (don't really need to do this, but consistent with step3a)
  lst<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_BSFRF_SBS,dfrYID_BSFRF_SBS,resampleIndivs=FALSE);
  dfrRSD_BSFRF_SBS<-lst$dfrSDr;#--rows have been reordered in ascending order
  dfrRHD_BSFRF_SBS<-lst$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
  dfrRID_BSFRF_SBS<-lst$dfrIDr;#--newHAULJOINs match those in dfrRHD_BSFRF_SBS, newHAULJOINS must replace HAULJOINS below
  rm(lst);
  #----DON'T resample NMFS hauls and individuals
  lst<-extractHaulsAndIndivs(dfrRSD_SBS,dfrYHD_NMFS_SBS,dfrYID_NMFS_SBS,resampleIndivs=FALSE);
  dfrRSD_NMFS_SBS<-lst$dfrSDr;#--rows have been reordered in ascending order
  dfrRHD_NMFS_SBS<-lst$dfrHDr;#--newHAULJOINs indicate row in dfrRSD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
  dfrRID_NMFS_SBS<-lst$dfrIDr;#--newHAULJOINs match those in dfrRHD_NMFS_SBS, newHAULJOINS must replace HAULJOINS below
  rm(lst);

  #--dfrRSD_BSFRF_SBS and dfrRSD_NMFS_SBS are identical
  #--one needs to replace dfrRSD_SBS so re-ordering is correct
  dfrRSD_SBS<-dfrRSD_BSFRF_SBS;
  rm(dfrRSD_BSFRF_SBS,dfrRSD_NMFS_SBS);

  #--now replace HAULJOINS with newHAULJOINS
  dfrRHD_BSFRF_SBS %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) |> dplyr::select(-newHAULJOIN);
  dfrRID_BSFRF_SBS %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) |> dplyr::select(-newHAULJOIN);
  dfrRHD_NMFS_SBS  %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) |> dplyr::select(-newHAULJOIN);
  dfrRID_NMFS_SBS  %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) |> dplyr::select(-newHAULJOIN);
  #--rename GIS_STATIONS to avoid duplicates
  dfrRSD_SBS$GIS_STATION      <-as.character(1:nS_SBS);
  dfrRHD_BSFRF_SBS$GIS_STATION<-as.character(1:nS_SBS);
  dfrRHD_NMFS_SBS$GIS_STATION <-as.character(1:nS_SBS);

  #--create single stratum for SBS stations
  dfrRSD_SBS %<>% dplyr::mutate(STRATUM="1", STRATUM_CODE="1",
                                STRATUM_AREA=wtsUtilities::Sum(STATION_AREA),
                                STRATUM_AREA_BYSTATION=STRATUM_AREA);
  # dfrRSD_SBS$STRATUM<-"1";
  # dfrRSD_SBS$STRATUM_CODE<-"1";
  # dfrRSD_SBS$STRATUM_AREA          <-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
  # dfrRSD_SBS$STRATUM_AREA_BYSTATION<-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
  dfrRHD_BSFRF_SBS$STRATUM<-"1";
  dfrRHD_NMFS_SBS$STRATUM <-"1";

  #--recast data into form required by selfisher
  #----join area swept and sampling factor info to indiv data
  recastIndivData1<-function(idfr){
    qry<-"select
            HAULJOIN, SEX, SIZE, SAMPLING_FACTOR,
            sum(numIndivs) as numIndivs
          from idfr as i
          group by
            HAULJOIN, SEX, SIZE, SAMPLING_FACTOR;";
    res<-sqldf::sqldf(qry);
    return(res);
  }
  dfrRCI_NMFS <-recastIndivData1(dfrRID_NMFS_SBS);
  dfrRCI_BSFRF<-recastIndivData1(dfrRID_BSFRF_SBS);
  uHs<-unique(unique(dfrRCI_BSFRF$HAULJOIN),unique(dfrRCI_NMFS$SHAULJOIN));
  uXs<-unique(unique(dfrRCI_BSFRF$SEX),     unique(dfrRCI_NMFS$SEX));
  uZs<-unique(unique(dfrRCI_BSFRF$SIZE),    unique(dfrRCI_NMFS$SIZE));
  #--expand indiv info dataframes to all hauljoins, sexes, and sizes
  expandIndivData<-function(idfr,uHs,uXs,uZs){
    uHs<-data.frame(HAULJOIN=uHs,stringsAsFactors=FALSE);
    uXs<-data.frame(SEX=uXs,stringsAsFactors=FALSE);
    uZs<-data.frame(SIZE=uZs,stringsAsFactors=FALSE);
    qry<-"select HAULJOIN,SEX,SIZE from uHs,uXs,uZs;";
    uHXZs<-sqldf::sqldf(qry);
    qry<-"select
            u.HAULJOIN,u.SEX,u.SIZE,
            i.SAMPLING_FACTOR,i.numIndivs
          from uHXZs as u left join idfr as i
          on
            u.HAULJOIN=i.HAULJOIN and
            u.SEX=i.SEX and
            u.SIZE=i.SIZE;";
    tmp<-sqldf::sqldf(qry);
    idx<-is.na(tmp$numIndivs);
    tmp$SAMPLING_FACTOR[idx]<-1;
    tmp$numIndivs[idx]      <-0;
    return(tmp);
  }
  dfrREI_NMFS <-expandIndivData(dfrRCI_NMFS, uHs,uXs,uZs);
  dfrREI_BSFRF<-expandIndivData(dfrRCI_BSFRF,uHs,uXs,uZs);
  #--merge haul+sediment data with indiv data
  #----note: the expansion factors are defined such that 
  #        CPUE [number/(unit area)] = (number sampled)*expFactor
  #----and expFactor = SAMPLING_FACTOR/AREA_SWEPT
  #----note: the SAMPLING_FACTOR = 1/(sampling fraction) so
  #----    expFactor = 1/(AREA_SWEPT * (sampling fraction))
  mergeHaulData<-function(hdfr,idfr){
    qry<-"select
            h.HAULJOIN, BOTTOM_DEPTH, GEAR_TEMPERATURE, phi, sorting,
            AREA_SWEPT_VARIABLE,
            SEX, SIZE, SAMPLING_FACTOR,
            SAMPLING_FACTOR/AREA_SWEPT_VARIABLE as expFACTOR,
            numIndivs
          from hdfr as h, idfr as i
          where
            h.HAULJOIN = i.HAULJOIN";
    res<-sqldf::sqldf(qry);
    return(res);
  }
  dfrRHI_NMFS <-mergeHaulData(dfrRHD_NMFS_SBS, dfrREI_NMFS);
  dfrRHI_BSFRF<-mergeHaulData(dfrRHD_BSFRF_SBS,dfrREI_BSFRF);
  #--merge BSFRF and NMFS RHI dataframes
  qry<-"select
          n.HAULJOIN,n.BOTTOM_DEPTH, n.GEAR_TEMPERATURE,n.phi,n.sorting,
          n.SEX,n.SIZE,
          n.numIndivs as numNMFS,
          n.AREA_SWEPT_VARIABLE as aswNMFS,
          n.SAMPLING_FACTOR as sfNMFS,
          n.expFACTOR as expfNMFS,
          b.numIndivs as numBSFRF,
          b.AREA_SWEPT_VARIABLE as aswBSFRF,
          b.SAMPLING_FACTOR as sfBSFRF,
          b.expFACTOR as expfBSFRF,
          n.numIndivs+b.numIndivs as numTot,
          n.numIndivs/(n.numIndivs+b.numIndivs) as propNMFS,
          n.AREA_SWEPT_VARIABLE/b.AREA_SWEPT_VARIABLE as ratioAS
        from dfrRHI_NMFS as n, dfrRHI_BSFRF as b
        where
          n.HAULJOIN=b.HAULJOIN and
          n.SEX=b.SEX and
          n.SIZE=b.SIZE;";
  dfrProps<-sqldf::sqldf(qry) |>
              dplyr::filter(!is.na(propNMFS)) |> 
              dplyr::mutate(q=expfBSFRF/expfNMFS,  #--correct sense of ratio for selfisher: logit(propNMFS) = log[r(z)]+log(q)
                            obsLnR=log(propNMFS/(1-propNMFS))-log(q),
                            obsR=exp(obsLnR),
                            YEAR=as.character(iY));

  lstPropsAllYears[[as.character(iY)]] = dfrProps;
}#--iY
dfrPropsAllYears = dplyr::bind_rows(lstPropsAllYears) |> 
                     dplyr::mutate(HAULJOIN=paste0(YEAR,"-",HAULJOIN),
                                   YEAR=as.character(YEAR));
##--save dfrPropsAllYears----
wtsUtilities::saveObj(dfrPropsAllYears,file.path(dirThs,"rda_Step1b_dfrPropsAllSBSYears.RData"));
  
#--clean up a bit
rm(lstPropsAllYears,dfrProps,qry,
   dfrRHI_BSFRF,dfrRHI_NMFS,dfrREI_BSFRF,dfrREI_NMFS,dfrRCI_BSFRF,dfrRCI_NMFS,
   dfrRSD_SBS,dfrRHD_BSFRF_SBS,dfrRID_BSFRF_SBS,dfrRHD_NMFS_SBS,dfrRID_NMFS_SBS,
   dfrYSD_SBS,dfrYHD_BSFRF_SBS,dfrYID_BSFRF_SBS,dfrYHD_NMFS_SBS,dfrYID_NMFS_SBS);

#--calculate haul-level mean annual NMFS proportion of SBS catch----
#----and survey-level mean annual proportion based on CPUE (swept area-scaled values)
qry<-"select
        SEX,YEAR,SIZE,
        sum(numTot) as numTot,
        avg(propNMFS) as mnPropNMFS,
        avg(numNMFS*expfNMFS) as mnCPUE_NMFS,
        avg(numBSFRF*expfBSFRF) as mnCPUE_BSFRF
      from dfrPropsAllYears
      group by SEX,YEAR,SIZE;";
dfrMnPropsAllYears<-sqldf::sqldf(qry) |> 
                      dplyr::mutate(mnLnR       = log(mnPropNMFS/(1-mnPropNMFS)),
                                    mnCPUE_Prop = mnCPUE_NMFS/(mnCPUE_NMFS+mnCPUE_BSFRF),
                                    mnCPUE_lnR  = log(mnCPUE_Prop/(1-mnCPUE_Prop)),
                                    );
##--save dfrMnPropsAllYears----
wtsUtilities::saveObj(dfrMnPropsAllYears,file.path(dirThs,"rda_Step1c_dfrMnPropsAllSBSYears.RData"));
