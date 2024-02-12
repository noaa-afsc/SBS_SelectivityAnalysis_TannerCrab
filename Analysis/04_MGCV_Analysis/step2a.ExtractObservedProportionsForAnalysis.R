#--extract size-specific empirical proportions from SBS + sed data
#--for use in estimating selectivity

#--load required packages
require(magrittr);
#--other required packages
#----dplyr
#----sqldf
#----wtsUtilities
#--source required functions
source("r_extractHaulsAndIndivs.R");

#--load SBS data with added sediment data
load("data.SBS_DataWithSedData.RData")

#--determine unique years
uYs<-unique(dfrSD_SBS$YEAR);
#uYs<-c(2016,2017);#--Tanner-specific, west of Bristol Bay

#--set up processing parameters
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
center<-as.numeric(c((maxSize["MALE"]+minSize["MALE"])/2,(maxSize["FEMALE"]+minSize["FEMALE"])/2)); names(center)<-c("MALE","FEMALE");
width <-as.numeric(c((maxSize["MALE"]-minSize["MALE"]),  (maxSize["FEMALE"]-minSize["FEMALE"])));   names(width) <-c("MALE","FEMALE");

#--drop "missing" sex data
nmiss_BSFRF<-dfrID_BSFRF_SBS %>% subset(SEX=="MISSING") %>% nrow();
nmiss_NMFS <-dfrID_NMFS_SBS  %>% subset(SEX=="MISSING") %>% nrow();
dfrID_BSFRF_SBS<-dfrID_BSFRF_SBS %>% subset(SEX!="MISSING");
dfrID_NMFS_SBS <-dfrID_NMFS_SBS  %>% subset(SEX!="MISSING");

#--apply cut points to sizes
dfrID_BSFRF_SBS$SIZE<-cutpts[cut(dfrID_BSFRF_SBS$SIZE,breaks=cutpts,include.lowest=TRUE)]+delta/2;
dfrID_NMFS_SBS$SIZE <-cutpts[cut(dfrID_NMFS_SBS$SIZE, breaks=cutpts,include.lowest=TRUE)]+delta/2;

#--extract data for analysis
dfrZCs<-NULL;
dfrESs<-NULL;
#nBs<-0;
#for (iB in 0){
#  cat("Processing",iB,"of",nBs,"\n");
  dfrPropsAllYears<-NULL;
  for (iY in uYs){
    cat("Processing",iY,"\n");
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
    # dfrRHD_BSFRF_SBS$HAULJOIN<-dfrRHD_BSFRF_SBS$newHAULJOIN; dfrRHD_BSFRF_SBS %<>% wtsUtilities::deleteCol("newHAULJOIN");
    # dfrRID_BSFRF_SBS$HAULJOIN<-dfrRID_BSFRF_SBS$newHAULJOIN; dfrRID_BSFRF_SBS %<>% wtsUtilities::deleteCol("newHAULJOIN");
    # dfrRHD_NMFS_SBS$HAULJOIN <-dfrRHD_NMFS_SBS$newHAULJOIN;  dfrRHD_NMFS_SBS  %<>% wtsUtilities::deleteCol("newHAULJOIN");
    # dfrRID_NMFS_SBS$HAULJOIN <-dfrRID_NMFS_SBS$newHAULJOIN;  dfrRID_NMFS_SBS  %<>% wtsUtilities::deleteCol("newHAULJOIN");
    dfrRHD_BSFRF_SBS %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) %>% dplyr::select(-newHAULJOIN);
    dfrRID_BSFRF_SBS %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) %>% dplyr::select(-newHAULJOIN);
    dfrRHD_NMFS_SBS  %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) %>% dplyr::select(-newHAULJOIN);
    dfrRID_NMFS_SBS  %<>% dplyr::mutate(HAULJOIN=newHAULJOIN) %>% dplyr::select(-newHAULJOIN);
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
    uXs<-unique(unique(dfrRCI_BSFRF$SEX),unique(dfrRCI_NMFS$SEX));
    uZs<-unique(unique(dfrRCI_BSFRF$SIZE),unique(dfrRCI_NMFS$SIZE));
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
            n.expFACTOR as expfNMFS,
            b.numIndivs as numBSFRF,
            b.expFACTOR as expfBSFRF,
            n.numIndivs+b.numIndivs as numTot,
            n.numIndivs/(n.numIndivs+b.numIndivs) as propNMFS
          from dfrRHI_NMFS as n, dfrRHI_BSFRF as b
          where
            n.HAULJOIN=b.HAULJOIN and
            n.SEX=b.SEX and
            n.SIZE=b.SIZE;";
    dfrProps<-sqldf::sqldf(qry);
    dfrProps %<>% subset(!is.na(propNMFS)) %>% 
                  dplyr::mutate(q=expfBSFRF/expfNMFS,  #--correct sense of ratio for selfisher: logit(propNMFS) = log[r(z)]+log(q)
                                lnR=log(propNMFS/(1-propNMFS))-log(q),
                                YEAR=as.character(iY));

    dfrPropsAllYears<-rbind(dfrPropsAllYears,dfrProps);
    # #--calculate empirical selectivity
    # lst<-calcEmpiricalSelectivity(dfrRSD_SBS,
    #                               dfrRHD_BSFRF_SBS,
    #                               dfrRID_BSFRF_SBS,
    #                               dfrRHD_NMFS_SBS,
    #                               dfrRID_NMFS_SBS,
    #                               aggBySex=aggBySex,
    #                               aggByMaturity=aggByMaturity,
    #                               aggByShellCondition=aggByShellCondition,
    #                               cutpts=cutpts,
    #                               truncate.low=truncate.low,
    #                               truncate.high=truncate.high,
    #                               showPlot=FALSE,
    #                               verbosity=0);
    # lst$dfrZCs[["iB"]]<-iB;
    # lst$dfrESs[["iB"]]<-iB;
    # dfrZCs<-rbind(dfrZCs,lst$dfrZCs);
    # dfrESs<-rbind(dfrESs,lst$dfrESs);
    # rm(lst);
  }#--iY
  dfrPropsAllYears %<>% dplyr::mutate(HAULJOIN=paste0(YEAR,"-",HAULJOIN),
                                      YEAR=as.character(YEAR));
  wtsUtilities::saveObj(dfrPropsAllYears,"dfrPropsAllYears.RData");
  
  #--clean up a bit
  rm(dfrProps,qry,
     dfrRHI_BSFRF,dfrRHI_NMFS,dfrREI_BSFRF,dfrREI_NMFS,dfrRCI_BSFRF,dfrRCI_NMFS,
     dfrRSD_SBS,dfrRHD_BSFRF_SBS,dfrRID_BSFRF_SBS,dfrRHD_NMFS_SBS,dfrRID_NMFS_SBS,
     dfrYSD_SBS,dfrYHD_BSFRF_SBS,dfrYID_BSFRF_SBS,dfrYHD_NMFS_SBS,dfrYID_NMFS_SBS);

  #--calculate annual NMFS proportion of SBS catch based on expanded values
  qry<-"select
          SEX,YEAR,SIZE,
          sum(numNMFS/expfNMFS) as expNumNMFS,
          sum(numNMFS/expfNMFS+numBSFRF/expfBSFRF) as expNumTot,
          sum(numTot) as numTot,
          avg(propNMFS) as mnPropNMFS
        from dfrPropsAllYears
        group by SEX,YEAR,SIZE;";
  dfrMnPropsAllYears<-sqldf::sqldf(qry);
  dfrMnPropsAllYears %<>% dplyr::mutate(propNMFS = expNumNMFS/expNumTot,
                                        lnR      = log(propNMFS/(1-propNMFS)));
  wtsUtilities::saveObj(dfrMnPropsAllYears,"dfrMnPropsAllYears.RData");




