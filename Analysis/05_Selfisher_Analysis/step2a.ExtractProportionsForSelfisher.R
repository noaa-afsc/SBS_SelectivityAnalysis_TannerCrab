#--extract size-specific empirical proportions from SBS data
#--for use in estimating selectivity with selfisher

require(magrittr);
#require(selfisher); require(bbmle);
source("extractHaulsAndIndivs.R");
#source("CalcFitPredictAndPlot.R")
load("data.SBS_Data.RData")

#--determine unique years
uYs<-unique(dfrSD_SBS$YEAR);
uYs<-c(2016,2017);#--Tanner-specific, west of Bristol Bay

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

#--extract data for selfisher analysis
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
    dfrRHD_BSFRF_SBS$HAULJOIN<-dfrRHD_BSFRF_SBS$newHAULJOIN; dfrRHD_BSFRF_SBS<-wtsUtilities::deleteCol(dfrRHD_BSFRF_SBS,"newHAULJOIN");
    dfrRID_BSFRF_SBS$HAULJOIN<-dfrRID_BSFRF_SBS$newHAULJOIN; dfrRID_BSFRF_SBS<-wtsUtilities::deleteCol(dfrRID_BSFRF_SBS,"newHAULJOIN");
    dfrRHD_NMFS_SBS$HAULJOIN <-dfrRHD_NMFS_SBS$newHAULJOIN;  dfrRHD_NMFS_SBS <-wtsUtilities::deleteCol(dfrRHD_NMFS_SBS, "newHAULJOIN");
    dfrRID_NMFS_SBS$HAULJOIN <-dfrRID_NMFS_SBS$newHAULJOIN;  dfrRID_NMFS_SBS <-wtsUtilities::deleteCol(dfrRID_NMFS_SBS, "newHAULJOIN");
    #--rename GIS_STATIONS to avoid duplicates
    dfrRSD_SBS$GIS_STATION      <-as.character(1:nS_SBS);
    dfrRHD_BSFRF_SBS$GIS_STATION<-as.character(1:nS_SBS);
    dfrRHD_NMFS_SBS$GIS_STATION <-as.character(1:nS_SBS);

    #--create single stratum for SBS stations
    dfrRSD_SBS$STRATUM<-"1";
    dfrRSD_SBS$STRATUM_CODE<-"1";
    dfrRSD_SBS$STRATUM_AREA          <-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
    dfrRSD_SBS$STRATUM_AREA_BYSTATION<-sum(dfrRSD_SBS$STATION_AREA,na.rm=TRUE);
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
    #--merge haul data with inidiv data
    mergeHaulData<-function(hdfr,idfr){
      qry<-"select
              h.HAULJOIN, BOTTOM_DEPTH, GEAR_TEMPERATURE, AREA_SWEPT_VARIABLE,
              SEX, SIZE, SAMPLING_FACTOR,
              AREA_SWEPT_VARIABLE/SAMPLING_FACTOR as expFACTOR,
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
            n.HAULJOIN,n.BOTTOM_DEPTH, n.GEAR_TEMPERATURE,
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
    dfrProps<-dfrProps[!is.na(dfrProps$propNMFS),];
    dfrProps$q<-dfrProps$expfBSFRF/dfrProps$expfNMFS;#--correct sense of ratio for selfisher logit(phi) = log[r(z)]+log(q)
    dfrProps$logit_resp<-log(dfrProps$propNMFS/(1-dfrProps$propNMFS))+log(dfrProps$q);
    dfrProps$YEAR<-as.character(iY);

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
  dfrPropsAllYears$HAULJOIN<-paste0(dfrPropsAllYears$YEAR,"-",dfrPropsAllYears$HAULJOIN);
  dfrPropsAllYears$YEAR<-as.character(dfrPropsAllYears$YEAR);
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
          sum(numNMFS) as numNMFS,
          sum(numTot) as numTot
        from dfrPropsAllYears
        group by SEX,YEAR,SIZE;";
  dfrMnPropsAllYears<-sqldf::sqldf(qry);
  dfrMnPropsAllYears$propNMFS   <-dfrMnPropsAllYears$expNumNMFS/dfrMnPropsAllYears$expNumTot;
  dfrMnPropsAllYears$logit_resp <-log(dfrMnPropsAllYears$propNMFS/(1-dfrMnPropsAllYears$propNMFS));
  dfrMnPropsAllYears$propNMFS0  <-dfrMnPropsAllYears$numNMFS/dfrMnPropsAllYears$numTot;
  dfrMnPropsAllYears$logit_resp0<-log(dfrMnPropsAllYears$propNMFS0/(1-dfrMnPropsAllYears$propNMFS0));
  wtsUtilities::saveObj(dfrMnPropsAllYears,"dfrMnPropsAllYears.RData");

  #plot proportions
  require(ggplot2);
  for (x in uXs){
    pdf(file=paste0("plots2a1.ObservedProportions_",x,".pdf"),width=9,height=6.5);
    tmp0<-dfrPropsAllYears[dfrPropsAllYears$SEX==x,];
    tmp1<-dfrMnPropsAllYears[dfrMnPropsAllYears$SEX==x,];
    p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=propNMFS,size=numTot,colour=YEAR));
    p <- p + geom_point(alpha=0.2,stat="identity",position=position_jitter(width=2.0)) + scale_size_area();
    p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=propNMFS,colour=YEAR),alpha=0.4,size=2);
    p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=propNMFS0,colour=YEAR),alpha=0.4,size=1);
    p <- p + facet_grid(rows=YEAR~.);
    p <- p + labs(x="size (mm CW)",y="proportion (NMFS/[NMFS+BSFRF])",colour="year",size="total\nnumber");
    print(p);
    p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=logit_resp,size=numTot,colour=YEAR));
    p <- p + geom_point(alpha=0.2,stat="identity",position=position_jitter(width=2.0)) + scale_size_area();
    p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=logit_resp,colour=YEAR),alpha=0.4,size=2);
    p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=logit_resp0,colour=YEAR),alpha=0.4,size=1);
    p <- p + facet_grid(rows=YEAR~.);
    p <- p + labs(x="size (mm CW)",y="logit(NMFS/[NMFS+BSFRF])+log(q)",colour="year",size="total\nnumber");
    print(p);
    dev.off();
  }
  pdf(file=paste0("plots2a2.LnQ.pdf"),width=9,height=4.5);
  p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=log(q),size=numTot,colour=YEAR));
  p <- p + geom_point(alpha=0.4,stat="identity",position=position_jitter(width=0.5)) + scale_size_area();
  p <- p + labs(x="size (mm CW)",y="log(q)",colour="sex",size="total/nnumber",shape="year");
  print(p);
  dev.off();
  rm(p,tmp0,tmp1,x);



