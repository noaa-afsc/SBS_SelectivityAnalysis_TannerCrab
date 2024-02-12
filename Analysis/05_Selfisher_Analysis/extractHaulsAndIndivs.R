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
  cat("Resampling from",nH,"hauls\n");
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
