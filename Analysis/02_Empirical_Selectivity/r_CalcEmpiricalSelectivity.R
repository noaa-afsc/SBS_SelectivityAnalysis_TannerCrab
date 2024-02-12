
calcEmpiricalSelectivity<-function(dfrSD_SBS,
                                    dfrHD_BSFRF_SBS,
                                    dfrID_BSFRF_SBS,
                                    dfrHD_NMFS_SBS,
                                    dfrID_NMFS_SBS,
                                    aggBySex=FALSE,
                                    aggByMaturity=TRUE,
                                    aggByShellCondition=TRUE,
                                    cutpts=seq(from=25,to=185,by=5),
                                    truncate.low=TRUE,
                                    truncate.high=FALSE,
                                    showPlot=FALSE,
                                    verbosity=0){
  #----calculate ZCs for BSFRF SBS data
  lstZCs_BSFRF_SBS<-doCalcs_ZCs(dfrSD_SBS,
                               dfrHD_BSFRF_SBS,
                               dfrID_BSFRF_SBS,
                               calcByEW166=FALSE,
                               aggBySex=aggBySex,
                               aggByMaturity=aggByMaturity,
                               aggByShellCondition=aggByShellCondition,
                               cutpts=cutpts,
                               truncate.low=truncate.low,
                               truncate.high=truncate.high,
                               dropLevels=list(SEX=c("MISSING","HERMAPHRODITIC")),
                               verbosity=verbosity);
  dfrZCs_BSFRF_SBS<-lstZCs_BSFRF_SBS$EBS;
  dfrZCs_BSFRF_SBS$survey<-"BSFRF";

  #----calculate ZCs for NMFS SBS data
  lstZCs_NMFS_SBS<-doCalcs_ZCs(dfrSD_SBS,
                               dfrHD_NMFS_SBS,
                               dfrID_NMFS_SBS,
                               calcByEW166=FALSE,
                               aggBySex=aggBySex,
                               aggByMaturity=aggByMaturity,
                               aggByShellCondition=aggByShellCondition,
                               cutpts=cutpts,
                               truncate.low=truncate.low,
                               truncate.high=truncate.high,
                               dropLevels=list(SEX=c("MISSING","HERMAPHRODITIC")),
                               verbosity=verbosity);
  dfrZCs_NMFS_SBS<-lstZCs_NMFS_SBS$EBS;
  dfrZCs_NMFS_SBS$survey<-"NMFS";

  #----combine size compositions
  dfrZCs<-rbind(dfrZCs_BSFRF_SBS,dfrZCs_NMFS_SBS);
  names(dfrZCs)<-tolower(names(dfrZCs));
  dfrZCs<-dfrZCs[,c("survey","year","sex","maturity","shell_condition","size","numindivs","totabundance")];
  names(dfrZCs)<-c("fleet","y","x","m","s","z","n","val");
  dfrZCs$x<-tolower(dfrZCs$x);
  dfrZCs$m<-tolower(dfrZCs$m);
  dfrZCs$s<-tolower(dfrZCs$s);
  dfrZCs$type<-"observed";

  #--Calculate empirical selectivity as ratio of NMFS to BSFRF size compositions
  dfrZCsp1<-reshape2::dcast(dfrZCs,"type+y+x+z~fleet",fun.aggregate=wtsUtilities::Sum,value.var="n");
  names(dfrZCsp1)[5:6]<-paste0("n_",names(dfrZCsp1)[5:6]);
  dfrZCsp2<-reshape2::dcast(dfrZCs,"type+y+x+z~fleet",fun.aggregate=wtsUtilities::Sum,value.var="val");
  dfrZCsp2[["emp_sel"]]<-dfrZCsp2$NMFS/dfrZCsp2$BSFRF;
  idx<-is.infinite(dfrZCsp2[["emp_sel"]]);
  dfrZCsp2[["emp_sel"]][idx]<-NaN;
  dfrZCsp<-cbind(dfrZCsp1,BSFRF=dfrZCsp2$BSFRF,NMFS=dfrZCsp2$NMFS,emp_sel=dfrZCsp2$emp_sel);

  return(list(dfrESs=dfrZCsp,dfrZCs=dfrZCs));
}
