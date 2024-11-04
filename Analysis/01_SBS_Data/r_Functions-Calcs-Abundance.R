#--various calculation functions

#--function to calculate CPUE by SBS haul/sex/size class
##--lst is the list object from Step1
calcCPUEs<-function(lst,
                    bySex=TRUE,
                    byMaturity=FALSE,
                    byShellCondition=FALSE,
                    bySize=TRUE,
                    cutpts=seq(from=0,to=185,by=5),
                    truncate.low=TRUE,
                    truncate.high=FALSE,
                    verbosity=0){
  require(tcsamSurveyData);
  #----calculate CPUE by size class for BSFRF SBS data
  dfrCPUE_BSFRF_SBS<-tcsamSurveyData::calcCPUE.ByHaul(
                       lst$dfrHD_BSFRF_SBS,
                       lst$dfrID_BSFRF_SBS,
                       bySex=bySex,
                       byMaturity=byMaturity,
                       byShellCondition=byShellCondition,
                       bySize=bySize,
                       cutpts=cutpts,
                       truncate.low=truncate.low,
                       truncate.high=truncate.high,
                       verbosity=verbosity) |> 
                      dplyr::filter(!(SEX %in% c("MISSING","HERMAPHRODITIC"))) |> 
                      dplyr::mutate(survey="BSFRF");

  #----calculate CPUE by size class for NMFS SBS data
  dfrCPUE_NMFS_SBS <-tcsamSurveyData::calcCPUE.ByHaul(
                       lst$dfrHD_NMFS_SBS,
                       lst$dfrID_NMFS_SBS,
                       bySex=bySex,
                       byMaturity=byMaturity,
                       byShellCondition=byShellCondition,
                       bySize=bySize,
                       cutpts=cutpts,
                       truncate.low=truncate.low,
                       truncate.high=truncate.high,
                       verbosity=verbosity) |> 
                      dplyr::filter(!(SEX %in% c("MISSING","HERMAPHRODITIC"))) |> 
                      dplyr::mutate(survey="NMFS");

  #----combine size compositions
  dfrCPUE<-rbind(dfrCPUE_BSFRF_SBS,dfrCPUE_NMFS_SBS);
  names(dfrCPUE)<-tolower(names(dfrCPUE));
  dfrCPUE<-dfrCPUE[,c("survey","year","gis_station","sampling_factor","area_swept_variable",
                      "sex","maturity","shell_condition","size",
                      "numindivs","numcpue")];
  names(dfrCPUE)<-c("fleet","y","gis_station","sampling_factor","area_swept_variable","x","m","s","z","n","val");
  dfrCPUE$x<-tolower(dfrCPUE$x);
  dfrCPUE$m<-tolower(dfrCPUE$m);
  dfrCPUE$s<-tolower(dfrCPUE$s);
  dfrCPUE$type<-"observed";

  return(dfrCPUE);
}


#--function to calculate expanded abundance size comps by SBS haul/sex/size class
##--lst is the list object from Step1
calcTotalSizeComps<-function(lst,
                              aggBySex=FALSE,
                              aggByMaturity=TRUE,
                              aggByShellCondition=TRUE,
                              cutpts=seq(from=25,to=185,by=5),
                              truncate.low=TRUE,
                              truncate.high=FALSE,
                              verbosity=0){
  require(tcsamSurveyData);
  #--regard annual set of paired hauls as comprising a single stratum
  dfrSD_SBS = lst$dfrSD_SBS |> dplyr::mutate(STRATUM="SBS");
  #----calculate ZCs for BSFRF SBS data
  lstZCs_BSFRF_SBS<-tcsamSurveyData::doCalcs_ZCs(
                      dfrSD_SBS,
                      lst$dfrHD_BSFRF_SBS |> dplyr::mutate(STRATUM="SBS"),
                      lst$dfrID_BSFRF_SBS,
                      useStratumArea=FALSE,
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
  lstZCs_NMFS_SBS<-tcsamSurveyData::doCalcs_ZCs(
                      dfrSD_SBS,
                      lst$dfrHD_NMFS_SBS |> dplyr::mutate(STRATUM="SBS"),
                      lst$dfrID_NMFS_SBS,
                      useStratumArea=FALSE,
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
  dfrZCs<-dfrZCs[,c("survey","year",
                    "sex","maturity","shell_condition","size",
                    "numindivs","totabundance")];
  names(dfrZCs)<-c("fleet","y","x","m","s","z","n","val");
  dfrZCs$x<-tolower(dfrZCs$x);
  dfrZCs$m<-tolower(dfrZCs$m);
  dfrZCs$s<-tolower(dfrZCs$s);
  dfrZCs$type<-"observed";

  return(dfrZCs);
}
