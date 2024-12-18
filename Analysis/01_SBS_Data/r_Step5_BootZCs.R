#--calculate survey size comps

doStep5<-function(){
  require(ggplot2);
  require(rlang);
  require(tables);
  require(wtsSizeComps);
  source("r_Functions-Calcs-Abundance.R")
  source("r_Functions-Figures-Abundance.R")
  
  #--create output list
  out = list();
  
  #--get haul-level data----
  lst1 = wtsUtilities::getObj("rda_Step1_SBS_RawData.RData");

  #--define size bins----
  cutpts = seq(from=24.5,to=184.5,by=5);
  zBs    = wtsUtilities::calcMidpoints(cutpts);
  out = c(out,list(cutpts=cutpts,zBs=zBs));
  
  #--do resampling, if necessary----
  doBootstraps = FALSE;
  if (doBootstraps){
    doBootstraps<-function(){
      nB = 1000; #--number of bootstrap repetitions/year
      lstZCs = list();#--temporary list to accumulate bootsrapped ZCs
      #--calculate "observed" size comps (bootrep=0)
      bootrep=0;
      dfrZCs_BSFRF_SBS = tcsamSurveyData::calcSizeComps.ByStratum(
                           lst1$dfrSD_SBS |> dplyr::mutate(STRATUM="SBS"),
                           tbl_hauls=lst1$dfrHD_BSFRF_SBS |> dplyr::mutate(STRATUM="SBS"),
                           tbl_indivs=lst1$dfrID_BSFRF_SBS |> dplyr::mutate(STRATUM="SBS"),
                           useStratumArea=FALSE, #--use sum of STATION_AREAs
                           bySex=TRUE,
                           byMaturity=FALSE,
                           byShellCondition=FALSE,
                           cutpts=cutpts,
                           truncate.low=TRUE,
                           truncate.high=FALSE,
                           verbosity=0) |> 
                         wtsUtilities::dropLevels(dropLevels=list(SEX=c('MISSING',"HERMAPHRODITIC"))) |> 
                         dplyr::mutate(survey="BSFRF",
                                       bootrep=bootrep);
      dfrZCs_NMFS_SBS  = tcsamSurveyData::calcSizeComps.ByStratum(
                           lst1$dfrSD_SBS |> dplyr::mutate(STRATUM="SBS"),
                           tbl_hauls=lst1$dfrHD_NMFS_SBS |> dplyr::mutate(STRATUM="SBS"),
                           tbl_indivs=lst1$dfrID_NMFS_SBS |> dplyr::mutate(STRATUM="SBS"),
                           useStratumArea=FALSE, #--use sum of STATION_AREAs
                           bySex=TRUE,
                           byMaturity=FALSE,
                           byShellCondition=FALSE,
                           cutpts=cutpts,
                           truncate.low=TRUE,
                           truncate.high=FALSE,
                           verbosity=0) |> 
                         wtsUtilities::dropLevels(dropLevels=list(SEX=c('MISSING',"HERMAPHRODITIC"))) |> 
                         dplyr::mutate(survey="NMFS",
                                       bootrep=bootrep);
      lstZCs[[paste("ALL","+",bootrep)]] = dplyr::bind_rows(dfrZCs_BSFRF_SBS,dfrZCs_NMFS_SBS);
  
      #--do bootstrapping for size comps
      uYs = unique(lst1$dfrSD_SBS$YEAR);
      rng_seed = 111111;
      set.seed(rng_seed);
      colsSD  = names(lst1$dfrSD_SBS);
      colsHDB = names(lst1$dfrHD_BSFRF_SBS);
      colsHDN = names(lst1$dfrHD_NMFS_SBS);
      colsIDB = names(lst1$dfrID_BSFRF_SBS);
      colsIDN = names(lst1$dfrID_NMFS_SBS);
      for (uY in uYs) {
        #--testing: uY = 2018;
        #--select data for YEAR, collapse to single STRATUM 
        ##--need to set useStratumArea=FALSE to calc area from individual station areas where appropriate
        dfrSD  = lst1$dfrSD_SBS |> dplyr::filter(YEAR==uY)       |> dplyr::mutate(STRATUM="SBS");
        dfrHDB = lst1$dfrHD_BSFRF_SBS |> dplyr::filter(YEAR==uY) |> dplyr::mutate(STRATUM="SBS");
        dfrHDN = lst1$dfrHD_NMFS_SBS  |> dplyr::filter(YEAR==uY) |> dplyr::mutate(STRATUM="SBS");
        dfrIDB = lst1$dfrID_BSFRF_SBS |> dplyr::inner_join(dfrHDB,by="HAULJOIN");
        dfrIDN = lst1$dfrID_NMFS_SBS  |> dplyr::inner_join(dfrHDN,by="HAULJOIN");
        nH     = nrow(dfrSD);#--number of hauls
        for (bootrep in 1:nB){
          #--testing: bootrep = 1;
          cat("Processing",uY,"bootrep",bootrep,"of",nB,"\n");
          #--resample stations with replacement
          dfrSDp  = dfrSD  |> dplyr::slice_sample(prop=1.00,replace=TRUE) |> 
                              dplyr::mutate(new_stn=dplyr::row_number());#--create new station index
          #--extract corresponding hauls
          dfrHDBp = dfrSDp |> dplyr::inner_join(dfrHDB,by=c("YEAR","STRATUM","GIS_STATION"));
          dfrHDNp = dfrSDp |> dplyr::inner_join(dfrHDN,by=c("YEAR","STRATUM","GIS_STATION"));
          #--extract corresponding individual data and resample individuals within (new) hauls and pop categories (just sex, here)
          dfrIDBp = dfrHDBp |> dplyr::left_join(dfrIDB,by="HAULJOIN",relationship="many-to-many") |> 
                               dplyr::group_by(new_stn,HAULJOIN,SEX,SIZE) |> 
                               dplyr::slice_sample(prop=1.00,replace=TRUE) |> 
                               dplyr::ungroup() |> 
                               dplyr::mutate(HAULJOIN=new_stn) |>          #--use new_stn as HAULJOIN index
                               dplyr::select(tidyselect::all_of(colsIDB)); #--keep original columns
          dfrIDNp = dfrHDNp |> dplyr::left_join(dfrIDN,by="HAULJOIN",relationship="many-to-many") |> 
                               dplyr::group_by(new_stn,HAULJOIN,SEX,SIZE) |> 
                               dplyr::slice_sample(prop=1.00,replace=TRUE) |> 
                               dplyr::ungroup() |> 
                               dplyr::mutate(HAULJOIN=new_stn) |>          #--use new_stn as HAULJOIN index
                               dplyr::select(tidyselect::all_of(colsIDN)); #--keep original columns
          #--re-index GIS_STATION and HAULJOIN to use new_stn in dfrSDp and dfrHD*p's
          dfrSDp = dfrSDp |> dplyr::mutate(GIS_STATION=new_stn) |> 
                             dplyr::select(tidyselect::all_of(colsSD));
          dfrHDBp = dfrHDBp |> dplyr::mutate(GIS_STATION=new_stn,
                                             HAULJOIN=new_stn) |> 
                               dplyr::select(tidyselect::all_of(colsHDB));
          dfrHDNp = dfrHDNp |> dplyr::mutate(GIS_STATION=new_stn,
                                             HAULJOIN=new_stn) |> 
                               dplyr::select(tidyselect::all_of(colsHDN));
          #--calculate size comps
          dfrZCs_BSFRF_SBS = tcsamSurveyData::calcSizeComps.ByStratum(
                                 dfrSDp,
                                 tbl_hauls=dfrHDBp,
                                 tbl_indivs=dfrIDBp,
                                 useStratumArea=FALSE, #--use sum of STATION_AREAs
                                 bySex=TRUE,
                                 byMaturity=FALSE,
                                 byShellCondition=FALSE,
                                 cutpts=cutpts,
                                 truncate.low=TRUE,
                                 truncate.high=FALSE,
                                 verbosity=0) |> 
                               wtsUtilities::dropLevels(dropLevels=list(SEX=c('MISSING',"HERMAPHRODITIC"))) |> 
                               dplyr::mutate(survey="BSFRF",
                                             bootrep=bootrep);
          dfrZCs_NMFS_SBS  = tcsamSurveyData::calcSizeComps.ByStratum(
                                 dfrSDp,
                                 tbl_hauls=dfrHDNp,
                                 tbl_indivs=dfrIDNp,
                                 useStratumArea=FALSE, #--use sum of STATION_AREAs
                                 bySex=TRUE,
                                 byMaturity=FALSE,
                                 byShellCondition=FALSE,
                                 cutpts=cutpts,
                                 truncate.low=TRUE,
                                 truncate.high=FALSE,
                                 verbosity=0) |> 
                               wtsUtilities::dropLevels(dropLevels=list(SEX=c('MISSING',"HERMAPHRODITIC"))) |> 
                               dplyr::mutate(survey="NMFS",
                                             bootrep=bootrep);
          lstZCs[[paste(uY,"+",bootrep)]] = dplyr::bind_rows(dfrZCs_BSFRF_SBS,dfrZCs_NMFS_SBS);
        }#--i
      }#--uY
      dfrBootZCs = dplyr::bind_rows(lstZCs) |> 
                     dplyr::rename_with(tolower) |> 
                     dplyr::rename(tidyselect::all_of(c(gear="survey",
                                                        y="year",
                                                        x="sex",
                                                        m="maturity",
                                                        s="shell_condition",
                                                        z="size",
                                                        area="stratum_area",    #--in sq. nmi
                                                        abd="totabundance",     #--in millions of crab
                                                        bio="totbiomass"))) |>  #--in 1,000's mt
                     dplyr::mutate(x=tolower(x),
                                   m=tolower(m),
                                   s=tolower(s),
                                   z=z+3.0,      #--shift to middle of size bin
                                   abd=1000*abd);#--now in 1,000's of crab
      wtsUtilities::saveObj(dfrBootZCs,"rda_Step5_BootZCs.RData");
      return(dfrBootZCs)
    }
    dfrBootZCs = doBootstraps();
  } else {
    dfrBootZCs = wtsUtilities::getObj("rda_Step5_BootZCs.RData");
  }
  dfrBootZCs = dfrBootZCs |> dplyr::mutate(abd=abd/area); #--NOW in 1,000's/sq. nmi. [move into above!!]
  
  nB = max(dfrBootZCs$bootrep);
  out = c(out,list(dfrBootZCs=dfrBootZCs,nB=nB));

  #--do figures for boostrapped ZCs----
  #| label: fig-BootZCs-BSFRF
  cap = paste0('"Observed" (thick line, dots) and bootstrapped (thin lines) size compositions from the BSFRF gear. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times, with size compositions for both gear types calculated each time.");
  p = plotBZCs(dfr=dfrBootZCs,gear_types="BSFRF",sexs=c("female","male"),factor=x);
  out = c(out,list(`fig-BootZCs-BSFRF`=list(p=p,cap=cap)));
  
  #| label: fig-BootZCs-NMFS
  cap = paste0('"Observed" (thick line, dots) and bootstrapped (thin lines) size compositions from the NMFS gear. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times, with size compositions for both gear types calculated each time.");
  p = plotBZCs(dfr=dfrBootZCs,gear_types="NMFS",sexs=c("female","male"),factor=x);
  out = c(out,list(`fig-BootZCs-NMFS`=list(p=p,cap=cap)));
  
  #| label: fig-BootZCs-Ms
  cap = paste0('"Observed" (thick line, dots) and bootstrapped (thin lines) size compositions for males. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times, with size compositions for both gear types calculated each time.");
  p = plotBZCs(dfr=dfrBootZCs,gear_types=c("BSFRF","NMFS"),sexs=c("male"),factor=gear);
  out = c(out,list(`fig-BootZCs-Ms`=list(p=p,cap=cap)));
  
  #| label: fig-BootZCs-Fs
  cap = paste0('"Observed" (thick line, dots) and bootstrapped (thin lines) size compositions for females. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times, with size compositions for both gear types calculated each time.");
  p = plotBZCs(dfr=dfrBootZCs,gear_types=c("BSFRF","NMFS"),sexs=c("female"),factor=gear,xlims=c(25,130));
  out = c(out,list(`fig-BootZCs-Fs`=list(p=p,cap=cap)));
  
  #--do tables----
  #--tbl.SurveyLevelZCsMales}
  cap = "Survey-level size compositions (1,000s/sq. nmi.) for males from the SBS studies, by 5-mm size bin and gear type."
  tblr = tabular(Factor(z,name="size")~Factor(y)*Factor(gear)*abd*sum,
                 data=dfrBootZCs |> dplyr::filter(x=="male",bootrep==0));
  colLabels(tblr) = colLabels(tblr)[c(2,4),];
  kbl = tblr |> wtsQMD::convert_tblr_to_kbl(c(1,1+c(2,4,6,8,10,12)),
                                            isPDF=wtsQMD::isOutputPDF());
  #print(kbl);
  out = c(out,list(tblrSurveyLevelZCsMales=tblr));
  rm(cap,tblr,kbl);

  #--tbl.SurveyLevelZCsFemales}
  cap = "Survey-level size compositions (1,000s/sq. nmi.) for females from the SBS studies, by 5-mm size bin and gear type."
  tblr = tabular(Factor(z,name="size")~Factor(y)*Factor(gear)*abd*sum,,
                 data=dfrBootZCs |> dplyr::filter(x=="female",bootrep==0));
  colLabels(tblr) = colLabels(tblr)[c(2,4),];
  kbl = tblr |> wtsQMD::convert_tblr_to_kbl(c(1,1+c(2,4,6,8,10,12)),
                                            isPDF=wtsQMD::isOutputPDF());
  #print(kbl);
  out = c(out,list(tblrSurveyLevelZCsFemales=tblr));
  rm(cap,tblr,kbl);

  #--compute statistics for bootstrapped size comps----
  dfrStatsBZCs = dfrBootZCs |> dplyr::filter(bootrep>0) |> 
                   dplyr::group_by(gear,y,x,z) |> 
                   dplyr::summarize(mn=mean(abd),
                                    md=median(abd),
                                    vr=var(abd),
                                    sd=sqrt(vr),
                                    l90=quantile(abd,0.05),
                                    u90=quantile(abd,0.95)) |> 
                   dplyr::ungroup();

  #--do plots for bootstrap statistics----
  #| label: fig-StatsBootZCs-Ms
  cap = paste0('"Observed" (black line, dots) and bootstrapped median (lines) and 90 percent confidence intervals by gear type for male size compositions. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times by study year, with size compositions for both gear types calculated each time.");
  p = plotStatsBZCs(dfrBootZCs |> dplyr::filter(bootrep==0),
                    dfrStatsBZCs,gear_types=c("BSFRF","NMFS"),sexs=c("male"),factor=gear);
  out = c(out,list(`fig-StatsBootZCs-Ms`=list(p=p,cap=cap)));
  
  #| label: fig-StatsBootZCs-Fs
  cap = paste0('"Observed" (black line, dots) and bootstrapped median (lines) and 90 percent confidence intervals by gear type for female size compositions. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times by study year, with size compositions for both gear types calculated each time.");
  p = plotStatsBZCs(dfrBootZCs |> dplyr::filter(bootrep==0),
                    dfrStatsBZCs,gear_types=c("BSFRF","NMFS"),sexs=c("female"),factor=gear,xlims=c(25,135));
  out = c(out,list(`fig-StatsBootZCs-Fs`=list(p=p,cap=cap)));
  
  #| label: fig-CVsBootZCs-Ms
  cap = paste0('CVs from bootstrapped male size compositions. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times by study year, with size compositions for both gear types calculated each time.");
  p = plotCVsCPUE(dfrStatsBZCs |> dplyr::filter(x=="male"),mn,sd,factor_=gear,ylab="Size Composition CVs",xlims=c(25,185));
  out = c(out,list(`fig-CVsBootZCs-Ms`=list(p=p,cap=cap)));
  
  #| label: fig-CVsBootZCs-Fs
  cap = paste0('CVs from bootstrapped female size compositions. ',
               "Bootstrapping consisted of hierarchical resampling of paired hauls and measured crab within hauls ",
               nB," times by study year, with size compositions for both gear types calculated each time.");
  p = plotCVsCPUE(dfrStatsBZCs |> dplyr::filter(x=="female"),mn,sd,factor_=gear,ylab="Size Composition CVs",xlims=c(25,135));
  out = c(out,list(`fig-CVsBootZCs-Fs`=list(p=p,cap=cap)));
  
  wtsUtilities::saveObj(out,"rda_Step5_BootZCs_AllResults.RData");
  return(out);
}
out = doStep5();
#rm(out);