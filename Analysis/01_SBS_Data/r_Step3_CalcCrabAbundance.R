#--calculate various aspects of crab abundance

calcCrabAbundance<-function(){
  #--calculate "raw" size comps for numbers caught by station and gear type
  #----DON'T use tcsamSurveyData functions for this (they scale by area swept)
  require(ggplot2);
  require(rlang);
  require(tables);
  require(wtsSizeComps);
  source("r_Functions-Figures-Abundance.R")
  source("r_Functions-Calcs-Abundance.R")
  
  #--create list for output----
  out = list();

  #--get haul-level data----
  lst = wtsUtilities::getObj("rda_Step1_SBS_RawData.RData");
  
  #--define size bins----
  cutpts = seq(-0.5,204.5,5);
  zBs    = wtsUtilities::calcMidpoints(cutpts);
  out = c(out,list(cutpts=cutpts,zBs=zBs));
  
  #--calculate numbers caught by haul----
  ##--combine individual data from gear types
  dfrID_SBS = dplyr::bind_rows(
                lst$dfrID_BSFRF_SBS |> 
                  dplyr::inner_join(lst$dfrHD_BSFRF_SBS,by="HAULJOIN") |> 
                  dplyr::select(YEAR,station=GIS_STATION,SEX,MATURITY,SHELL_CONDITION,SIZE,
                                numIndivs,SAMPLING_FACTOR) |> 
                  dplyr::mutate(gear="BSFRF"),
                lst$dfrID_NMFS_SBS |> 
                  dplyr::inner_join(lst$dfrHD_NMFS_SBS,by="HAULJOIN") |> 
                  dplyr::select(YEAR,station=GIS_STATION,SEX,MATURITY,SHELL_CONDITION,SIZE,
                                numIndivs,SAMPLING_FACTOR) |> 
                  dplyr::mutate(gear="NMFS")) |> 
               dplyr::rename_with(tolower,.cols=!numIndivs) |> 
               dplyr::filter(sex %in% c("MALE","FEMALE"));#--drop other categories
  #--determine unique year/station combinations and all size bins
  dfrYSZ = dfrID_SBS |> dplyr::distinct(year,station) |> 
             dplyr::cross_join(tibble::tibble(size=zBs));
  #--calc size comps based on numbers sampled----
  dfrZCs_Num = wtsSizeComps::calcSizeComps(dfrID_SBS,"size","numIndivs",
                                           c("gear","year","station","sex"),
                                           cutpts,expandToAllFactorCombos=FALSE) |> 
                 dplyr::mutate(size=size+2.5); #--adjust to center of size bins
  #----expand size comps to all year/station/size bin combinations
  dfrZCs_Num = dfrYSZ |> dplyr::left_join(dfrZCs_Num,by=c("year","station","size")) |> 
                 dplyr::mutate(sampled=ifelse(is.na(numIndivs),0,numIndivs)) |> 
                 dplyr::arrange(year,station,sex,size,gear) |> 
                 dplyr::select(!numIndivs);
  #--calc size comps based on estimated numbers caught----
  dfrZCs_Cat = wtsSizeComps::calcSizeComps(dfrID_SBS,"size","sampling_factor",
                                           c("gear","year","station","sex"),
                                           cutpts,expandToAllFactorCombos=FALSE) |> 
                 dplyr::mutate(size=size+2.5); #--adjust to center of size bins
  #----expand size comps to all year/station/size bin combinations
  dfrZCs_Cat = dfrYSZ |> dplyr::left_join(dfrZCs_Cat,by=c("year","station","size")) |> 
                 dplyr::mutate(caught=ifelse(is.na(sampling_factor),0,sampling_factor)) |> 
                 dplyr::arrange(year,station,sex,size,gear) |> 
                 dplyr::select(!sampling_factor);
  #--combine tables
  dfrZCs_RawByStn = dplyr::bind_cols(dfrZCs_Num,caught=dfrZCs_Cat$caught);
  rm(dfrZCs_Cat,dfrZCs_Num);
  
  dfrZCs_RawTot = dfrZCs_RawByStn |> 
                    dplyr::group_by(year,sex,size,gear) |> 
                    dplyr::summarize(ss=wtsUtilities::Sum(ss),
                                     sampled=wtsUtilities::Sum(sampled),
                                     caught=wtsUtilities::Sum(caught)) |>
                    dplyr::ungroup() |> 
                    tidyr::pivot_longer(c(sampled,caught),names_to="type",values_to="total") |>
                    dplyr::mutate(type=factor(type,levels=c("sampled","caught")));
  out = c(out,list(dfrZCs_RawTot=dfrZCs_RawTot));

#| label: fig-RawNumbersMales
  cap = "Total numbers of males sampled (measured) in the SBS studies, by 5-mm size bin and gear type (solid lines). Also shown are the estimated numbers caught prior to any sub-sampling, by 5-mm size bin and gear type (dotted lines). The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotTotRawZCs(dfrZCs_RawTot,"MALE",total,"male numbers");
  out = c(out,list(figRawNumbersMales=list(p=p,cap=cap)));
  
#| label: fig-RawNumbersFemales
  cap = "Total numbers of females sampled (measured) in the SBS studies, by 5-mm size bin and gear type (solid lines). Also shown are the estimated numbers caught prior to any sub-sampling, by 5-mm size bin and gear type (dotted lines). The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotTotRawZCs(dfrZCs_RawTot,"FEMALE",total,"female numbers");
  out = c(out,list(figRawNumbersFemales=list(p=p,cap=cap)));
  
#| label: fig-PctNon0sCPUEMales
  cap = "The percentage of non-zero catch hauls for male CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotPctNon0sCPUE(dplyr::filter(dfrStatsCPUE,x=='male'));
  out = c(out,list(figPctNon0sCPUEMales=list(p=p,cap=cap)));

#| label: fig-PctNon0sCPUEFemales
  cap = "The percentage of non-zero catch hauls for female CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotPctNon0sCPUE(dplyr::filter(dfrStatsCPUE,x=='female'));
  out = c(out,list(figPctNon0sCPUEFemales=list(p=p,cap=cap)));

#| label: fig-PctNon0sCPUEBSFRF
  cap = "The percentage of non-zero catch hauls for the BSFRF gear in the SBS studies, by 5-mm size bin and sex. The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotPctNon0sCPUE(dplyr::filter(dfrStatsCPUE,fleet=='BSFRF'),factor=x);
  out = c(out,list(figPctNon0sCPUEBSFRF=list(p=p,cap=cap)));

#| label: fig-PctNon0sCPUENMFS
  cap = "The percentage of non-zero catch hauls for the NMFS gear in the SBS studies, by 5-mm size bin and sex. The dotted vertical line marks the lower limit of the size bins used in the assessment."
  p = plotPctNon0sCPUE(dplyr::filter(dfrStatsCPUE,fleet=='NMFS'),factor=x);
  out = c(out,list(figPctNon0sCPUENMFS=list(p=p,cap=cap)));

#--calculate CPUE by haul----
  dfrCPUE = calcCPUEs(lst,
                      bySex=TRUE,
                      byMaturity=FALSE,
                      byShellCondition=FALSE,
                      bySize=TRUE,
                      cutpts=cutpts);
  avgASw = dplyr::distinct(dfrCPUE,fleet,gis_station,area_swept_variable) |> 
           dplyr::group_by(fleet) |> 
           dplyr::summarize(mean_as=mean(area_swept_variable,na.rm=TRUE)) |> 
           dplyr::ungroup();
  dfrStatsCPUE = dplyr::filter(dfrCPUE,!is.na(val)) |> 
                 dplyr::group_by(fleet,y,x,z) |> 
                 dplyr::summarize(n=dplyr::n(),
                                  nn0=sum(val>0),
                                  mn=mean(val),
                                  min=min(val),
                                  max=max(val),
                                  vr=var(val),
                                  sd=sqrt(vr),
                                  l90=quantile(val,0.05),
                                  u90=quantile(val,0.95)) |> 
                 dplyr::ungroup() |> 
                 dplyr::filter(mn>0);
  out = c(out,list(dfrCPUE=dfrCPUE,avgASw=avgASw,dfrStatsCPUE=dfrStatsCPUE));

#| label: fig-StatsCPUEMales
  cap = "Male CPUE in the SBS studies, by 5-mm size bin and gear type. Solid lines: mean; shaded area: empirical 90% CI. The dotted vertical line marks the lower limit of the size bins used in the assessment. The y-axis is on a log scale to facilitate comparison across the range of values."
  p = plotStatsCPUE(dplyr::filter(dfrStatsCPUE,x=='male'),mn,l90,u90);
  out = c(out,list(figStatsCPUEMales=list(p=p,cap=cap)));

#| label: fig-StatsCPUEFemales
  cap = "Female CPUE in the SBS studies, by 5-mm size bin and gear type. Solid lines: mean; shaded area: empirical 90% CI. The dotted vertical line marks the lower limit of the size bins used in the assessment. The y-axis is on a log scale to facilitate comparison across the range of values."
  p = plotStatsCPUE(dplyr::filter(dfrStatsCPUE,x=='female'),mn,l90,u90);
  out = c(out,list(figStatsCPUEFemales=list(p=p,cap=cap)));

#| label: fig-CVsCPUEMales
  cap = "CVs for male CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment. The horizontal line indicates a CV of 1."
  p = plotCVsCPUE(dplyr::filter(dfrStatsCPUE,x=='male'),mn,sd);
  out = c(out,list(figCVsCPUEMales=list(p=p,cap=cap)));

#| label: fig-CVsCPUEFemales
  cap = "CVs for female CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment. The horizontal line indicates a CV of 1."
  p = plotCVsCPUE(dplyr::filter(dfrStatsCPUE,x=='female'),mn,sd);
  out = c(out,list(figCVsCPUEFemales=list(p=p,cap=cap)));

  #--save objects
  wtsUtilities::saveObj(out,"rda_Step3_SBS_CrabAbundance.RData");
  return(out)

}
#--run function
out = calcCrabAbundance();
#--clean up
rm(out,calcCrabAbundance);

