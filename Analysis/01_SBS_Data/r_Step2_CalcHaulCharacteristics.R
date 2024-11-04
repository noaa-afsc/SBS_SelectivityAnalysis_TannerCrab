#--calculate haul characteristics
##--define function to create objects
calcHaulCharacteristics<-function(){
  require(ggplot2)
  require(tables);
  #--create list for output
  out = list();

  #--get haul-level areas swept, depth, and temperature data
  lst = wtsUtilities::getObj("rda_Step1_SBS_RawData.RData");
  
  #--extract areas swept, depth, and temperature----
  dfrADT = dplyr::bind_rows(
                lst$dfrHD_BSFRF_SBS |> 
                  dplyr::select(YEAR,station=GIS_STATION,`area swept`=AREA_SWEPT_VARIABLE,
                                depth=BOTTOM_DEPTH,temp=GEAR_TEMPERATURE) |> 
                  dplyr::rename_with(tolower,everything()) |> 
                  dplyr::mutate(gear="BSFRF"),
                lst$dfrHD_NMFS_SBS |> 
                  dplyr::select(YEAR,station=GIS_STATION,`area swept`=AREA_SWEPT_VARIABLE,
                                depth=BOTTOM_DEPTH,temp=GEAR_TEMPERATURE) |> 
                  dplyr::rename_with(tolower,everything()) |> 
                  dplyr::mutate(gear="NMFS"));
  out = c(out,list(dfrADT=dfrADT));
  
  #--create summary table of haul characteristics
  tblrADT = tabular(Factor(year) ~ (n=1)*Factor(gear) + 
                                   (`area swept`+
                                      Heading("bottom depth")*depth+
                                      Heading("bottom temp.")*temp)*Factor(gear)*(min+max),
                    data=dfrADT);
  colLabels(tblrADT) = colLabels(tblrADT)[c(1,3,4),];
  colLabels(tblrADT)[1:3,1] = c("","","N");
  tblrADT = tblrADT[,-2];
  out = c(out,list(tblrADT=tblrADT));
  
  #--print table (needs to be done again in qmd)
  kblADT = tblrADT |> wtsQMD::convert_tblr_to_kbl(c(1,2,4,6,8,10,12),
                                                  isPDF=wtsQMD::isOutputPDF());
  # print(kblADT);
  rm(kblADT);

  #--calc stats for areas swept----
  dfrAS = dfrADT |> 
          tidyr::pivot_wider(id_cols=c("year","station"),
                             names_from="gear",values_from=c("area swept")) |> 
          dplyr::mutate(ratio=`NMFS`/`BSFRF`);
  maxASR = max(dfrAS$ratio);
  minASR = max(dfrAS$ratio);
  mnAS = mean(dfrAS$ratio);
  out = c(out,list(dfrAS=dfrAS,maxASR=maxASR,minASR=minASR,mnAS=mnAS));

  #--calc stats for haul depths----
  dfrD = dfrADT |> 
          tidyr::pivot_wider(id_cols=c("year","station"),names_from="gear",values_from="depth");
  minDiffD = min(dfrD$NMFS-dfrD$BSFRF);
  maxDiffD = max(dfrD$NMFS-dfrD$BSFRF);
  mnDiffD = mean(dfrD$NMFS-dfrD$BSFRF);
  out = c(out,list(dfrD=dfrD,maxDiffD=maxDiffD,minDiffD=minDiffD,mnDiffD=mnDiffD));

  #--calc stats for bottom temperatures----
  dfrT = dfrADT |> 
          tidyr::pivot_wider(id_cols=c("year","station"),names_from="gear",values_from="temp") |> 
          dplyr::mutate(diff=NMFS-BSFRF);
  dfrMaxDiffT = dfrT |> dplyr::filter(abs(diff)==max(abs(diff),na.rm=TRUE));
  mxDiffT = dfrMaxDiffT$diff;
  minDiffTp = min((dfrT |> dplyr::filter(abs(diff)<abs(mxDiffT)))$diff);
  maxDiffTp = max((dfrT |> dplyr::filter(abs(diff)<abs(mxDiffT)))$diff);
  mnDiffTp  = mean((dfrT |> dplyr::filter(abs(diff)<abs(mxDiffT)))$diff);
  out = c(out,list(dfrT=dfrT,mxDiffT=mxDiffT,maxDiffTp=maxDiffTp,minDiffTp=minDiffTp,mnDiffTp=mnDiffTp));

#| label: fig-AreasSwept
  cap = "Comparison of areas swept by the NMFS (x-axis) and BSFRF (y-axis) gear types. Each histogram is normalized by its maximum to show relative frequencies across bins."
  xmx = max(dfrAS$`NMFS`,na.rm=TRUE);
  ymx = max(dfrAS$`BSFRF`,na.rm=TRUE);
  p = wtsPlots::ggMarginal_Hist2D(dfrAS,`NMFS`,`BSFRF`,
                                  xlab="NMFS area swept (sq. nmi.)",
                                  ylab="BSFRF area swept (sq. nmi.)",
                                  xparams=list(limits=c(0,xmx)),
                                  yparams=list(limits=c(0,ymx)),
                                  legend.position=c(0.01,0.99),
                                  legend.justification=c(0,1));
  out = c(out,list(figASs=list(p=p,cap=cap)));
  
#| label: fig-AreasSweptRatios
  cap = "Ratios of haul-specific area swept by the NMFS and BSFRFs gears."
  p = ggplot(dfrAS,aes(x=ratio)) + 
        geom_histogram() + 
        geom_vline(xintercept=mnAS,linetype=2) + 
        scale_x_continuous(breaks=seq(0,20,2)) + 
        labs(x="area swept ratios (NMFS/BSFRF)") + 
        wtsPlots::getStdTheme();
  out = c(out,list(figASRs=list(p=p,cap=cap)));

#| label: fig-HaulDepths
  cap = "Comparison of the difference (NMFS-BSFRF) paired-haul depths, by NMFS haul depth. Colour indicates the year each haul was conducted. The smooth line (and shading) indcates the overall trend with depth. The dotted line indicates the mean difference."
  p = ggplot(dfrD,aes(x=NMFS,y=NMFS-BSFRF,colour=as.character(year))) + 
        geom_smooth(colour="black") + 
        geom_point(shape=22) + geom_hline(yintercept=0,linetype=2) + 
        geom_hline(yintercept=mnDiffD,linetype=3) + 
        labs(x="NMFS haul depth (m)",y="NMFS-BSFRF haul depth (m)",colour="year") + 
        wtsPlots::getStdTheme() + theme(legend.position=c(0.01,0.99),legend.justification=c(0,1));
  out = c(out,list(figHDs=list(p=p,cap=cap)));

#| label: fig-BottomTemps
  cap = "Comparison of measured bottom temperatures by paired haul, with the largest difference excluded to show detail. The difference (NMFS-BSFRF) between the paired haul depths is plotted against the NMFS haul depth. Colour indicates the year each haul was conducted. The smooth line (and shading) indcates the overall trend with depth. The inset shows the same, but with the largest difference included."
  p1 = ggplot(dfrT,
             aes(x=NMFS,y=NMFS-BSFRF,colour=as.character(year))) + 
#        geom_smooth(colour="black") + 
        geom_point(shape=22) + geom_hline(yintercept=0,linetype=2) + 
        labs(x="NMFS bottom temperature",
             y="paired haul bottom temperature differences \n(NMFS-BSFRF) ",
             colour="year") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="none",legend.justification=c(0,0),
              axis.title=element_blank());
  tbl = tibble::tibble(x=0.99,y=0.99,plot=list(p1));
  p2 = ggplot(dfrT |> dplyr::filter(abs(NMFS-BSFRF)<2),
             aes(x=NMFS,y=NMFS-BSFRF,colour=as.character(year))) + 
        ggpp::geom_plot_npc(data=tbl,mapping=aes(npcx=x,npcy=y,label=plot)) + 
        geom_smooth(colour="black") + 
        geom_point(shape=22) + geom_hline(yintercept=0,linetype=2) + 
                labs(x="NMFS bottom temperature",
                     y="paired haul bottom temperature differences \n(NMFS-BSFRF) ",
                     colour="year") +  
        wtsPlots::getStdTheme() + theme(legend.position=c(0.01,0.99),legend.justification=c(0,1));
  out = c(out,list(figBTs=list(p=p2,cap=cap)));
  
  #--save results
  wtsUtilities::saveObj(out,"rda_Step2_SBS_ADTs.RData");
  return(out);
}
#--run function
out = calcHaulCharacteristics();
#--clean up
rm(out,calcHaulCharacteristics);

