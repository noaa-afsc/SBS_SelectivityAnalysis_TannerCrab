#--sediment characteristics results

getSedResults<-function(){
  out = list();
  
  #--get haul data with appended sediment characteristics (1975-2023)
  #--see r_Step1_GetSedDataAtHaulLocations.R
  lst = wtsUtilities::getObj("rda_Step1_SedDataAtAllNMFSHaulLocations.RData");
  
  
  #--get normalized temp, phi, and sorting as a function of depth for SBS years
  dfrHDp = lst$dfrHD |> sf::st_drop_geometry() |> 
                     dplyr::filter(dplyr::between(YEAR,2013,2018)) |> 
                     dplyr::select(BOTTOM_DEPTH,GEAR_TEMPERATURE,phi,sorting) |>
                     dplyr::rename(depth=BOTTOM_DEPTH,tmp=GEAR_TEMPERATURE,srt=sorting) |>
                     dplyr::mutate(dtmp=(tmp-min(tmp,na.rm=TRUE))/(max(tmp,na.rm=TRUE)-min(tmp,na.rm=TRUE)),
                                   dphi=(phi-min(phi,na.rm=TRUE))/(max(phi,na.rm=TRUE)-min(phi,na.rm=TRUE)),
                                   dsrt=(srt-min(srt,na.rm=TRUE))/(max(srt,na.rm=TRUE)-min(srt,na.rm=TRUE))) |>
                     dplyr::select(depth,tmp,phi,srt,dtmp,dphi,dsrt);
  dfrHDpp = dfrHDp |>
               dplyr::select(depth,dtmp,dphi,dsrt) |>
               dplyr::rename(`bottom temperature`=dtmp,phi=dphi,sorting=dsrt) |>
               tidyr::pivot_longer(c(`bottom temperature`,phi,sorting),
                                   names_to="variable",
                                   values_to="values");
  dfrHDp = dfrHDp |> dplyr::select(depth,tmp,phi,srt) |>
                     dplyr::rename(`bottom temperature`=tmp,phi=phi,sorting=srt);


  #--make pairs plots 
#| label: fig-SedRes-PairsPlot
  cap = "Pairs plot showing correlations among depth, bottom temperature, and interpolated mean grain size ('phi', in phi units) and sorting coefficient sediment characteristics at NMFS SBS haul locations (2013-2018)."
  p = GGally::ggpairs(dfrHDp);
  out = c(out,list(`fig-SedRes-PairsPlot`=list(p=p,cap=cap)));
  rm(dfrHDp,p);


#| label: fig-SedRes-TPSvD
  cap = "Relative distributions of bottom temperature, interpolated mean grain size ('phi', in phi units), and sorting coefficient with haul depth at NMFS SBS haul locations (2013-2018)."
  p = ggplot2::ggplot(dfrHDpp,mapping=aes(x=depth,y=values,colour=variable,fill=variable)) +
            ggplot2::geom_point(alpha=0.2) +
            ggplot2::geom_smooth(alpha=0.4) + 
            ggplot2::labs(x="bottom depth (m)",y="normalized value",colour="quantity",fill="quantity") +
            ggplot2::theme(legend.position="bottom");
  out = c(out,list(`fig-SedRes-TPSvD`=list(p=p,cap=cap)));
  rm(dfrHDpp,p);
  
  wtsUtilities::saveObj(out,"rda_Step2_SedDataResults.RData");
  return(out);
}
out = getSedResults();
rm(out);
