#--export sel curves to csv files----
require(ggplot2);

#--define grid for output----
##--grid for years----
grd_y = c(1982:2019,2021:2023);
dfrY = tibble::tibble(y=grd_y);
##--grid for sizes----
grd_z = seq(27.5,182.5,5);
dfrZ = tibble::tibble(z=grd_z);
##--combined grid----
dfrYZ = dfrY |> dplyr::cross_join(dfrZ);

#--mean selectivity curve----
dfrMnSel = wtsUtilities::getObj("rda_GaussianModels_BestEstMnSels.RData") |> 
             dplyr::select(z,est=mnR3,lb=lower3,ub=upper3);
selMxZ = (dfrMnSel |> dplyr::filter(z==max(z)))[["est"]];

p = ggplot(dfrMnSel,aes(x=z,y=est)) + 
      geom_ribbon(aes(ymin=lb,ymax=ub),alpha=0.3) + 
      geom_line() + 
      scale_x_continuous(limits=c(NA,max(grd_z))) + 
      wtsPlots::getStdTheme();
print(p);

dfrMnSlx = dfrYZ |> dplyr::left_join(dfrMnSel |> dplyr::select(z,est),
                                     by=c("z")) |> 
             dplyr::mutate(est=ifelse(is.na(est),selMxZ,est));
p = p + geom_line(data=dfrMnSlx,colour="blue");
print(p);

dfrExp = dfrMnSlx |> 
           tidyr::pivot_wider(id_cols=c("y"),names_from="z",values_from="est") |> 
           dplyr::mutate(`#--year`=paste("#--",y)) |> 
           dplyr::select(!y);
readr::write_delim(dfrExp,"txt_BestEstSels-MnR3.txt",quote="none");

#--selectivity curves for all years----
dfrSels = wtsUtilities::getObj("rda_GaussianModels_BestEstSelsByYr.RData") |> 
             dplyr::select(y,z,est=mnR3,lb=lower3,ub=upper3);
selMxZ = dfrSels |> 
            dplyr:: group_by(y) |> 
            dplyr::summarize(max=est[z==max(z)]) |> 
            dplyr::ungroup() |> 
            dplyr::mutate(y=as.numeric(levels(y)))

p = ggplot(dfrSels,aes(x=z,y=est,colour=y)) + 
      geom_ribbon(aes(ymin=lb,ymax=ub,fill=y),alpha=0.3) + 
      geom_line() + 
      geom_line(data=dfrMnSel,mapping=aes(x=z,y=est),color="black",inherit.aes=FALSE) + 
      scale_x_continuous(limits=c(NA,max(grd_z))) + 
      wtsPlots::getStdTheme() + 
      theme(legend.position="none");
print(p);

dfrSlx = dfrYZ |> dplyr::left_join(dfrSels |> dplyr::select(y,z,est) |> 
                                     dplyr::mutate(y=as.numeric(as.character(y))),
                                     by=c("y","z")) |> 
             dplyr::left_join(selMxZ,by=c("y")) |>
             dplyr::mutate(est=ifelse(is.na(est),max,est)) |> 
             dplyr::select(!max);
p = p + geom_line(data=(dfrSlx |> dplyr::mutate(y=factor(y))));
print(p);

dfrExp = dfrSlx |> 
           tidyr::pivot_wider(id_cols=c("y"),names_from="z",values_from="est") |> 
           dplyr::mutate(`#--year`=paste("#--",y)) |> 
           dplyr::select(!y);
readr::write_delim(dfrExp,"txt_BestEstSels-AllYrsR3.txt",quote="none");




