#--plot range of BSFRF haul-level CPUEs by size
require(ggplot2);

dirPrj = rstudioapi::getActiveProject();
dfrCPUE = wtsUtilities::getObj(file.path(dirPrj,"Analysis/01_SBS_Data/rda_Step3_SBS_CrabAbundance.RData"))$dfrCPUE;

dfrBSFRF = dfrCPUE |> dplyr::filter(fleet=="BSFRF",val>0) |> 
             dplyr::mutate(y=as.factor(y)) |> 
             dplyr::select(y, gis_station, sampling_factor, 
                           area_swept_variable, x, z, n, val);

p1 = ggplot(dfrBSFRF,aes(x=z,y=val,colour=y)) + 
      geom_point(alpha=0.2) + 
      facet_wrap(~x,ncol=1) + 
      wtsPlots::getStdTheme();
print(p1);

dfrStats = dfrBSFRF |> 
             dplyr::group_by(y,x,z) |> 
             dplyr::summarize(mn=mean(val),
                              sd=sd(val),
                              p05=quantile(val,probs=0.05),
                              p95=quantile(val,probs=0.95),
                              p10=quantile(val,probs=0.10),
                              p90=quantile(val,probs=0.90)) |> 
             dplyr::ungroup() |> 
             dplyr::mutate(cv=sd/mn,
                           r=sd^2/mn);
p2 = ggplot(dfrStats,aes(x=z,y=mn,ymin=p10,ymax=p90,colour=y,fill=y,group=y)) + 
      geom_ribbon(alpha=0.2,colour=NA) + 
      geom_line(alpha=1.0,linetype=2) + 
      geom_hline(yintercept=5000,linetype=3) + 
      geom_hline(yintercept=0,linetype=3) + 
      facet_wrap(~x,ncol=1) + 
      labs(x="size (mm CW)",y="estimated density (number/sq. nmi)",colour="year",fill="year") + 
      wtsPlots::getStdTheme();
print(p2);

p2a = ggplot(dfrStats,aes(x=z,y=mn,ymin=p10,ymax=p90,colour=y,fill=y,group=y)) + 
      geom_ribbon(alpha=0.2,colour=NA) + 
      geom_line(alpha=1.0,linetype=2) + 
      facet_wrap(~x,ncol=1) + 
      scale_y_log10() + 
      wtsPlots::getStdTheme();
print(p2a);

p3 = ggplot(dfrStats,aes(x=z,y=cv,colour=y,fill=y,group=y)) + 
      geom_line(alpha=1.0,linetype=1) + 
      facet_wrap(~x,ncol=1) + 
      wtsPlots::getStdTheme();
print(p3);

p4 = ggplot(dfrStats,aes(x=z,y=r,colour=y,fill=y,group=y)) + 
      geom_line(alpha=1.0,linetype=1) + 
      facet_wrap(~x,ncol=1) + 
      wtsPlots::getStdTheme();
print(p3);

