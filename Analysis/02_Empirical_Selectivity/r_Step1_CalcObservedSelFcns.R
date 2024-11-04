#--calculate empirical selectivity from size composition results

#--list for output results
out = list();

#--get list object with size composition results
lst = wtsUtilities::getObj(file.path("../01_SBS_Data/rda_Step5_BootZCs_AllResults.RData"));

#--calculate empirical selectivity from observed and bootstrapped size comps
dfrBootEmpSels = lst$dfrBootZCs |> dplyr::select(x,y,z,bootrep,gear,abd) |>
                   tidyr::pivot_wider(names_from="gear",values_from=abd) |> 
                   dplyr::mutate(emp_sel=NMFS/BSFRF,
                                 emp_prp=NMFS/(NMFS+BSFRF));

plotEmpSel<-function(dfr,qty,sex,ylab="Empirical Selectivity (NMFS/BSFRF)",xlim=185,ylim=1){
  require(ggplot2)
  p = ggplot(dfr |> dplyr::filter(x %in% sex,bootrep>0,is.finite({{qty}})),
             aes(x=z,y={{qty}},colour=x,group=paste0(y,"+",x,"+",bootrep))) + 
        geom_line(linewidth=0.2,alpha=0.1) + 
        geom_line(data=dfr |> dplyr::filter(x %in% sex,bootrep==0),
                  mapping=aes(x=z,y={{qty}},group=paste0(y,"+",x,"+",bootrep)),colour="black",inherit.aes=FALSE) + 
        geom_point(data=dfr |> dplyr::filter(x %in% sex,bootrep==0),
                   mapping=aes(x=z,y={{qty}},group=paste0(y,"+",x,"+",bootrep)),colour="black",inherit.aes=FALSE) + 
        geom_hline(yintercept=c(0.5,1.0),linetype=3,alpha=0.5) + 
        scale_x_continuous(limits=c(0,xlim)) + 
        scale_y_continuous(limits=c(0,ylim),oob=scales::squish) + 
        labs(x="size (mm CW)",y=ylab) + 
        facet_wrap(~y,ncol=2,scales="free_y") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="none");
  return(p);
}

plotEmpSel(dfrBootEmpSels,emp_sel,"male",  xlim=185,ylim=3);
plotEmpSel(dfrBootEmpSels,emp_sel,"female",xlim=135,ylim=3);
plotEmpSel(dfrBootEmpSels,emp_prp,"male",  xlim=185,ylim=1,"Empirical Proportion (NMFS/(NMFS+BSFRF))");
plotEmpSel(dfrBootEmpSels,emp_prp,"female",xlim=135,ylim=1,"Empirical Proportion (NMFS/(NMFS+BSFRF))");

  #--compute statistics for bootstrapped empriical selectivities and proportions----
calcStats<-function(qty){
  dfr = dfrBootEmpSels |> dplyr::filter(bootrep>0) |> 
                   dplyr::group_by(x,y,z) |> 
                   dplyr::summarize(mn=mean({{qty}},na.rm=TRUE),
                                    md=median({{qty}},na.rm=TRUE),
                                    vr=var({{qty}},na.rm=TRUE),
                                    sd=sqrt(vr),
                                    l90=quantile({{qty}},0.05,na.rm=TRUE),
                                    u90=quantile({{qty}},0.95,na.rm=TRUE)) |> 
                   dplyr::ungroup();
  return(dfr);
}
dfrStatsEmpSel = calcStats(emp_sel);#--lots of infinities
dfrStatsEmpPrp = calcStats(emp_prp);

#--plot stats
plotStats<-function(sex,xlim=185){
  require(ggplot2)
  p = ggplot(dfrStatsEmpPrp |> dplyr::filter(x=={{sex}}),aes(x=z)) + 
        geom_ribbon(aes(ymin=l90,ymax=u90),alpha=0.4) + 
        geom_line(aes(y=md)) + 
        geom_line(aes(y=mn),linetype=3) + 
        geom_hline(yintercept=0.5,linetype=2) + 
        facet_wrap(~y,ncol=2,scales="free_y") + 
        scale_x_continuous(limits=c(25,xlim)) +
        scale_y_continuous(limits=c(0,1)) + 
        labs(x="size (mm CW)",y="Proportion") + 
        wtsPlots::getStdTheme() + 
          theme(legend.title=element_blank(),
                legend.position=c(0.99,0.30),
                legend.justification=c(1,1));
  return(p)
}
plotStats("male",xlim=185)
plotStats("female",xlim=135)

#--weights??: number of crab sampled in BSFRF