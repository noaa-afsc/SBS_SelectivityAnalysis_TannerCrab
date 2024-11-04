### Bootstrapping analysis

```{r SurvSel-GetBootstrapping,eval=FALSE}
  load(file.path(dirThs,"rda_Step3_EmpiricalSelectivityFromBootstrapping.RData"));
```

```{r SurvSel-DefFcn-PlotSurvSel,eval=FALSE}
plotSurvSel<-function(dfr,colour,z_lim=c(22.5,187.5),n_min=0,wrap=TRUE){
  p = ggplot(dfr |> dplyr::filter(n_BSFRF+n_NMFS >= n_min),
             aes(x=z,y=emp_sel,colour={{colour}},fill={{colour}},group=paste(y,iB)));
  p = p + geom_line(alpha=0.1,size=0.1) + 
          geom_point(alpha=0.1,size=0.1,position=position_jitter(0.5)) + 
          geom_line(data=dfr |> dplyr::filter(iB==1),alpha=1,linewidth=1);
  if (wrap){
    p = p + geom_smooth(mapping=aes(group=NULL),
                        colour="yellow",fill="yellow",alpha=1,linewidth=1.0,
                        method="gam",formula=y~s(x,bs="cs",k=7));
  } else {
    p = p + geom_smooth(mapping=aes(x=z,y=emp_sel),
                        colour="yellow",fill="yellow",alpha=1,linewidth=1.5,
                        method="gam",formula=y~s(x,bs="cs",k=7),
                        inherit.aes=FALSE);
  }
  p = p + geom_hline(yintercept=c(0,1),linetype=2) + 
          geom_hline(yintercept=c(0.5),linetype=3) + 
          labs(x="size (mm CW)",y="empirical selectivity",colour="year",fill="year") + 
          scale_y_continuous(limits=c(0,3),oob=scales::squish) + 
          scale_x_continuous(limits=z_lim,oob=scales::squish) + 
          scale_colour_discrete(aesthetics=c("colour","fill")) + 
          guides(colour=guide_legend(override.aes=list(alpha=1))) + 
          wtsPlots::getStdTheme() + 
          theme(panel.grid.major=element_line(colour="gray75"),
                panel.grid.minor=element_line(colour="white"));
  if (wrap) p = p + facet_wrap(~y,nrow=2);
  return(p);
}
```

```{r, eval=FALSE}
#| label: fig-SurvSel-PlotBESsF
#| fig-cap: "Boostrapped empirical selectivities for female Tanner crab in the SBS studies. Upper plot: Annual SBS-derived selectivities (heavy line is mean value, yellow curve is smoothed fit); lower plot: all SBS years combined, yellow curve is smoothed fit to all."
  #--females
  dfrESsp<-dfrBESs |> dplyr::filter((x=="female"),dplyr::between(z,25,125));
  p1 = plotSurvSel(dfrESsp,colour=factor(y),z_lim=c(22.5,127.5),n_min=0,wrap=TRUE) + 
         theme(legend.position="none",axis.title.x=element_blank());
  p2 = plotSurvSel(dfrESsp,colour=factor(y),z_lim=c(22.5,127.5),n_min=0,wrap=FALSE) +
         theme(legend.position=c(0.01,0.99),legend.justification=c(0,1));
  pg = cowplot::plot_grid(p1,p2,ncol=1);
  lstFigs = c(lstFigs,wtsQMD::printGGplot(pg));
  rm(dfrESsp,p1,p2,pg);
```

```{r, eval=FALSE}
#| label: fig-SurvSel-PlotBESsM
#| fig-cap: "Boostrapped empirical selectivities for male Tanner crab in the SBS studies. Upper plot: Annual SBS-derived selectivities (heavy line is mean value, yellow is smoothed fit); lower plot: all SBS years combined, yellow is smoothed fit to all."
  #--males
  dfrESsp<-dfrBESs |> dplyr::filter((x=="male"),dplyr::between(z,25,185));
  p1 = plotSurvSel(dfrESsp,colour=factor(y),z_lim=c(22.5,187.5),n_min=0,wrap=TRUE) + 
         theme(legend.position="none",axis.title.x=element_blank());
  p2 = plotSurvSel(dfrESsp,colour=factor(y),z_lim=c(22.5,187.5),n_min=0,wrap=FALSE) +
         theme(legend.position=c(0.01,0.99),legend.justification=c(0,1));
  pg = cowplot::plot_grid(p1,p2,ncol=1);
  lstFigs = c(lstFigs,wtsQMD::printGGplot(pg));
  rm(dfrESsp,p1,p2,pg);
```

```{r SurvSel-CalcStats, eval=FALSE}
#--calculate bootstrap means by sex, year, size
dfrMeans<-reshape2::dcast(dfrBESs,"x+y+z~.",fun.aggregate=mean,na.rm=TRUE,value.var="emp_sel");
names(dfrMeans)[4]<-"mean";
dfrStds <-reshape2::dcast(dfrBESs,"x+y+z~.",fun.aggregate=var, na.rm=TRUE,value.var="emp_sel");
dfrStds[["."]]<-sqrt(dfrStds[["."]]);
names(dfrStds)[4]<-"stddev";
cis<-wtsUtilities::calcCIs(dfrMeans$mean,sdvs=dfrStds$stddev,pdfType="normal",ci=0.8);
dfrStats<-cbind(dfrMeans,stddev=dfrStds[["stddev"]],lci=cis$lci,uci=cis$uci);
dfrStats$y<-factor(dfrStats$y)
```

```{r, eval=FALSE}
#| label: fig-SurvSel-BStats-MnCIs
#| fig-cap: "Mean and 80% confidence intervals for annual empirical selectivity from bootstrapping analysis. Upper plot: females; lower plot: males. Results for different SBS study years are shown with different colors."
  pl<-wtsPlots::plotMDFR.XY(dfrStats,x="z",value.var="mean",
                            plotLines=TRUE, plotPoints=FALSE, 
                            colour="y", facet_grid="x~.", showPlot=FALSE)
  pl<-pl + ggplot2::geom_ribbon(mapping=aes(ymin=lci,ymax=uci,fill=y),alpha=0.4,colour=NA);
  pl<-pl + ggplot2::labs(x="size (mm CW)",y="empirical selectivity",colour="year",fill="year");
  pl<-pl + ggplot2::geom_hline(yintercept=c(0,1),linetype=2,colour="grey50");
  pl<-pl + ggplot2::scale_y_continuous(limits=c(-0.5,2),oob=scales::squish);
  pl<-pl+wtsPlots::getStdTheme();
  lstFigs = c(lstFigs,wtsQMD::printGGplot(pl));
```
