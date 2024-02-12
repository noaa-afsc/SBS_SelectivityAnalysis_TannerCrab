plotEmpiricalSelectivity<-function(dfrZCs,
                                   dfrESs,
                                   plotPoints=FALSE,
                                   points=list(alpha=0.2,size=0.5,dodge=2.5),
                                   plotLines=FALSE,
                                   plotViolins=TRUE,
                                   violins=list(scale="width",stat="ydensity",dodge=2.5,alpha=0.3),
                                   plotSmooths=TRUE,
                                   smooths=list(method="auto",formula=y~x),
                                   showPlot=FALSE){
  require(ggplot2);
  uXs<-unique(dfrZCs$x);
  plots = list();
  #------estimated total abundance, by size
  for (uX in uXs){
    pByX = list()
    xlim<-c(25,185)
    if (uX=="female") xlim<-c(25,125);
    tmp <- dfrZCs[dfrZCs$x==uX,]
    pl<-wtsPlots::plotMDFR.XY(tmp,
                               x="z",
                               value.var="val",
                               colour="fleet",
                               fill="fleet",
                               plotPoints=plotPoints,
                               points=points,
                               plotLines=plotLines,
                               plotViolins=plotViolins,
                               violins=violins,
                               agg.formula="type+fleet+y+x+z+iB",
                               agg.function=wtsUtilities::Sum,
                               facet_grid=y~x,
                               scales="free_y",
                               xlim=xlim,
                               xlab="size (mm CW)",
                               ylab="abundance (millions)",
                               showPlot=FALSE) + wtsPlots::getStdTheme();
    pByX[["ZCs"]] = pl;
    if (showPlot) print(pl);
    pl <- pl + ggplot2::scale_y_log10();
    pl <- pl + ggplot2::ylab("abundance (millions, log10-scale)");
    pByX[["lnZCs"]] = pl;
    if (showPlot) print(pl);
    if (plotSmooths) {
      pl <- pl + ggplot2::scale_y_continuous(trans="identity");
      pl <- pl + ggplot2::ylab("abundance (millions)");
      pl <- pl + ggplot2::geom_smooth(alpha=0.3,linetype=2,
                                      method=smooths$method,formula=smooths$formula,
                                      mapping=ggplot2::aes_string(colour="fleet",fill="fleet"));
      print(pl);
      if (showPlot) pByX[["smZCs"]] = pl;
    }
  }
  #
  #------empirical selectivity, by size
  for (uX in uXs){
    tmp <- dfrESs[dfrESs$x==uX,];
    pl <- wtsPlots::plotMDFR.XY(tmp,
                                 x="z",
                                 value.var="emp_sel",
                                 colour="x",
                                 fill="x",
                                 plotPoints=plotPoints,
                                 points=points,
                                 plotLines=plotLines,
                                 plotViolins=plotViolins,
                                 violins=violins,
                                 facet_grid=y~type+x,
                                 ylim=c(0,1.2),
                                 xlab="size (mm CW)",
                                 ylab="empirical selectivity",
                                 guideTitleColour="sex",
                                 guideTitleFill="sex",
                                 showPlot=FALSE) + wtsPlots::getStdTheme();
    if (plotSmooths) pl <- pl + ggplot2::geom_smooth(alpha=0.3,linetype=2,
                                                     method=smooths$method,formula=smooths$formula,
                                                     colour="blue",fill="blue");
    pByX[["ESs"]] = pl;
    if (showPlot) print(pl);
    plots[[uX]] = pByX;
  }
  return(plots);
}
