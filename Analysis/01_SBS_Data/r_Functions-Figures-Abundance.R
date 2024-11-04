#--various abundance plotting functions

require(ggplot2);
require(rlang);
  
plotTotRawZCs<-function(dfr_,sex_,col,label){
  p = ggplot(dfr_ |> dplyr::filter(sex==sex_),
             aes(x=size,y={{col}},colour=gear,linetype=type)) + 
        geom_line() + 
        geom_vline(xintercept=25,linetype=3) + 
        facet_wrap(~year,ncol=2,scales="free_y") + 
        labs(x="size (mm CW)",y=label) + 
        wtsPlots::getStdTheme() + 
        theme(legend.position=c(0.99,0.30),
              legend.justification=c(1,1),
              legend.box="horizontal",
              legend.title=element_blank());
  return(p);
}

  plotBZCs<-function(dfr,gear_types,sexs,factor_=x,xlims=c(25,180)){
    dfrb = dplyr::filter(dfr,gear %in% gear_types,x %in% sexs,bootrep>0)
    dfro = dplyr::filter(dfr,gear %in% gear_types,x %in% sexs,bootrep==0)
    p = ggplot(dfrb,
               aes(x=z,y=abd,colour={{factor_}},group=factor(paste0(y,"+",bootrep,x,gear)))) + 
          geom_line(linewidth=0.2,alpha=0.1) + 
          geom_line(data=dfro,linewidth=1) + 
          geom_line(data=dfro,colour="black",linewidth=0.1) + 
          geom_point(data=dfro,colour="black",size=0.5) + 
          facet_wrap(~y,scales="free_y",ncol=2) + 
          scale_x_continuous(limits=xlims) +
          labs(x="size (mm CW)",y="Size Composition (1,000s/sq. nmi)") +
          wtsPlots::getStdTheme() + 
          theme(legend.title=element_blank(),
                legend.position=c(0.99,0.30),
                legend.justification=c(1,1));
    return(p);
  }
  
  plotStatsBZCs<-function(dfrOrig,dfrBstats,gear_types,sexs,factor_=x,xlims=c(25,180)){
    dfrb = dplyr::filter(dfrBstats,gear %in% gear_types,x %in% sexs);
    dfro = dplyr::filter(dfrOrig,  gear %in% gear_types,x %in% sexs);
    p = ggplot(dfrb,
               aes(x=z,y=md,ymin=l90,ymax=u90,
                   colour={{factor_}},fill={{factor_}},
                   group=factor(paste0(y,"+",x,"+",gear)))) + 
          geom_ribbon(colour=NA,alpha=0.2) + 
          geom_line() + 
          geom_line(data=dfro,aes(x=z,y=abd,group=factor(paste0(y,"+",x,"+",gear))),
                    colour="black",linewidth=0.1,inherit.aes=FALSE) + 
          geom_point(data=dfro,aes(x=z,y=abd,group=factor(paste0(y,"+",x,"+",gear))),
                     colour="black",size=0.5,inherit.aes=FALSE) + 
          facet_wrap(~y,ncol=2,scales="free_y") + 
          scale_x_continuous(limits=xlims) +
          scale_y_log10() + 
          labs(x="size (mm CW)",y="Size Composition (1,000s/sq. nmi.)") + 
          wtsPlots::getStdTheme() + 
            theme(legend.title=element_blank(),
                  legend.position=c(0.99,0.30),
                  legend.justification=c(1,1));
    return(p);
  }
  
  plotCPUE<-function(dfr_,sex_,col,label){
    p = ggplot(dfr_ |> dplyr::filter(x==sex_),
               aes(x=z,y={{col}},colour=fleet,group=z)) + 
          geom_point(alpha=0.1,position="jitter") + 
          gghalves::geom_half_violin() + 
          geom_vline(xintercept=25,linetype=3) + 
          facet_wrap(~y,nrow=2,scales="free_y") + 
          labs(x="size (mm CW)",y=label) + 
          theme(legend.position.inside=c(0.99,0.48),
                legend.justification=c(1,1)) +
          wtsPlots::getStdTheme();
    return(p);
  }

plotStatsCPUE<-function(dfr=dfrStatsCPUE,mn,lci,uci){
  p = ggplot(dfr,aes(x=z,y={{mn}},ymin={{lci}},ymax={{uci}},colour=fleet,fill=fleet)) + 
        geom_vline(xintercept=25,linetype=3,alpha=0.6) + 
        geom_ribbon(colour=NA,alpha=0.2) + 
        geom_line() + 
        scale_x_continuous(limits=c(0,NA)) + 
        scale_y_log10() + 
        facet_wrap(~y,ncol=2,scales="free_y") + 
        labs(x="size (mm CW)",y="CPUE (no./sq. nmi.)") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="inside",
              legend.position.inside=c(0.99,0.30),
              legend.justification.inside=c(1,1),
              legend.title=element_blank())
  return(p)
}

plotCVsCPUE<-function(dfr,mn_=mn,sd_=sd,factor_=fleet,ylab="CV of CPUE",xlims=c(0,NA)){
  p = ggplot(dfr,aes(x=z,y={{sd_}}/{{mn_}},colour={{factor_}},fill={{factor_}})) + 
        geom_vline(xintercept=25,linetype=3,alpha=0.6) + 
        geom_hline(yintercept=1, linetype=3,alpha=0.6) + 
        geom_line() + 
        scale_x_continuous(limits=xlims) + 
        scale_y_continuous(limits=c(0,NA)) + 
        facet_wrap(~y,ncol=2,scales="fixed") + 
        labs(x="size (mm CW)",y=ylab) + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="inside",
              legend.position.inside=c(0.99,0.30),
              legend.justification.inside=c(1,1),
              legend.title=element_blank())
  return(p)
}

plotVMRsCPUE<-function(dfr,mn,vr){
  p = ggplot(dfr,aes(x=z,y={{vr}}/{{mn}},colour=fleet,fill=fleet)) + 
        geom_vline(xintercept=25,linetype=3,alpha=0.6) + 
        geom_hline(yintercept=1, linetype=3,alpha=0.6) + 
        geom_line() + 
        scale_x_continuous(limits=c(0,NA)) + 
        scale_y_log10() + 
        facet_wrap(~y,ncol=2,scales="free_y") + 
        labs(x="size (mm CW)",y="VMR") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="inside",
              legend.position.inside=c(0.99,0.30),
              legend.justification.inside=c(1,1),
              legend.title=element_blank())
  return(p)
}

plotPctNon0sCPUE<-function(dfr,factor=fleet){
  p = ggplot(dfr,aes(x=z,y=nn0/n,colour={{factor}},fill={{factor}})) + 
        geom_vline(xintercept=25,linetype=3,alpha=0.6) + 
        geom_line() + 
        scale_x_continuous(limits=c(0,NA)) + 
        scale_y_continuous(limits=c(0,1)) + 
        facet_wrap(~y,ncol=2,scales="free_y") + 
        labs(x="size (mm CW)",y="% non-0 hauls") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="inside",
              legend.position.inside=c(0.99,0.30),
              legend.justification.inside=c(1,1),
              legend.title=element_blank())
  return(p)
}

plotVarVsMean<-function(dfr=dfrStatsCPUE,facet=fleet,factor=x,label="sex"){
  p = ggplot(dfr, aes(x=log(mn),y=log(vr))) + 
        facet_wrap(vars({{facet}}),ncol=1) + 
        geom_abline(slope=1) + 
        geom_point(aes(shape={{factor}})) + 
        #geom_point(aes(shape={{factor}},size=n)) + 
        scale_size_area() + 
        labs(colour="size (mm CW)",shape=label,x="ln(mean)",y="ln(variance)",size="count") + 
        #ggnewscale::new_scale_colour() + 
        geom_smooth(aes(x=log(mn),y=log(vr),colour={{factor}},fill={{factor}}),
                    inherit.aes=FALSE,method="gam") + 
        geom_smooth(aes(x=log(mn),y=log(vr),colour={{factor}}),fill=NA,
                    inherit.aes=FALSE,method="lm",linetype=3) + 
        labs(colour=label,fill=label) + 
        wtsPlots::getStdTheme() + 
        theme(legend.position=c(0.01,0.99),
              legend.justification=c(0,1)); 
  return(p);
}
  