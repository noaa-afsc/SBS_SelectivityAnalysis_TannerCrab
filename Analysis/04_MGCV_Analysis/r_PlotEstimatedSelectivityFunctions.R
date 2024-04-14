plotWgtdSels<-function(dfrRbyH,dfrMnR,n,lbl_size){
  iwt = paste0("invwgt",n);
  mnR = paste0("mnR",n);
  lwr = paste0("lower",n);
  upr = paste0("upper",n);
  p = ggplot(dfrRbyH,
        aes(x=z,y=prdR,size=1/.data[[iwt]])) + 
        geom_point() + 
        scale_size_area() + 
        geom_ribbon(aes(x=z,ymin=.data[[lwr]],ymax=.data[[upr]]),
                    data=dfrMnR,colour=NA,inherit.aes=FALSE,fill="red",alpha=0.2) + 
        geom_line(aes(x=z,y=.data[[mnR]]),
                  data=dfrMnR,colour="red",inherit.aes=FALSE) + 
        geom_hline(yintercept = c(0.5,1.0),linetype=3) + 
        scale_y_continuous(limits=c(0,2),oob=scales::squish) + 
        labs(x="size (mm CW)",y="predicted selectivity",size=lbl_size) + 
        wtsPlots::getStdTheme();
  return(p);
}

plotWgtdSelsWithY<-function(dfrRbyH,dfrMnR,n,lbl_size){
  iwt = paste0("invwgt",n);
  mnR = paste0("mnR",n);
  lwr = paste0("lower",n);
  upr = paste0("upper",n);
  p = ggplot(dfrRbyH,
        aes(x=z,y=prdR,size=1/.data[[iwt]])) + 
        geom_point() + 
        scale_size_area() + 
        geom_ribbon(aes(x=z,ymin=.data[[lwr]],ymax=.data[[upr]],fill=y),
                    data=dfrMnR,inherit.aes=FALSE,colour=NA,alpha=0.2) + 
        geom_line(aes(x=z,y=.data[[mnR]],colour=y),
                  data=dfrMnR,inherit.aes=FALSE) + 
        geom_hline(yintercept = c(0.5,1.0),linetype=3) + 
        scale_y_continuous(limits=c(0,2),oob=scales::squish) + 
        scale_fill_viridis_d(aesthetics=c("fill","colour")) +
        labs(x="size (mm CW)",y="predicted selectivity",size=lbl_size) + 
        wtsPlots::getStdTheme() +
        theme(axis.text.y=element_blank());
  return(p)
}

plotWgtdSelsByY<-function(dfrRbyH,dfrMnRbyY,n,lbl_size){
  iwt = paste0("invwgt",n);
  mnR = paste0("mnR",n);
  lwr = paste0("lower",n);
  upr = paste0("upper",n);
  ps = list();
  for (yd_ in levels(dfrMnRbyY$yd)){
    dfrMnRbyYp = dfrMnRbyY  |> dplyr::filter(yd==yd_);
    ps[[yd_]]=ggplot(dfrRbyH |> dplyr::filter(yd==yd_),
                     aes(x=z,y=prdR,ymin=lower,ymax=upper,size=1/.data[[iwt]])) + 
                geom_point() + 
                scale_size_area() + 
                geom_ribbon(aes(x=z,ymin=.data[[lwr]],ymax=.data[[upr]],fill=y),
                            data=dfrMnRbyYp,inherit.aes=FALSE,colour=NA,alpha=0.2) + 
                geom_line(aes(x=z,y=.data[[mnR]],colour=y),                
                          data=dfrMnRbyYp,inherit.aes=FALSE) + 
                geom_hline(yintercept = c(0.5,1.0),linetype=3) + 
                scale_y_continuous(limits=c(0,2),oob=scales::squish) + 
                scale_fill_manual(aesthetics=c("fill","colour"),values=viridis::inferno(10)) + 
                facet_wrap(~yd) + 
                labs(x="size (mm CW)",y="predicted selectivity",size=lbl_size) + 
                wtsPlots::getStdTheme() + 
                theme(legend.position="none",
                      axis.text.y=element_blank());
  }
  p = patchwork::wrap_plots(ps,nrow=2) + patchwork::plot_layout(axis_titles="collect");
  return(p)
}

