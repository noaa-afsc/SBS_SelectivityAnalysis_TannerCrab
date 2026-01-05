  #--load packages
  require(ggplot2);
  require(spatstat);#--note that this loads a number of spatstat-related packages
  require(rSimSpatialProcessSampling);

plotSpatialSimSBS<-function(i=100, #--density (number/km^2)
                            s=0.5, #--NMFS selectivity
                            cases=NULL,
                            nreps=1,
                            lstSamplingInfo=NULL,
                            seed=NULL){
  #--set random seed
  set.seed(seed);

  #--load SBS sampling info from rSimSpatialProcessSampling
  if (is.null(lstSamplingInfo)) 
    eval(data("lstSamplingInfo",package="rSimSpatialProcessSampling",envir=));
  
  tdBSFRF = lstSamplingInfo$tdBSFRF;#--typical tow distance
  nwBSFRF = lstSamplingInfo$nwBSFRF;#--typical net width
  asBSFRF = lstSamplingInfo$asBSFRF;#--typical area swept
  tdNMFS = lstSamplingInfo$tdNMFS;#--typical tow distance
  nwNMFS = lstSamplingInfo$nwNMFS;#--typical net width
  asNMFS = lstSamplingInfo$asNMFS;#--typical area swept
  offset = lstSamplingInfo$offset;#--offset distance between haul tracks
  
  # cat("tdBSFRF =",tdBSFRF,"\n")
 
  #----set the NMFS tow distance (window width)
  wid = tdNMFS;                
  #----set the max width (NMFS+BSFRF net widths + offset distance; window height)
  hgt = nwNMFS+offset+nwBSFRF;
  #--set the sampling "window"s
  win_all   = spatstat.geom::owin(xrange=c(0,wid),    yrange=              c(0,hgt),    unitname="km");
  win_NMFS  = spatstat.geom::owin(xrange=c(0,tdNMFS), yrange=              c(0,nwNMFS), unitname="km");
  win_BSFRF = spatstat.geom::owin(xrange=c(0,tdBSFRF),yrange=nwNMFS+offset+c(0,nwBSFRF),unitname="km");
  owin2tbbl<-function(win,type){
    tbbl = tibble::tibble_row(type=type,
                              x=NA,y=NA,
                              xmin=win$xrange[1],xmax=win$xrange[2],
                              ymin=win$yrange[1],ymax=win$yrange[2]);
    return(tbbl);
  }
  dfrWins = dplyr::bind_rows(owin2tbbl(win_all,"rect_seascape"),
                             owin2tbbl(win_NMFS,"rect_NMFS"),
                             owin2tbbl(win_BSFRF,"rect_BSFRF"));
              

  #--generate spatial patterns
  lst = list();
  if (is.null(cases)) {cases = tibble::tibble(i=i,s=s);}
  for (nc in 1:nrow(cases)){
    case = cases[nc,];
    i = case$i; s = case$s;
    case_ = paste0("CPUE = ",i,"/km^2; NMFS Sel = ",s);
    for (rep in 1:nreps){
      pat = spatstat.random::rpoispp(i,win=win_all,nsim=1,drop=TRUE);
      #--do BSFRF sampling
      as_BSFRF = asBSFRF;
      in_BSFRF = inside.owin(pat$x,pat$y,win_BSFRF); #--number inside area swept
      n_BSFRF  = sum(in_BSFRF);
      #--do NMFS sampling, including selectivity ratio
      as_NMFS = asNMFS;
      in_NMFS1 = inside.owin(pat$x,pat$y,win_NMFS); #--number inside area swept
      n_NMFS1  = sum(in_NMFS1);
      in_NMFS2 = rbinom(n_NMFS1,1,s)==TRUE;        #--apply selectivity
      n_NMFS2  = sum(in_NMFS2);
      #--put it all together
      dfrA  = tibble::tibble(type="seascape",x=pat$x,y=pat$y);
      dfrB  = tibble::tibble(type="BSFRF",x=pat$x[in_BSFRF],y=pat$y[in_BSFRF]);
      dfrNa = tibble::tibble(type="NMFS available",x=pat$x[in_NMFS1],y=pat$y[in_NMFS1]);
      if (n_NMFS2>0) {
        dfrNs = dfrNa[in_NMFS2,]; dfrNs$type="NMFS selected";
      } else {
        dfrNs = tibble::tibble(type="NMFS selected",x=NA,y=NA);
      }
      idx = paste0(rep," BSFRF: ",n_BSFRF,"\n NMFS: ",n_NMFS2);
      dfr = dplyr::bind_rows(dfrWins,dfrA,dfrB,dfrNa,dfrNs) |> 
              dplyr::mutate(case=case_,idx=idx);
      lst[[paste(case_,rep)]] = dfr;
    }
  }
  dfrAll = dplyr::bind_rows(lst); 
  rm(lst);
  plt = ggplot(dfrAll) + 
          geom_rect(data=dfrAll |> dplyr::filter(type=="rect_seascape"),
                    aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    alpha=0.2,colour="black",fill="gray50") +
          geom_rect(data=dfrAll |> dplyr::filter(type=="rect_NMFS"),
                    aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    alpha=0.2,colour="blue",fill="blue") +
          geom_rect(data=dfrAll |> dplyr::filter(type=="rect_BSFRF"),
                    aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
                    alpha=0.2,colour="green",fill="green") +
          geom_point(aes(x=x,y=y),dfrAll |> dplyr::filter(type=="seascape"),size=0.2,colour="black") + 
          geom_point(aes(x=x,y=y),dfrAll |> dplyr::filter(type=="NMFS available"),size=0.2,colour="blue") + 
          geom_point(aes(x=x,y=y),dfrAll |> dplyr::filter(type=="NMFS selected"),size=0.4,shape=21,colour="cyan",fill="cyan") + 
          geom_point(aes(x=x,y=y),dfrAll |> dplyr::filter(type=="BSFRF"),size=0.4,shape=21,colour="green",fill="green") + 
          labs(x="along tow (km)",y="cross tow (km)") + 
          scale_x_continuous(expand=c(0,0)) +
          scale_y_continuous(expand=c(0,0)) +
          facet_grid(vars(idx),vars(case),drop=TRUE) + 
          wtsPlots::getStdTheme() + 
          theme(aspect.ratio=2*win_all$yrange[2]/win_all$xrange[2],
                panel.spacing.x=unit(1,"mm"),
                panel.spacing.y=unit(1,"mm"));
  return(plt);
}
  plotSpatialSimSBS(i=10,s=0.5,nreps=5,lstSamplingInfo=lstSamplingInfo);
if (FALSE){
  ggsave("fig_PairedTowExample.pdf",width=9,height=4)  
  ggsave("fig_PairedTowExample.png",width=9,height=4) 
  }
if (FALSE){
  data("lstSamplingInfo",package="rSimSpatialProcessSampling");
  thmx = theme(axis.title=element_blank(),plot.margin=margin(0,0,0,0),
               strip.text=element_text(size=rel(0.5)));
  p1 = plotSpatialSimSBS(i= 100,s=0.5,nreps=5,lstSamplingInfo=lstSamplingInfo); 
  p2 = plotSpatialSimSBS(i= 100,s=0.9,nreps=5,lstSamplingInfo=lstSamplingInfo);
  p3 = plotSpatialSimSBS(i=1000,s=0.5,nreps=5,lstSamplingInfo=lstSamplingInfo);
  p4 = plotSpatialSimSBS(i=1000,s=0.9,nreps=5,lstSamplingInfo=lstSamplingInfo);
  cowplot::plot_grid(p1+thmx,p2+thmx,p3+thmx,p4+thmx,
                     nrow=1,rel_widths=c(1,1,1,1));
  ggsave("fig_PairedTowExamples.pdf",width=9,height=0.376*9)  
  ggsave("fig_PairedTowExamples.png",width=9,height=0.376*9) 
  cases = tibble::tribble(  ~i, ~s,
                           100,0.5,
                           100,0.9,
                          1000,0.5,
                          1000,0.9);
  plotSpatialSimSBS(cases=cases,nreps=5,lstSamplingInfo=lstSamplingInfo);
}  
  
  
  