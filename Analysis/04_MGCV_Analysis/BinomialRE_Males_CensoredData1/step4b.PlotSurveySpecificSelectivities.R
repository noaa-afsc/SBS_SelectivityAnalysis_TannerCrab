#--plot predicted catchabilities/catchability ratios

#--load required packages
require(ggplot2);
require(magrittr);
#--other packages
#----cowplot
#----dplyr
#----ggdist
#----wtsUtilities

#--set up plotting output
device = "png";
pltctr = 1;
plotFN<-function(n,dev=device){paste0("plots4b_Males_",wtsUtilities::formatZeros(n),".",dev)}

#--plot haul-specific R's
dfrRbyH = wtsUtilities::getObj("./results_RData/dfrRbyH.RData");#--load previously-calculated results
# p = ggplot(dfrRbyH,mapping=aes(x=z,y=prdR)) + 
#       geom_point(position=position_jitter(width=1)) +
#       geom_line(mapping=aes(group=h),alpha=0.2) +
#       labs(x="size (mm CW)",y="predicted haul-specific  R");
# print(p);
# ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;

# --plot results by haul in groups of years
seFactor = 2;
intv = 6;#--number of years to plot in groups
plts1 = list();
for (iy in seq(1982,2019,intv)){
  ylvls = iy+(1:intv)-1;
  tmp1 = dfrRbyH %>% dplyr::filter(dplyr::between(y,min(ylvls),max(ylvls)));#,(prdR<10));
  tmp2 = tmp1 %>% subset(!(is.na(t)|is.na(f)|is.na(s))) %>%
                  dplyr::arrange(y,z,h) %>%
                  dplyr::group_by(y,z) %>% 
                  dplyr::summarize(mnR=mean(prdR,na.rm=TRUE),
                                   seR=sd(prdR,na.rm=TRUE),
                                   lower=mnR-seFactor*seR,
                                   upper=mnR+seFactor*seR) %>%
                  dplyr::ungroup() %>%
                  dplyr::select(y,z,mnR,seR,lower,upper);
  tmp1$y = factor(tmp1$y,levels=ylvls);
  tmp2$y = factor(tmp2$y,levels=ylvls);
  #--plot annual prdR averaged over haul-specific prdRs for each year (no weighting)
  p1 = ggplot2::ggplot(data=tmp2,mapping=aes(x=z,y=mnR,ymin=lower,ymax=upper,colour=y,fill=y)) +
        ggdist::geom_lineribbon(alpha=0.3) + ggplot2::geom_line(size=1) +
        #scale_colour_viridis_d(option="plasma",aesthetics=c("colour","fill")) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(y="mean predicted catchability (ratio)",colour="year",fill="year")+
        theme(legend.position = "bottom");
  pL = cowplot::get_legend(p1); p1 = p1+theme(legend.position = "none");
  p2 = ggplot2::ggplot() +
        ggplot2::geom_point(data=tmp1,mapping=aes(x=z,y=prdR,colour=y,group=h),alpha=0.3) +
        ggplot2::geom_line(data=tmp1,mapping=aes(x=z,y=prdR,colour=y,group=h),alpha=0.3) +
        ggdist::geom_lineribbon(data=tmp2,mapping=aes(x=z,y=mnR,ymin=lower,ymax=upper),colour="black",fill="black",alpha=0.3) +
        ggplot2::geom_line(data=tmp2,mapping=aes(x=z,y=mnR,colour=y),colour="black",alpha=0.3,size=1) +
        #scale_colour_viridis_d(option="plasma",aesthetics=c("colour","fill")) +
        coord_cartesian(ylim=c(0,NA)) +
        labs(x="size (mm CW)",y="predicted catchability (ratio)",colour="year",fill="year") +
        ggplot2::facet_wrap(vars(y),ncol=2) +
        theme(legend.position="none");
  pg = cowplot::plot_grid(cowplot::plot_grid(pL,p1,ncol=1,rel_heights=c(1.5,10)),
                          p2,rel_widths = c(5,5),nrow=1);
  plts1[[paste0(min(ylvls),":",max(ylvls))]] = pg;
  #print(pg);   
  #ggsave(plotFN(pltctr),width=11,height=6); pltctr %<>% +1;
}

#--plot year-specific survey catchability by multi-year intervals
dfrRbyY = wtsUtilities::getObj("./results_RData/dfrRbyY.RData");
intv = 6;#--number of intervals per plot
tmp1 = dfrRbyY %>% 
         dplyr::mutate(group=floor((y-min(y))/intv+1),
                       lbl=paste0((group-1)*intv+min(y)," : ",group*intv+min(y)-1),
                       numI=ifelse(numIndivs>4,">4",ifelse(numIndivs==0,"=0","1-4")));
tmp1$numI = factor(tmp1$numI,levels=c("=0","1-4",">4"));
ymin = min(tmp1$y);
uGs = unique(tmp1$group);
plts2 = list();
for (uG in uGs){
  tmp = tmp1 %>% subset(group==uG);
  tmp$y = as.factor(tmp$y);
  p = ggplot2::ggplot() +
        ggdist::geom_lineribbon(data=tmp,mapping=aes(x=z,y=mnR,ymin=lower,ymax=upper,colour=y,fill=y),alpha=0.3) +
        ggplot2::geom_line(data=tmp,mapping=aes(x=z,y=mnR,colour=y)) + 
        ggplot2::geom_point(data=tmp,mapping=aes(x=z,y=mnR,shape=numI,colour=y),size=3,alpha=1.0) +
        ggplot2::facet_grid(lbl~.) +
        ggplot2::coord_cartesian(ylim=c(0,NA)) +
        ggplot2::labs(x="size (mm CW)",y="CPUE/inverse variance-averaged\nannual catchability (ratio)",
                      colour="year",fill="year",shape="# crab");
  # print(p);
  # ggsave(plotFN(pltctr),width=11,height=6); pltctr %<>% +1;
  plts2[[paste0((uG-1)*intv+ymin,":",uG*intv+ymin-1)]]=p;
}

for (i in 1:length(plts1)){
  p1 = plts1[[i]];
  p2 = plts2[[i]];
  pg = cowplot::plot_grid(p1,p2,ncol=1);
  #print(pg);
  ggsave(plotFN(pltctr),width=8,height=6); pltctr %<>% +1;
}

tmp = dfrRbyY %>% dplyr::mutate(numI=ifelse(numIndivs>4,">4",ifelse(numIndivs==0,"=0","1-4")))
p =  ggplot2::ggplot(tmp,aes(x=z,y=mnR,colour=as.factor(y),shape=numI)) + 
        ggplot2::geom_line(alpha=0.3) +
        ggplot2::geom_point(size=3,alpha=1.0) +
        ggplot2::coord_cartesian(ylim=c(0,NA)) +
        ggplot2::labs(x="size (mm CW)",y="CPUE/inverse variance-averaged\nannual catchability (ratio)",
                      colour="year",fill="year",shape="# crab") +
       ggplot2::guides(shape=guide_legend(direction="horizontal"),
                       colour=guide_legend(direction="vertical",ncol=3));
ggsave(paste0("plot4b_Males_AllAnnualCatchabilityCurves.",device),plot=p,width=8,height=5);


