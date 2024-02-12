#--plot SBS empirical selectivity bootstrapping results

dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

# source(file.path(dirThs,"r_PlotEmpiricalSelectivity.R"));

load(file=file.path(dirThs,"rda_Step3a_EmpiricalSelectivityFromBootstrapping.RData"));

plotEmpSel<-function(dfr,colour,z_lim=c(22.5,187.5),n_min=0,wrap=TRUE){
  p = ggplot(dfr |> dplyr::filter(n_BSFRF+n_NMFS >= n_min),
             aes(x=z,y=emp_sel,colour={{colour}},fill={{colour}},group=paste(y,iB)));
  p = p + geom_line(alpha=0.1,size=0.1) + 
          geom_point(alpha=0.1,size=0.1,position=position_jitter(0.5)) + 
          geom_line(data=dfr |> dplyr::filter(iB==1),alpha=1,size=1);
  if (wrap){
    p = p + geom_smooth(mapping=aes(group=NULL),colour=NA,fill="gray25",
                        method="gam",formula=y~s(x,bs="cs",k=5));
  } else {
    p = p + geom_smooth(mapping=aes(x=z,y=emp_sel),colour=NA,fill="gray25",inherit.aes=FALSE,
                        method="gam",formula=y~s(x,bs="cs",k=5));
  }
  p = p + geom_hline(yintercept=c(0,1),linetype=2) + 
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

#--pdf("EmpiricalSelectivityFromBootstrapping.pdf",width=8,height=6);
#--females
dfrZCsp<-dfrZCs[(dfrZCs$x=="female")&dplyr::between(dfrZCs$z,25,125),];
dfrESsp<-dfrESs[(dfrESs$x=="female")&dplyr::between(dfrZCs$z,25,125),];
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(22.5,127.5),n_min=0,wrap=TRUE);
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(22.5,127.5),n_min=0,wrap=FALSE);
# plotEmpiricalSelectivity(dfrZCsp,dfrESsp,
#                           plotPoints=FALSE,
#                           points=list(alpha=0.2,size=0.5,dodge=2.5),
#                           plotLines=FALSE,
#                           plotViolins=TRUE,
#                           plotSmooths=TRUE,
#                           smooths=list(method="gam",formula=y~s(x,bs="cs",k=5),knots=c(25,50,75,100,125)));
#--males
dfrZCsp<-dfrZCs[(dfrZCs$x=="male")&(dfrZCs$z<=185),];
dfrESsp<-dfrESs[(dfrESs$x=="male")&(dfrESs$z<=185),];
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(22.5,187.5),n_min=0,wrap=TRUE);
plotEmpSel(dfrESsp,colour=factor(y),z_lim=c(22.5,187.5),n_min=0,wrap=FALSE);
# plotEmpiricalSelectivity(dfrZCsp,dfrESsp,
#                           plotPoints=FALSE,
#                           points=list(alpha=0.2,size=0.2,dodge=2.5),
#                           plotLines=FALSE,
#                           plotViolins=TRUE,
#                           plotSmooths=TRUE,
#                           smooths=list(method="gam",formula=y~s(x,bs="cs",k=7),knots=c(25,50,75,100,125,150,175)));
#--dev.off();
