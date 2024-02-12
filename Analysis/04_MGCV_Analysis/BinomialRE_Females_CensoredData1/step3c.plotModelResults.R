#--plot model results

#--load required packages
require(dplyr);
require(ggplot2);
require(magrittr);
require(mgcv);
require(mgcViz);
#--other required packages
#----ggdist
#----wtsUtilities


#--get model results
x = "FEMALE"; maxZ = 130;
mdl = wtsUtilities::getObj(paste0("./results_RData/selModels.BinomRE.Best.",tolower(x),"s.RData"));

#--get original data used in model fit
dfrDat = wtsUtilities::getObj("./results_RData/dfrDatRE.RData");

#--set up plotting output
device = "png";
pltctr = 1;
plotFN<-function(n,dev=device){paste0("plots3c_females_",wtsUtilities::formatZeros(n),".",dev)}

#--set up covariates for predictions
grd_z = seq(5,maxZ,5)+2.5; med_z = 77.5;
med_d = median(dfrDat$d,na.rm=TRUE); rng_d = range(dfrDat$d,na.rm=TRUE); grd_d = seq(from=rng_d[1],rng_d[2],length.out=50);
med_t = median(dfrDat$t,na.rm=TRUE); rng_t = range(dfrDat$t,na.rm=TRUE); grd_t = seq(from=rng_t[1],rng_t[2],length.out=50);
med_f = median(dfrDat$f,na.rm=TRUE); rng_f = range(dfrDat$f,na.rm=TRUE); grd_f = seq(from=rng_f[1],rng_f[2],length.out=50);
med_s = median(dfrDat$s,na.rm=TRUE); rng_s = range(dfrDat$s,na.rm=TRUE); grd_s = seq(from=rng_s[1],rng_s[2],length.out=50);
grds = list(z=grd_z,d=grd_d,t=grd_t,f=grd_f,s=grd_s);
newDFR<-function(vars,grds){
  #str = paste0(vars,"=grds[['",vars,"']]",collapse=",")
  gnms = names(grds);
  for (i in 1:length(gnms)){
    var = gnms[i];
    if (i==1){
      if (var %in% vars) 
        {str = paste0(var,"=grds[['",var,"']]");}
      else
        {str = paste0(var,"=NA");}
    } else {
      if (var %in% vars) 
        {str = paste0(str,",",var,"=grds[['",var,"']]");}
      else
        {str = paste0(str,",",var,"=NA");}
    }
  }
  eval(parse(text=paste0("dfrPrd1 = tidyr::expand_grid(",str,")")));
  return(dfrPrd1);
}

#--get model intercept
int = mdl$coefficients[1];

#--get model smooths
smths = mdl$smooth;
ns    = length(smths);

# smth = smths[[1]];
# cs = mdl$coefficients[smth$first.para:smth$last.para];
# dfrPrd1 = newDFR("z",grds);
# pm = mgcv::PredictMat(smth,dfrPrd1);#--prediction matrix (values of basis functions for each row of data)
# prd = int + pm %*% cs;
# dfrPrd2 = tibble::tibble(z=grd_z,`ln(R)`=prd,R=exp(prd));
# ggplot(data=dfrPrd2,mapping=aes(x=z,y=R)) +
#   geom_line();
# 
# smth = smths[[4]];
# trms = smth$term;
# cs = mdl$coefficients[smth$first.para:smth$last.para];
# dfrPrd1 = newDFR(trms,grds);
# pm = mgcv::PredictMat(smth,dfrPrd1);#--prediction matrix (values of basis functions for each row of data)
# prd = pm %*% cs;
# dfrPrd2 = cbind(dfrPrd1,tibble::tibble(R=prd));
# p = ggplot(dfrPrd2,mapping=aes(x=z,y=d,z=R,fill=R)) + geom_raster() + scale_fill_viridis_c(option="plasma");
# print(p)
# require(rayshader)
# plot_gg(p)
# 
# ps = list();   #--plot list
# dfrPrds = NULL;#--dataframe with predicted values
# for (s in 1:ns){
#   smth = smths[[s]];
#   if (!inherits(smth,"random.effect")){
#     trms = smth$term;
#     lbl  = smth$label;
#     cs   = mdl$coefficients[smth$first.para:smth$last.para];
#     dfrPrd1 = newDFR(trms,grds);
#     pm = mgcv::PredictMat(smth,dfrPrd1);#--prediction matrix (values of basis functions for each row of data)
#     prd = pm %*% cs;                    #--predicted logit-scale value for each row of data
#     dfrPrd2 = cbind(dfrPrd1,tibble::tibble(lgt_val=prd,label=lbl,terms=paste0(terms,collapse=",")));
#     dfrPrds = rbind(dfrPrds,dfrPrd2);
#     if (length(trms)==1){
#       p = ggplot(data=dfrPrd2,mapping=aes_string(x=trms[1],y="lgt_val")) +
#             geom_line() +
#             labs(y="ln-scale additive effect",subtitle=lbl);
#     } else {
#       p = ggplot(data=dfrPrd2,mapping=aes_string(x=trms[1],y=trms[2],fill="lgt_val")) +
#             geom_raster() + scale_fill_viridis_c(option="plasma") +
#             labs(fill="ln-scale\nadditive effect",subtitle=lbl);
#     }
#     print(p);
#     ps[[lbl]] = p;
#   }
# }

#--use mgcVis to plot iterms
viz = mgcViz::getViz(mdl);#get mgcViz object based on the model
# #----plot terms
# ps = plot(viz,allTerms = TRUE);
# for (p in ps$plots){
#   print(p);
#   ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;
# }
# 
# #----check model
o <- check.gamViz(viz);
print(o); #--print diagnostic plots
pg = cowplot::plot_grid(plotlist=o,ncol=2);#--need to do this to use ggsave
ggsave(paste0("plot3b_females_DiagnosticPlots.",device),plot=pg,width=8,height=5);

#--plot smooths again
os1 = list();
os2 = list();
for (s in 1:ns){
  smth = smths[[s]];
  if (!inherits(smth,"random.effect")){
    if (length(smth$term)==1){
      #--plot ith term, if a 1-variable smooth and not a random effect
      o = plot( sm(viz, s) ) +
            l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
            l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
            l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic();
      # print(o);
      # ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;
      os1[[smth$label]] = o;
    } else {
      o = plot( sm(viz, s) ) + 
            l_fitRaster() +
            l_fitContour() + 
            l_points(shape=19,size=1,alpha=0.1) +
            l_rug(alpha=0.8)
      # print(o);
      os2[[smth$label]] = o;
      # o3 = plotRGL(sm(viz, 9),se=TRUE,seWithMean=TRUE,unconditional=TRUE);
      # o3 %<>% filledContour3d();
      # aspect3d(1,1,0.1)
    }
  }
}

ncol = 1*(length(os1)<4)+2*(4<=length(os1));
pg = gridPrint(grobs=os1,ncol=ncol);
ggsave(paste0("smooths_females_1d.",device),pg,units="in",width=6.5,height=5);

ncol = 1*(length(os2)<4)+2*(4<=length(os2));
pg = gridPrint(grobs=os2,ncol=ncol);
ggsave(paste0("smooths_females_2d.",device),pg,units="in",width=6.5,height=8);

#--plot predictions for R from observations on arithmetic scale
int = mdl$coefficients[1];
seFactor=2;
#----calculate estimated values of R for observed data
prd = predict(mdl,dfrDat,se.fit=TRUE,unconditional=TRUE);
dfrPrd = dfrDat %>% dplyr::mutate(prdLnR=prd$fit,  #--offset is lnq, and is 0 identically
                                  seLnR=prd$se.fit,#--standard error
                                  prdR=exp(prdLnR),#--predicted selectivity ratio
                                  lower=exp(prdLnR-seFactor*seLnR),
                                  upper=exp(prdLnR+seFactor*seLnR));
#----calculate selectivity ratio as a function of size only
dfrPrdR = cbind(newDFR("z",grds),h=NA);
prdR    = predict(mdl,newdata=dfrPrdR,type="terms",terms="ti(z)",se.fit=TRUE,unconditional=TRUE,newdata.guaranteed=TRUE);
dfrPrdR %<>% dplyr::mutate(prdLnR=int+prdR$fit,  #--offset is lnq, and is 0 identically
                           seLnR=prdR$se.fit,#--standard error
                           prdR=exp(prdLnR),#--predicted selectivity ratio
                           lower=exp(prdLnR-seFactor*seLnR),
                           upper=exp(prdLnR+seFactor*seLnR));
p = ggplot(dfrPrd,mapping=aes(x=z,y=prdR)) +
      geom_line(mapping=aes(group=h),alpha=0.5) +
      geom_point(mapping=aes(size=n),alpha=0.5) + scale_size_area() +
      geom_ribbon(data=dfrPrdR,mapping=aes(x=z,ymin=lower,ymax=upper),fill="red",alpha=0.4,inherit.aes=FALSE) +
      geom_line(data=dfrPrdR,mapping=aes(x=z,y=prdR),colour="red",inherit.aes=FALSE) +
      labs(x="size (mm CW)",y="estimated selectvity ratio",size="n");
print(p);
ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;

#----calculate sample-weighted average R
dfrAvg = dfrPrd %>% 
           select(y,x,h,z,n,prdR) %>%
           arrange(y,x,z,h) %>%
           group_by(y,x,z) %>%
           mutate(avgR=sum(n*prdR)/sum(n)) %>%
           mutate(seR=sqrt(sum(n*(prdR-avgR)^2)/sum(n))) %>% 
           mutate(np=sum(n)) %>% 
           ungroup() %>%
           select(y,x,z,n=np,avgR,seR) %>%
           distinct() %>%
           arrange(y,x,z) %>%
           mutate(lower=avgR-seFactor*seR,upper=avgR+seFactor*seR);
dfrAvg$y=factor(dfrAvg$y);
lvls = levels(dfrAvg$y);
ny = length(lvls);
#for (iy in 1:ny){
  iy=ny;
  tmp = dfrAvg %>% subset(y %in% lvls[1:iy]);
  p = ggplot(tmp,mapping=aes(x=z,y=avgR,ymin=lower,ymax=upper,colour=y,fill=y)) +
        ggdist::geom_lineribbon(alpha=0.4,step=FALSE) +
        geom_line() + geom_point(mapping=aes(size=n),alpha=0.5) + scale_size_area() +
        ggdist::geom_lineribbon(data=dfrPrdR,mapping=aes(x=z,y=prdR,ymin=lower,ymax=upper),
                                colour="red",fill="red",inherit.aes=FALSE,alpha=0.3) +
        facet_wrap(vars(y),ncol=2) +
        labs(x="size (mm CW)",y="estimated selectvity ratio",size="n",colour="year",fill="year");
  print(p);
  ggsave(plotFN(pltctr),width=9,height=9); pltctr %<>% +1;
#}


