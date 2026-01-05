#--simulate SBS sampling
require(ggplot2);

plotObs<-function(values,exp,xtitle){
  dfr = tibble::tibble(lambda=lambda,sNMFS=sNMFS,sBSFRF=sBSFRF,vExp=exp,value=values)
  dfrStats = dfr |> dplyr::group_by(lambda,sNMFS,sBSFRF,vExp) |> 
                 dplyr::filter(is.finite(value)) |>
                 dplyr::summarize(mn=mean(value,na.rm=TRUE),
                                  lci=quantile(value,probs=0.1,na.rm=TRUE),
                                  uci=quantile(value,probs=0.9,na.rm=TRUE)) |> 
                 dplyr::ungroup();
  
  p = ggplot(dfr,aes(x=value,y=after_stat(ndensity))) + geom_histogram() + 
         geom_vline(xintercept=exp,linetype=3,linewidth=1) + 
         geom_vline(xintercept=dfrStats$lci,linetype=3,linewidth=1,colour="blue") + 
         geom_vline(xintercept=dfrStats$mn, linetype=3,linewidth=1,colour="red") + 
         geom_vline(xintercept=dfrStats$uci,linetype=3,linewidth=1,colour="blue") + 
         labs(x=xtitle) + 
         wtsPlots::getStdTheme() + 
         theme(axis.title.y=element_blank());
  return(list(dfr=dfr,dfrStats=dfrStats,p=p))
}

plotResults<-function(values){
  dfr      = tibble::tibble(lambda=lambda,sNMFS=sNMFS,sBSFRF=sBSFRF,value=values);
  dfrStats = dfr |> dplyr::group_by(lambda,sNMFS,sBSFRF) |> 
                 dplyr::filter(is.finite(value)) |>
                 dplyr::summarize(mn=mean(value,na.rm=TRUE),
                                  lci=quantile(value,probs=0.1,na.rm=TRUE),
                                  uci=quantile(value,probs=0.9,na.rm=TRUE)) |> 
                 dplyr::ungroup();
  
  p = ggplot(dfr,aes(x=value,y=after_stat(ndensity))) + geom_histogram() + 
         geom_vline(xintercept=sExp,linetype=3,linewidth=1) + 
         geom_vline(xintercept=dfrStats$lci,linetype=3,linewidth=1,colour="blue") + 
         geom_vline(xintercept=dfrStats$mn, linetype=3,linewidth=1,colour="red") + 
         geom_vline(xintercept=dfrStats$uci,linetype=3,linewidth=1,colour="blue") + 
         geom_vline(xintercept=1,linetype=3,linewidth=1,colour="green") + 
         scale_x_continuous(limits=c(0,NA)) + 
         labs(x="selectivity ratio") + 
         wtsPlots::getStdTheme() + 
         theme(axis.title.y=element_blank());
  return(list(dfr=dfr,dfrStats=dfrStats,p=p))
}

#--set up study----
n = 5000;#--number of samples to generate

#--set up gear and haul characteristics----
rQs    = 6;   #--ratio of NMFS area swept to BSFRF area swept
sBSFRF = 1.0;#--BSFRF selectivity
sNMFS  = 0.5; #--NMFS selectivity

#--expected values----
sExp = sNMFS/sBSFRF;                #--expected selectivity ratio
rExp = sNMFS*rQs/sBSFRF;            #--expected catch ratio
pExp = sNMFS*rQs/(sNMFS*rQs+sBSFRF);#--expected NMFS proportion of catch

#--loop over crab densities----
lst = list();
lambdas = c(5:19,seq(25,100,5)); #--crab density (number of crab in BSFRF area swept)
for (lambda in lambdas){
  ##----generate catches----
  nBSFRF = rpois(n,sBSFRF*lambda);
  nNMFS  = rpois(n,sNMFS*rQs*lambda);
  
  ##--calculate observed values----
  rObs = nNMFS/nBSFRF;        #--observed catch ratio
  pObs = nNMFS/(nNMFS+nBSFRF);#--observed NMFS proportion of catch
  
  ##--plot catch ratios and proportions NMFS----
  rObsLst = plotObs(rObs,rExp,"catch ratio");
  pObsLst = plotObs(pObs,pExp,"proprtion NMFS");
  
  #--calculate estimated selectivity ratios----
  sObsR = (nNMFS/nBSFRF)*(1/rQs);          #--selectivity ratio based on catch ratio
  sObsP = exp(log(pObs/(1-pObs))-log(rQs));#--selectivity ratio based on NMFS proportion
  
  ##--plot results----
  lbl = paste0("E[nBSFRF] = ",lambda,"\n",
               "sel(NMFS) = ",sNMFS)
  sLstR = plotResults(sObsR);
  sLstP = plotResults(sObsP)
  pg = cowplot::plot_grid(rObsLst$p,
                           pObsLst$p,
                           cowplot::ggdraw(sLstR$p) + 
                             cowplot::draw_plot_label(lbl,x=0.98,y=0.98,hjust=1.02,vjust=1.05,size=10),
                           sLstP$p,
                           ncol=2,byrow=FALSE);
  ##--save to lst----
  lst[[paste(lambda)]] = list(lambda=lambda,sExp=sExp,rExp=rExp,pExp=pExp,
                              rStats=rObsLst$dfrStats,
                              pStats=pObsLst$dfrStats,
                              sStatsR = sLstR$dfrStats,
                              sStatsP = sLstP$dfrStats,
                              plt=pg)
}

lstStats = list();
for (i in 1:length(lst)) {
  lstStats[[i]] = lst[[i]]$sStatsR;
  #print(lst[[i]]$plt);
}
dfrStats = dplyr::bind_rows(lstStats);
p = ggplot(dfrStats,aes(x=lambda,y=mn,ymin=lci,ymax=uci)) + 
      geom_ribbon(alpha=0.2) + geom_line() + geom_point() + 
      geom_hline(yintercept=sNMFS,linetype=3,linewidth=1) + 
      scale_y_continuous(limits=c(0,NA)) + 
      labs(x="Expected BSFRF Catch (numbers)",y="NMFS Selectivity") + 
      wtsPlots::getStdTheme();

wtsUtilities::saveObj(list(lst=lst,dfrStats=dfrStats,plt=p),"EstimatedOneHaulSamplingVariability.RData");


