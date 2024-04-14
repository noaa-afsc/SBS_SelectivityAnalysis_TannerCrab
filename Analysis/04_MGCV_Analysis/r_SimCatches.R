#--SBS simulations

#--get project path
dirPrj = rstudioapi::getActiveProject();

#--get SBS haul data
lstSBS = wtsUtilities::getObj(file.path(dirPrj,"Analysis/04_MGCV_Analysis/rda_SBS_DataWithSedData.RData"))

#--calculate mean areas swept
A_n = mean((lstSBS$dfrHD_NMFS_SBS |> dplyr::distinct(HAULJOIN,AREA_SWEPT_VARIABLE))$AREA_SWEPT_VARIABLE);
A_b = mean((lstSBS$dfrHD_BSFRF_SBS|> dplyr::distinct(HAULJOIN,AREA_SWEPT_VARIABLE))$AREA_SWEPT_VARIABLE);

#--calculate mean CPUE by size
dfrCPUE = tcsamSurveyData::calcCPUE.ByHaul(lstSBS$dfrHD_BSFRF_SBS,
                                           lstSBS$dfrID_BSFRF_SBS,
                                           bySex=TRUE,
                                           bySize=TRUE,
                                           byMaturity=FALSE,
                                           byShellCondition=FALSE,
                                           cutpts=seq(0,200,5)) |> 
            dplyr::filter(SEX!="MISSING")
dfrMnCPUE = dfrCPUE |> dplyr::select(sex=SEX,size=SIZE,numCPUE) |>
              dplyr::group_by(sex,size) |> 
              dplyr::summarize(mean=mean(numCPUE,na.rm=TRUE),
                               var=var(numCPUE,na.rm=TRUE),
                               inv_dsp=(var-mean)/(mean^2),
                               lci=quantile(numCPUE,probs=0.10,na.rm=TRUE,names=FALSE),
                               uci=quantile(numCPUE,probs=0.90,na.rm=TRUE,names=FALSE)) |> 
              dplyr::ungroup();
ggplot(dfrMnCPUE,aes(x=size,y=A_b*mean,ymin=A_b*lci,ymax=A_b*uci,colour=sex,fill=sex)) + 
  geom_ribbon(alpha=0.5) + geom_line() + 
  labs(x="size (mm CW)",y="mean BSFRF catch") + 
  wtsPlots::getStdTheme();
ggplot(dfrMnCPUE,aes(x=size,y=inv_dsp,colour=sex)) + geom_line() + 
  scale_y_log10() + 
  labs(x="size (mm CW)",y="1/(dispersion parameter)") + wtsPlots::getStdTheme();

#--For crab with CW "z" at location of SBS haul "h"
#----D_hz : local density 
#----N_bhz : actual BSFRF catch
#----N_nhz : actual NMFS catch
#----E<N_bhz> = A_b*D_hz      : expected catch in BSFRF haul
#----E<N_nhz> = r_hz*A_n*D_hz : expected catch in NMFS haul
#--rearranging terms gives
#----E<N_nhz> = r_hz * (A_n/A_b) * E<N_bhz>
#--and
#----r_hz = (A_b/A_n)*E<N_nhz>/E<N_bhz>
#----so ideally one would repeat paired hauls in same location to estimate E<.>'s
#----in order to estimate r_hz

#--to simulate catches, need 
#----1. probability distribution for catches: N ~ P(m,V(m)), where m = E<N>
#----2. mean-variance relationship

assemble<-function(r_hz,q_h,N_bhz,N_nhz,plot=FALSE){
  dfr = tibble::tibble(BSFRF=N_bhz,NMFS=N_nhz) |> 
          dplyr::mutate(p=N_nhz/(N_nhz+N_bhz),
                        lgtp = log(p/(1-p)),
                        ntot=NMFS+BSFRF,
                        obsR=(1/q_h)*NMFS/BSFRF,
                        lnq=log(q_h)
                        );
  p1 = wtsPlots::ggMarginal_Hist2D(dfr,BSFRF,NMFS,binwidths=1);
  p2 = ggplot(dfr,aes(x=p)) + geom_histogram()
  p3 = ggplot(dfr,aes(x=obsR)) + geom_histogram() + 
        geom_vline(xintercept=mean(dfr$obsR[dfr$BSFRF>0])) + 
        geom_vline(xintercept=r_hz,colour="red");
  p4 = ggplot(dfr) + 
        geom_histogram(aes(x=lgtp-log(q_h)),fill="lightblue") + 
        geom_histogram(aes(x=log(obsR)),colour="darkgreen",fill=NA) + 
        geom_vline(xintercept=mean(dfr$lgtp[is.finite(dfr$lgtp)])-log(q_h),      linewidth=2,colour="darkblue") + 
        geom_vline(xintercept=mean(log(dfr$obsR)[is.finite(log(dfr$obsR))],na.rm=TRUE),linewidth=2,colour="darkgreen",linetype=3) + 
        geom_vline(xintercept=log(r_hz),colour="red");
  if (plot) {print(p1); print(p2); print(p3); print(p4);}
  lst = list(dfr=dfr,p=p4,p1=p1,p2=p2,p3=p3,
             r=r_hz,
             mnLnR=mean(log(dfr$obsR)[is.finite(log(dfr$obsR))],na.rm=TRUE),
             mnLgtP=mean(dfr$lgtp[is.finite(dfr$lgtp)])-log(q_h));
  return(lst);
}

m_bhz = 10.5;
r_hz  = 0.5;
lnr_hz = log(r_hz);
#--q_h   = (A_n/A_b);
q_h = 1;
m_nhz = r_hz * q_h * m_bhz;
#--simulate a number of hauls
nHs = 10000;

require(mgcv);

#--POISSON distribution (spatially random: V(m) = m)
#----n = rpois(1,m) generates 1 RV with mean m, variance m
lstPois = assemble(r_hz,q_h,plot=TRUE,
                   rpois(nHs,m_bhz),
                   rpois(nHs,m_nhz));
dfrB = lstPois$dfr;
mnLgtP = mean(dfrB$lgtp[is.finite(dfrB$lgtp)]);
estLnR = mnLgtP - log(q_h);
mdlBin1a = mgcv::gam( formula=cbind(NMFS,BSFRF)~1,family=stats::binomial(link="logit"),data=dfrB);
mdlBin1b = stats::glm(formula=cbind(NMFS,BSFRF)~1,family=stats::binomial(link="logit"),data=dfrB);
mdlBin2a = mgcv::gam( formula=cbind(NMFS,BSFRF)~1,family=stats::binomial(link="logit"),data=dfrB,offset=dfrB$lnq);
mdlBin2b = stats::glm(formula=cbind(NMFS,BSFRF)~1,family=stats::binomial(link="logit"),data=dfrB,offset=dfrB$lnq);
mdlBin3a = mgcv::gam( formula=p~1,family=stats::binomial(link="logit"),data=dfrB,offset=dfrB$lnq,weights=ntot);
mdlBin3b = stats::glm(formula=p~1,family=stats::binomial(link="logit"),data=dfrB,offset=dfrB$lnq,weights=ntot);
#--coefficients are on the link scale, 
#----fitted values are on the response scale (no offsets included via offset=...)
intBin1a = coef(mdlBin1a); fitBin1a = fitted(mdlBin1a)[1];
intBin1b = coef(mdlBin1b); fitBin1b = fitted(mdlBin1b)[1];
intBin2a = coef(mdlBin2a); fitBin2a = fitted(mdlBin2a)[1];
intBin2b = coef(mdlBin2b); fitBin2b = fitted(mdlBin2b)[1];
intBin3a = coef(mdlBin3a); fitBin3a = fitted(mdlBin3a)[1];
intBin3b = coef(mdlBin3b); fitBin3b = fitted(mdlBin3b)[1];

#--NEGBINOMIAL distribution
#----n = rnbinom(1,size,prob,mu)
#------size: the dispersion parameter (gamma mixing distribution shape parameter) > 0
#------prob: probability of success in each trial (=size/(size+mu))
#------mu: the mean (alternative parameterization)
#--variance is mu + mu^2/size when parameterized by mu
lstNB = assemble(r_hz,q_h,plot=TRUE,
                 rnbinom(nHs,mu=m_bhz,size=0.1),
                 rnbinom(nHs,mu=m_nhz,size=0.1));
dfrNB = lstNB$dfr |> dplyr::filter(is.finite(log(obsR))); 
mnLnR = mean(log(dfrNB$obsR));
estR = exp(mnLnR);
mdlNB1a = mgcv::gam( formula=obsR~1,family=mgcv::nb(link="log"),data=dfrNB);
intNB1a = coef(mdlNB1a); fitNB1a = fitted(mdlNB1a)[1];

#--BETA regression ()
dfrBR = lstNB$dfr |> dplyr::filter(!is.nan(p));
mnLgtP = mean(dfrBR$lgtp[is.finite(dfrBR$lgtp)]);
mdlBR1a = mgcv::gam(formula=p~1,family=mgcv::betar(link="logit"),data=dfrBR);
mdlBR1b = mgcv::gam(formula=p~1,family=mgcv::betar(link="logit"),data=dfrBR,offset=dfrBR$lnq);
intBR1a = coef(mdlBR1a); fitBR1a = fitted(mdlBR1a)[1];
intBR1b = coef(mdlBR1b); fitBR1b = fitted(mdlBR1b)[1];

dfrBRG = dfrBR |> dplyr::filter(0<p,p<1)
mnP = mean(dfrBRG$p)
mnLgtP = mean(dfrBRG$lgtp);
mnP1 = exp(mean(dfrBRG$lgtp))/(1+exp(mean(dfrBRG$lgtp)))
mdlBRG = betareg::betareg(formula=p~1,link="logit",data=dfrBRG,offset=dfrBRG$lnq)
summary(mdlBRG)
coef(mdlBRG)
predict(mdlBRG)[1] 
 