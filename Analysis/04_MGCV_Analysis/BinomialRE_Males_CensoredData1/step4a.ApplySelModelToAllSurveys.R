#--apply selectivity model to NMFS survey data for single sex

#--load required packages
require(dplyr);
require(ggplot2);
require(magrittr);
require(mgcv);
#--other required packages
#----sf
#----tcsamSurveyData
#----tibble
#----tidyr

#--select sex
x = "male"; maxZ = 180;

#--get model
mdl = wtsUtilities::getObj(paste0("selModels.BinomRE.Best.",tolower(x),"s.RData"));

#--set up plotting output
pltctr = 1;
plotFN<-function(n){paste0("plots4a_",wtsUtilities::formatZeros(n),".pdf")}

#--determine valid interpolation ranges of covariates for predictions
dfrDat = mdl$model;
rng_d = range(dfrDat$d,na.rm=TRUE);
rng_t = range(dfrDat$t,na.rm=TRUE);
rng_f = range(dfrDat$f,na.rm=TRUE);
rng_s = range(dfrDat$s,na.rm=TRUE);
dfrRngs = tibble::tibble(d=rng_d,t=rng_t,f=rng_f,s=rng_s);

#--define grid for sizes
grd_z = seq(27.5,maxZ+2.5,5);
dfrZ = tibble::tibble(z=grd_z);

#--get survey haul data with spatial covariates depth, temperature, phi, and sorting
#--Limit to surveys with 83-112 net (1982+)
dfrHDwSCs = wtsUtilities::getObj("../dfrHD_NMFS_WithInterpolatedSedValues.RData") %>%
              dplyr::filter(YEAR>=1982) %>%
              sf::st_drop_geometry() %>%
              dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting);
#--check coverage of covariates relative to model fit
nHaulsByYr = dfrHDwSCs %>% group_by(y) %>% summarize(numHauls=n()) %>%
                          ungroup() %>% mutate(type="original");
p = ggplot(dfrHDwSCs,mapping=aes(x=d,y=t)) + geom_bin2d() +
      geom_vline(data=dfrRngs,mapping=aes(xintercept=d)) +
      geom_hline(data=dfrRngs,mapping=aes(yintercept=t)) +
      scale_fill_viridis_c(option="plasma") +
      labs(x="haul depth (m)",y="temperature (deg C)",fill="haul\ncount");
print(p);
ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;
p = ggplot(dfrHDwSCs,mapping=aes(x=f,y=s)) + geom_bin2d() +
      geom_vline(data=dfrRngs,mapping=aes(xintercept=f)) +
      geom_hline(data=dfrRngs,mapping=aes(yintercept=s)) +
      scale_fill_viridis_c(option="plasma") +
      labs(x="phi (ln-scale)",y="sorting",fill="haul\ncount");
print(p);
ggsave(plotFN(pltctr),width=8,height=5); pltctr %<>% +1;

#--filter table to restrict covariates to ranges sampled in SBS studies
filterCovariates = FALSE;
if (filterCovariates){
  dfrHDwSCs %<>%  dplyr::filter(between(d,rng_d[1],rng_d[2]),between(t,rng_t[1],rng_t[2]),
                                between(f,rng_f[1],rng_f[2]),between(s,rng_s[1],rng_s[2]));
  nHaulsByYrFiltered = dfrHDwSCs %>% 
                         group_by(y) %>% 
                         summarize(numHauls=n()) %>%
                         ungroup() %>% mutate(type="filtered");
  nHaulsByYr = rbind(nHaulsByYr,nHaulsByYrFiltered) %>%
                 tidyr::pivot_wider(id_cols=y,names_from=type,values_from=numHauls);
}

#--expand table with size grid info
dfrHDwSCs %<>% tidyr::expand_grid(dfrZ);

#--predict haul-specific selectivity ratio (year,hauljoin,size,R)
prd = predict(mdl,newdata=dfrHDwSCs,se.fit=TRUE,type="link");#--return predictions on logit-scale w/ no offset
seR = sqrt(exp(prd$se.fit^2)-1)*exp(prd$fit+(prd$se.fit^2)/2);
seFactor = 2;
dfrRbyH = dfrHDwSCs %>% 
            dplyr::mutate(prdLnR=prd$fit,  #--offset lnq is 0 identically
                          seLnR=prd$se.fit,#--standard error
                          prdR=exp(prdLnR),#--arithmetic-scale predicted selectivity ratio
                          seR=seR,         #--arithmetic-scale standard error
                          lower=exp(prdLnR-seFactor*seLnR),
                          upper=exp(prdLnR+seFactor*seLnR));
wtsUtilities::saveObj(dfrRbyH,"dfrRbyH.RData");

#--read in CPUE (abundance) by haul/station for selectivity-by-haul weighting
dfrCPUE = wtsUtilities::getObj(paste0("../dfrCPUE.",x,"s.RData"));

#--calculate year-specific survey selectivity by averaging over haul-specific selectivity
#---function to weights by n and inverse variance for averaging
calcWgts<-function(n,se){
  var = se^2;
  if (wtsUtilities::Sum(n)>0) {
      wgt=(n/var)/wtsUtilities::Sum(n/var);
    } else {
      wgt = (1/var)/wtsUtilities::Sum(1/var);
    }
  return(wgt)
}
dfrRbyY = dfrRbyH %>%  
            dplyr::left_join(dfrCPUE,by=c("y","h","z")) %>%
            subset(!(is.na(t)|is.na(f)|is.na(s))) %>%
            dplyr::arrange(y,z,h) %>%
            dplyr::group_by(y,z) %>% 
            dplyr::mutate(w = calcWgts(numCPUE,seR)) %>%
            dplyr::summarize(numIndivs=sum(numIndivs),
                             wgt=sum(w),
                             mnR=wtsUtilities::Sum(w*prdR),
                             seR=sqrt(wtsUtilities::Sum(w*(prdR-mnR)^2)),
                             lower=mnR-seFactor*seR,
                             upper=mnR+seFactor*seR) %>%
            dplyr::ungroup() %>%
            dplyr::select(y,z,numIndivs,wgt,mnR,seR,lower,upper);
wtsUtilities::saveObj(dfrRbyY,"dfrRbyY.RData");


  