#--apply selectivity model to NMFS survey data for single sex
require(dplyr);
require(ggplot2);
require(magrittr);
require(tcsamSurveyData);
#define sex
x = "male";

#--define grid for sizes
grd_z = seq(27.5,182.5,5);
dfrZ = tibble::tibble(z=grd_z);

#--get survey haul data with spatial covariates depth, temperature, phi, and sorting
#--(year,hauljoin,depth,temp,phi,sorting)
dfrHDwSCs = wtsUtilities::getObj("dfrHD_NMFS_WithInterpolatedSedValues.RData")

#--expand table with size grid info
dfrHDwSCs %<>% sf::st_drop_geometry() %>%
               dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting) %>%
               tidyr::expand_grid(dfrZ);

#--get model
mdl_res = wtsUtilities::getObj(paste0("selModels.Binom.Best.",tolower(x),"s.RData"));
mdl = mdl_res$model;

#--predict haul-specific selectivity ratio (year,hauljoin,size,R)
prd = predict(mdl,newdata=dfrHDwSCs,se.fit=TRUE,type="link");#--return predictions on logit-scale w/ no offset
seR = sqrt(exp(prd$se.fit^2)-1)*exp(prd$fit+(prd$se.fit^2)/2);
seFactor = 2;
dfrRbyH = dfrHDwSCs %>% 
            dplyr::mutate(prdLnR=prd$fit,  #--offset lnq is 0, identically
                          seLnR=prd$se.fit,#--standard error
                          prdR=exp(prdLnR),#--arithmetic-scale predicted selectivity ratio
                          seR=seR,         #--arithmetic-scale standard error
                          lower=exp(prdLnR-seFactor*seLnR),
                          upper=exp(prdLnR+seFactor*seLnR));

#--look at results by haul
wgtBySE<-function(x){
  if (wtsUtilities::Sum(x)>0) {
      wgt=x/wtsUtilities::Sum(x);
    } else {
      wgt = 1.0/length(x);
    }
  return(wgt)
}
tmp1 = dfrRbyH %>% subset(dplyr::between(y,2015,2019)&(prdR<10));
tmp2 = tmp1 %>% subset(!(is.na(t)|is.na(f)|is.na(s))) %>%
                dplyr::arrange(y,z,h) %>%
                dplyr::group_by(y,z) %>% 
                dplyr::summarize(mnR=mean(prdR,na.rm=TRUE),
                                 seR=sd(prdR,na.rm=TRUE),
                                 lower=mnR-seFactor*seR,
                                 upper=mnR+seFactor*seR) %>%
                dplyr::ungroup() %>%
                dplyr::select(y,z,mnR,seR,lower,upper);
                      
p = ggplot2::ggplot() +
      # ggplot2::geom_point(data=tmp1,mapping=aes(x=z,y=prdR,colour=as.factor(y),group=h),alpha=0.3) +
      # ggplot2::geom_line(data=tmp1,mapping=aes(x=z,y=prdR,colour=as.factor(y),group=h),alpha=0.3) +
      #ggplot2::geom_ribbon(data=tmp2,mapping=aes(x=z,ymin=lower,ymax=upper),colour="grey",alpha=0.3) +
      ggplot2::geom_line(data=tmp2,mapping=aes(x=z,y=mnR),colour="black",size=1) +
      #ggplot2::lims(y=c(0,2)) + 
      labs(y="predicted selectivity (ratio)",colour="year") +
      ggplot2::facet_grid(y~.);
print(p);   

#--read in CPUE (abundance) by haul/station for selectivity-by-haul weighting
dfrCPUE = wtsUtilities::getObj(paste0("dfrCPUE.",x,"s.RData"));

#--calculate year-specific survey selectivity by averaging over haul-specific selectivity
wgt<-function(x){
  if (wtsUtilities::Sum(x)>0) {
      wgt=x/wtsUtilities::Sum(x);
    } else {
      wgt = 1.0/length(x);
    }
  return(wgt)
}
dfrRbyY = dfrRbyH %>%  
            subset((prdR<10)) %>% 
            dplyr::left_join(dfrCPUE,by=c("y","h","z")) %>%
            subset(!(is.na(t)|is.na(f)|is.na(s))) %>%
            dplyr::arrange(y,z,h) %>%
            dplyr::group_by(y,z) %>% 
            dplyr::mutate(w = wgt(numCPUE)) %>%
            dplyr::summarize(mnR=wtsUtilities::Sum(w*prdR),
                             seR=sqrt(wtsUtilities::Sum(w*(prdR-mnR)^2)),
                             lower=quantile(prdR,probs=0.05,na.rm=TRUE),
                             upper=quantile(prdR,probs=0.95,na.rm=TRUE)) %>%
            dplyr::ungroup() %>%
            dplyr::select(y,z,mnR,seR,lower,upper);

#--plot year-specific survey selectivity by multi-year intervals
intv = 10;#--number of intervals per plot
tmp = dfrRbyY %>% dplyr::mutate(group=floor((y-min(y))/intv+1),lbl=paste0((group-1)*intv+min(y)," : ",group*intv+min(y)-1));
p = ggplot2::ggplot() +
      ggplot2::geom_ribbon(data=tmp,mapping=aes(x=z,y=mnR,ymin=lower,ymax=upper,group=as.factor(y)),alpha=0.3) +
      ggplot2::geom_line(data=tmp,mapping=aes(x=z,y=mnR,colour=as.factor(y),group=as.factor(y))) + 
      ggplot2::facet_grid(lbl~.) +
      ggplot2::labs(x="size (mm CW)",y="survey-averaged selectivity (ratio)",colour="year");
print(p);   

  