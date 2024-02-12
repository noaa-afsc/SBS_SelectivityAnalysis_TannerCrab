#--subset data before fitting models

#--load required packages
require(dplyr);
require(ggplot2);
require(magrittr);
#--other required packages
#----
#--source required files
source("r_subsetResponseData.R");

#--read in previously-calculated SBS proportions data
dfrProps_AYs =  wtsUtilities::getObj("dfrPropsAllYears.RData");

#--the logit-scale observed proportions "o" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_a is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

#--set up plotting output
pltctr = 1;
plotFN<-function(n){paste0("plots3a_",wtsUtilities::formatZeros(n),".pdf")}

#--extract proportions data with environmental covariates
dfrDat = dfrProps_AYs %>% 
         dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting,
                          z=SIZE,x=SEX,p=propNMFS,n=numTot,q=q,lnq=log(q),lgtp=log(p/(1-p))-lnq);

#--check for outliers in q (should be ~6 with 5 min. BSFRF tow, 30 min. NMFS tow)
qs = c(3,6,12);
lst1 = check_qs(dfrDat,qs);
ggsave(plotFN(pltctr),lst1$ps$p1,"pdf",width=8,height=10); pltctr %<>% +1;
ggsave(plotFN(pltctr),lst1$ps$p2,"pdf",width=8,height=10); pltctr %<>% +1;
ggsave(plotFN(pltctr),lst1$ps$p3,"pdf",width=8,height=10); pltctr %<>% +1;
wtsUtilities::saveObj(lst1$dfrDrop,paste0("dfrDrop.RData"));

#--drop cells with < n_min individuals
n_min = 5;
lst2 = selectSizeData(lst1$dfrDat,n_min=n_min);
ggsave(plotFN(pltctr),lst2$ps$p1,"pdf",width=8,height=5); pltctr %<>% +1;
ggsave(plotFN(pltctr),lst2$ps$p2,"pdf",width=8,height=5); pltctr %<>% +1;
dfrDat = lst2$dfrDat;
wtsUtilities::saveObj(dfrDat,"dfrTrimmedData.RData");
rm(dfrProps_AYs,lst1,lst2);




