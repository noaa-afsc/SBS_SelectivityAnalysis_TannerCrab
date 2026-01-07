#--ExtraPoissonVariability: Tweedie distribution
##--variance-mean relationship based on 'reduced' model fitting
##--log(V) ~ intercept + slope * log(mean) + fleet x sex
pwr = 1.5456; #--'lnmn' from reduced model
phi = 4.6743; #--intercept from reduced model

aswB = 0.002;#--standard area swept (sq nmi) BSFRF
aswN = 0.012;#--standard area swept (sq nmi) NMFS

cpues = c(100,500,1000,5000,10000);#--numbers/sq nmi
fcpus = ordered(cpues);
x = c(seq(5,1000,5),seq(2000,100000,1000));#--CPUE values at which to evaluate probabilities

lst   = list();
for (i in seq_along(cpues)){
  cpue = cpues[i];
  mu=cpue*aswB; #--expected number caught in BSFRF haul
  dt = tweedie::dtweedie(x*aswB,mu=cpue*aswB,phi=exp(phi+(2-pwr)*log(aswB)),power=pwr)
  lst[[fcpus[i]]] = tibble::tibble(mu=mu,cpue=fcpus[i],x=x*aswB,val=dt,type="tweedie");
}
dfrBSFRF = dplyr::bind_rows(lst) |> dplyr::mutate(fleet="BSFRF");
rm(lst);
  
lst   = list();
for (i in seq_along(cpues)){
  cpue = cpues[i];
  mu=cpue*aswN; #--expected number caught in BSFRF haul
  dt = tweedie::dtweedie(x*aswN,mu=cpue*aswN,phi=exp(phi+(2-pwr)*log(aswN)),power=pwr)
  lst[[fcpus[i]]] = tibble::tibble(mu=mu,cpue=fcpus[i],x=x*aswN,val=dt,type="tweedie");
}
dfrNMFS = dplyr::bind_rows(lst) |> dplyr::mutate(fleet="NMFS");
rm(lst);
  
dfr = dplyr::bind_rows(dfrBSFRF,dfrNMFS);

ggplot(dfr |> dplyr::filter(val>0.0001),aes(x=x,y=val,colour=fleet)) + 
  geom_line() + 
  geom_vline(aes(xintercept=mu,colour=fleet),linetype=3) + 
  geom_hline(yintercept=0,linetype=3) + 
  scale_y_continuous(limits=c(0,NA)) + 
  labs(x="catch numbers",y="probability") + 
  facet_grid(cpue~.) + 
  wtsPlots::getStdTheme();


##--simulate ratios based on tweedie
lst   = list();
for (i in seq_along(cpues)){
  cpue = cpues[i];
  rnB = tweedie::rtweedie(1000,mu=cpue*aswB,phi=exp(phi+(2-pwr)*log(aswB)),power=pwr)
  rnN = tweedie::rtweedie(1000,mu=cpue*aswN,phi=exp(phi+(2-pwr)*log(aswN)),power=pwr)
  lst[[fcpus[i]]] = tibble::tibble(cpue=fcpus[i],rnB=rnB,rnN=rnN,r=rnN/mn/(rn),type="tweedie");
}
dfrBSFRF = dplyr::bind_rows(lst) |> dplyr::mutate(fleet="BSFRF");
rm(lst);

dfr = dplyr::bind_rows(dfrBSFRF,dfrNMFS);
