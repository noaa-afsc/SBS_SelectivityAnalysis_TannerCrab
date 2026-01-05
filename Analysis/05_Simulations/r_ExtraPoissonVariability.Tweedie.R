#--ExtraPoissonVariability: Tweedie distribution
##--variance-mean relationship based on 'reduced' model fitting
##--log(V) ~ intercept + slope * log(mean) + fleet x sex
pwr = 1.5456; #--'lnmn' from reduced model
phi = 4.6743; #--intercept from reduced model

mnAswB = 0.002;#--sq nmi BSFRF
mnAswN = 0.012;#--sq nmi NMFS

cpues = c(100,500,1000);#--number/sq nmi

cpue = 1000;
mu=cpue*mnAswB;
x = c(seq(5,1000,5),seq(2000,100000,1000));
dt = tweedie::dtweedie(x*mnAswB,mu=cpue*mnAswB,phi=exp(phi+(2-pwr)*log(mnAswB)),power=pwr)
dfr = tibble::tibble(x=x*mnAswB,val=dt,type="tweedie")

ggplot(dfr,aes(x=x,y=val,colour=type)) + 
  geom_line() + 
  geom_vline(xintercept=mu,linetype=3) + 
  geom_hline(yintercept=0,linetype=3) + 
  scale_y_continuous(limits=c(0,NA)) + 
  labs(x="catch numbers",y="probability") + 
  wtsPlots::getStdTheme();
