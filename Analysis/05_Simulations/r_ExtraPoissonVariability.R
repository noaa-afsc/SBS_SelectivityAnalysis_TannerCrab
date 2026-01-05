#--extra poisson variance
require(ggplot2)

##   lambda = density
##   A      = area swept
##   n      = number caught (random variable)

##--poisson-distributed variable
###   n ~ P(lambda * A)
###   E[n] = mu = lambda * A
###   V[n] = mu = lambda * A

##--extra poisson variability
### if theta ~ Gamma(mu,1/phi)  (Gamma parameterized by mean mu and dispersion factor phi) 
### E[theta] = mu           (= a * b; a = shape parameter, b = scale parameter)
### V[theta] = mu^2*phi     (= a * b^2 = (ab)*b = mu * b = mu^2*phi => b = mu * phi )
### then E[n] = mu
###      V[n] = mu + (mu^2)*phi
###NOTE: rnbinom and other nbinom R functions can use mu and "size" = 1/phi
###      If Gamma parameterized in terms of a and b (shape and scale parameters),
###      then 
###        a = mu / b; 
###        mu^2 * phi = a*b^2 => 
###        mu^2 * phi = (mu/b) * b^2 => 
###        mu *phi = b =>
###        a = 1/phi     = size     (per rnbinom notation)
###        b = mu * phi  = mu/size  (per rnbinom notation)

lambda = 50; A = 1; mu = lambda *A;
x = (0:100);
dp = dpois(x=x,lambda=mu);
dfrp = tibble::tibble(x,val=dp,type="poisson",phi=0);
phis = c(0.1,0.25,0.50,0.75,1.0,1.25,1.5);
lst = list();
for (phi in phis){
  db = dnbinom(x=x,mu=mu,size=1/phi);
  lst[[as.character(phi)]] = tibble::tibble(x,val=db,type=paste0("extra ",phi),phi=phi);
}
dfrb = dplyr::bind_rows(lst);
dfr = dplyr::bind_rows(dfrp,dfrb);
dfr$type = factor(dfr$type,levels=c("poisson",paste0("extra ",phis)));
ggplot(dfr,aes(x=x,y=val,colour=type)) + 
  geom_point() + geom_line() + 
  geom_vline(xintercept=mu,linetype=3) + 
  #geom_vline(xintercept=mu+sqrt(mu)*c(-1,1),linetype=3) + 
  geom_vline(xintercept=mu,linetype=3) + 
  #geom_vline(xintercept=mu+sqrt(mu*(1+mu*phi))*c(-1,1),linetype=3,colour="orange") + 
  labs(x="number caught",y="probability") + 
  wtsPlots::getStdTheme()




