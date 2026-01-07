#--extra poisson variance effect on numbers caught
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

#--simulate catch numbers for various densities (lambdas), NMFS q's, and levels of extra variance
aswB = 0.002;#--standard area swept (sq nmi) BSFRF
aswN = 0.012;#--standard area swept (sq nmi) NMFS

N = 1000; #--number of samples to take
lambdas = c(1000,5000,10000,50000,100000); #--densities (numbers/sq nmi)
qs      = seq.default(0.1,1,0.1);          #--NMFS catchabilities
phis = c(0.1,0.25,0.50,0.75,1.0,1.25,1.5); #--extra dispersion factors
lst0 = list();
for (i in seq_along(lambdas)){
  lambda_ = lambdas[i];
  lst1 = list();
  for (j in seq_along(qs)){
    q_    = qs[j];
    lst2  = list();
    type_ = "poisson";
    dfrp  = tibble::tibble(type=type_,
                           lambda=lambda_,
                           q=q_,
                           phi=0.0,
                           nB=rpois(N,lambda=lambda_*aswB),
                           nN=rpois(N,lambda=lambda_*aswN*q_),
                           rQ=(nN/aswN)/(nB/aswB));
    lst2[[type_]] = dfrp; rm(dfrp);
    for (k in seq_along(phis)){
      phi_ = phis[k];
      type_ = paste0("extra ",phi_);
      dfrp  = tibble::tibble(type=type_,
                             lambda=lambda_,
                             q=q_,
                             phi=phi_,
                             nB=rnbinom(N,mu=lambda_*aswB,size=1/phi_),
                             nN=rnbinom(N,mu=lambda_*aswN*q_,size=1/phi_),
                             rQ=(nN/aswN)/(nB/aswB));
      lst2[[type_]] = dfrp; rm(dfrp);
    }#--k
    lst1[[as.character(q_)]] = dplyr::bind_rows(lst2);
    rm(lst2);
  }#--j
  lst0[[as.character(lambda_)]] = dplyr::bind_rows(lst1);
  rm(lst1);
}#--i
dfr = dplyr::bind_rows(lst0) |> 
        dplyr::mutate(type=factor(type,levels=c("poisson",paste0("extra ",phis))));
rm(lst0);

dfrStats = dfr |> dplyr::group_by(type,lambda,q,phi) |> 
             dplyr::summarize(mnB=mean(nB,na.rm=TRUE),
                              mnN=mean(nN,na.rm=TRUE),
                              nRs=sum(is.finite(rQ)),
                              mnQ=mean(rQ[is.finite(rQ)],na.rm=TRUE),
                              p05Q=quantile(rQ[is.finite(rQ)],probs=0.05,na.rm=TRUE),
                              p95Q=quantile(rQ[is.finite(rQ)],probs=0.95,na.rm=TRUE),
                              p10Q=quantile(rQ[is.finite(rQ)],probs=0.10,na.rm=TRUE),
                              p90Q=quantile(rQ[is.finite(rQ)],probs=0.90,na.rm=TRUE)) |> 
             dplyr::ungroup();

p = ggplot(dfrStats,
            aes(x=q,y=nRs/N)) + 
       geom_line() + 
       geom_hline(yintercept=1,linetype=2) + 
       facet_grid(type~lambda,scales="fixed") + 
       labs(x="true q",y="% successful haul pairs") + 
       wtsPlots::getStdTheme();
print(p);

p0 = ggplot(dfrStats |> dplyr::filter(type=="poisson"),
            aes(x=q,y=mnQ/q,ymin=p05Q/q,ymax=p95Q/q)) + 
       geom_ribbon(alpha=0.5,colour=NA,fill="blue") + 
       geom_ribbon(aes(ymin=p10Q/q,ymax=p90Q/q),alpha=0.25,colour=NA,fill="green") + 
       geom_line() + 
       geom_hline(yintercept=1,linetype=2) + 
       facet_grid(lambda~type) + 
       wtsPlots::getStdTheme();
print(p0);

p1 = ggplot(dfrStats |> dplyr::filter(lambda==5000),
            aes(x=q,y=mnQ/q,ymin=p05Q/q,ymax=p95Q/q)) + 
       geom_ribbon(alpha=0.5,colour=NA,fill="blue") + 
       geom_ribbon(aes(ymin=p10Q/q,ymax=p90Q/q),alpha=0.25,colour=NA,fill="green") + 
       geom_line() + 
       geom_hline(yintercept=1,linetype=2) + 
       facet_grid(type~lambda,scales="free_y") + 
       labs(x="true q",y="ratio: (est q)/(true q)") + 
       wtsPlots::getStdTheme();
print(p1);

p1 = ggplot(dfrStats,
            aes(x=q,y=mnQ/q,ymin=p05Q/q,ymax=p95Q/q)) + 
       geom_ribbon(alpha=0.5,colour=NA,fill="blue") + 
       geom_ribbon(aes(ymin=p10Q/q,ymax=p90Q/q),alpha=0.25,colour=NA,fill="green") + 
       geom_line() + 
       geom_hline(yintercept=1,linetype=2) + 
       facet_grid(type~lambda,scales="free_y") + 
       labs(x="true q",y="ratio: (est q)/(true q)") + 
       wtsPlots::getStdTheme();
print(p1);

for (i in seq_along(lambdas)){
  p1 = ggplot(dfrStats |> dplyr::filter(lambda==lambdas[i]),
              aes(x=q,y=mnQ/q,ymin=p05Q/q,ymax=p95Q/q)) + 
         geom_ribbon(alpha=0.5,colour=NA,fill="blue") + 
         geom_ribbon(aes(ymin=p10Q/q,ymax=p90Q/q),alpha=0.25,colour=NA,fill="green") + 
         geom_line() + 
         geom_hline(yintercept=1,linetype=2) + 
         facet_grid(type~lambda,scales="fixed") + 
         labs(x="true q",y="ratio: (est q)/(true q)") + 
         wtsPlots::getStdTheme();
  print(p1);
}

#--look at "average" % bias as function of density and extra variability
##---by averaging over q
dfrS1 = dfrStats |> dplyr::group_by(type,lambda,phi) |> 
          dplyr::summarize(mnPctBias=100*(mean(mnQ/q)-1)) |> 
          dplyr::ungroup();
p1 = ggplot(dfrS1,
            aes(x=phi,y=mnPctBias,colour=factor(lambda,levels=as.character(lambdas)))) + 
       geom_line() + 
       geom_hline(yintercept=0,linetype=2) + 
       #facet_grid(lambda~.,scales="free_y") + 
       labs(x="extra poisson variability",y="mean % bias",colour="density\n(number/sq. nmi.)") + 
       wtsPlots::getStdTheme() + 
       theme(legend.position="inside",
            legend.position.inside=c(0.01,0.99),
            legend.justification=c(0,1));
print(p1);


