#--functions to calculate log-likelihoods
#' 
#' @title Calculate the log-likelihood for binomially-distributed observations
#' @description Function to calculate the log-likelihood for observations based on 
#' predictions from a model assuming a binomial distribution.
#' @param p_s - vector of observed probability values
#' @param n_s - vector of associated number of trials
#' @param mu_s - predicted values (must be on response scale)
#' @param family - the binomial family object used in the fitted model
#' @param offsets - vector of logit-scale offsets (or NULL, the default)
#' @return a vector of the associated log-likelihood value for each observation
#' @details If the `offsets` vector is given,`family` must also be given. The
#' predicted values are then transformed to the link scale, the offsets are added, 
#' and the results are inverse-tranformed to the response scale.
#' 
#' @examplesIf FALSE
#' #--BINOMIAL regression  models for lnR----
#' famB = stats::binomial(link="logit");
#' ks=c(20,10);
#' k1 = ks[1]; k2 = ks[2];
#' frmla  = p~ti(z,bs="ts",k=k1)   +
#'            ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2);
#' fitp  = mgcv::gam(frmla,family=famB,data=dfrDatp,method="ML",fit=TRUE,offset=lnq,weights=n);
#' prdp = unname(predict(fitp,dfrDatp,"response"));
#' calcLogLike.binomial(dfrDatp$p,dfrDatp$n,prdp,famB,offsets=dfrDatp$lnq);
#' logLik(fitp); 
#' 
#' @export
#' 
calcLogLike.binomial<-function(p_s,n_s,mu_s,family=NULL,offsets=NULL){
  if (!is.null(offsets)){
    if (is.null(family)||(family$family!="binomial")){
      #--transform mu_s to response scale
      stop("calcLogLike.binomial: offsets are specified, so `family` must be the stats::binomial object used in the fitted model.")
    }
    #--mu_s are predicted values on response scale; offsets are on link scale
    mu_s = family$linkfun(mu_s);#--transform to link scale
    mu_s = mu_s + offsets;      #--add in offsets
    mu_s = family$linkinv(mu_s);#--transform to response scale
  }
    
  # if (!is.null(family)){
  #   #--transform mu_s to response scale
  #   if (family$family!="binomial")
  #     stop("calcLogLike.binomial: offsets are specified, so `fam` must be the stats::binomial object used in the fitted model.")
  #   mu_s = family$linkinv(mu_s)
  # }
  dbinom(p_s*n_s,n_s,mu_s,log=TRUE)
}

#--test BINOMIAL likelihood----
if (FALSE){
famB = stats::binomial(link="logit");
ks=c(20,10);
k1 = ks[1]; k2 = ks[2];
frmla  = p~ti(z,bs="ts",k=k1)   +
           ti(d,bs="ts",k=k2)   + ti(t,bs="ts",k=k2);
fitp  = mgcv::gam(frmla,family=famB,data=dfrDatp,method="ML",fit=TRUE,offset=lnq,weights=n);
prdp = unname(predict(fitp,dfrDatp,"response"));
sum(calcLogLike.binomial(dfrDatp$p,dfrDatp$n,prdp,famB,offsets=dfrDatp$lnq));
logLik(fitp);
}