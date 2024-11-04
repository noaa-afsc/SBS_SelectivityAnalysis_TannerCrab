prdMod<-function(mdl,trms,lst,type="response",keep=NULL,p=0.05){
  dfr = wtsMGCV::createGridTbl(lst);
  if (any(trms=="all")){
    #--add intercept and all smooth terms
    trmsp = "(Intercept)";
    trms = gratia::smooths(mdl);
    for (trm in trms) trmsp = c(trmsp,trm);
    trms = trmsp;
  }
  prd = dplyr::bind_cols(
            dfr,
            tibble::as_tibble(
              mgcv::predict.gam(mdl,dfr,type=type,terms=trms,se.fit=TRUE),
            ) |> 
            dplyr::mutate(type="fit",
                          lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                          uci=qnorm(p,fit,se.fit,lower.tail=FALSE),
                          terms=paste(trms,collapse=" + "))
        ) |> 
          dplyr::rename(emp_sel=fit);
  if (!is.null(keep)){
    prd = prd |> 
            dplyr::distinct(pick(tidyselect::any_of(keep),emp_sel,se.fit,lci,uci,terms));
    drop = names(lst)[!(names(lst) %in% keep)];
    for (drp in drop) prd[[drp]] = NA;
  }
  return(prd);
}
   mdl3m = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,y,bs="fs"),
                     family=tw,data=dfrESsp,method="REML",
                     weights=n_BSFRF/mean(n_BSFRF));
  wtsMGCV::gam.check.plots(mdl3m)
  #--on response scale
  prd3mb = prdMod(mdl3m,trms=c("(Intercept)","ti(z)"),lst=grids,type="response",keep="z");
  ggplot(prd3mb,aes(x=z,y=emp_sel,ymin=lci,ymax=uci)) + geom_line() + geom_ribbon(alpha=0.3);
  prd3m = prdMod(mdl3m,trms=NULL,lst=grids,type="response");
  ggplot(prd3m,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrESsp,mapping=aes(x=z,y=emp_sel,colour=y,size=n_BSFRF),inherit.aes=FALSE) +  
    scale_size_area() + 
    scale_y_continuous(limits=c(0,1.5),oob=scales::oob_squish) +
    geom_line(data=prd3mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    geom_ribbon(data=prd3mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 
  #--on link scale
  prd3mb = prdMod(mdl3m,trms=c("(Intercept)","ti(z)"),lst=grids,type="link",keep="z");
  prd3m = prdMod(mdl3m,trms=NULL,lst=grids,type="link");
  ggplot(prd3m,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrESsp,mapping=aes(x=z,y=log(emp_sel),colour=y,size=n_BSFRF),inherit.aes=FALSE) +  
    scale_size_area() + 
    scale_y_continuous(limits=c(-5,5),oob=scales::oob_squish) +
    geom_line(data=prd3mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    geom_ribbon(data=prd3mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 

  dfr = gratia:::smooth_estimates(mdl3m,select="ti(z,y)",unconditional=TRUE,overall_uncertainty=TRUE) |> 
          gratia::add_confint();
  dfrPRs = dplyr::bind_cols(dfrESsp |> dplyr::filter(!is.na(emp_sel)),
                            gratia::partial_residuals(mdl3m,select="ti(z,y)") |> dplyr::rename(residuals="ti(z,y)"));
  ggplot(dfr,aes(x=z,y=.estimate,ymin=.lower_ci,ymax=.upper_ci,colour=y,fill=y)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrPRs,mapping=aes(x=z,y=residuals,colour=y,size=n_BSFRF),inherit.aes=FALSE) +  
    scale_size_area(name="n (BSFRF)") + facet_wrap(~y); # + 
    # scale_y_continuous(limits=c(-5,5),oob=scales::oob_squish) +
    # geom_line(data=prd3mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    # geom_ribbon(data=prd3mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 
  dfrMB = gratia:::smooth_estimates(mdl3m,select="ti(z)",unconditional=TRUE,overall_uncertainty=TRUE) |> 
            gratia::add_confint();
  dfrMBPRs = dplyr::bind_cols(dfrESsp |> dplyr::filter(!is.na(emp_sel)),
                              gratia::partial_residuals(mdl3m,select="ti(z)") |> dplyr::rename(residuals="ti(z)"));
  ggplot(dfrMB,aes(x=z,y=.estimate,ymin=.lower_ci,ymax=.upper_ci)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrMBPRs,mapping=aes(x=z,y=residuals,size=n_BSFRF,colour=y,fill=y),inherit.aes=FALSE) +  
    scale_size_area() + facet_wrap(~y);# +
    # scale_y_continuous(limits=c(-5,5),oob=scales::oob_squish) +
    # geom_line(data=prd3mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    # geom_ribbon(data=prd3mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 

#--prediction intervals