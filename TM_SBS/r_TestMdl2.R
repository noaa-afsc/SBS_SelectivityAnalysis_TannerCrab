prdMod<-function(mdl,trms,lst,type="response",keep=NULL,p=0.05){
  dfr = wtsMGCV::createGridTbl(lst);
  if (any(trms=="all")){
    #--add intercept and other parametric terms, plus all smooth terms
    trmsp = "(Intercept)";
    trms = c(gratia::parametric_terms(mdl),gratia::smooths(mdl));
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
  mdl2m = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,by=y,bs="sz"),
                     family=tw,data=dfrESsp,method="REML",
                     weights=n_BSFRF/mean(n_BSFRF));
  wtsMGCV::gam.check.plots(mdl2m)
  #--on response scale
  prd2mb = prdMod(mdl2m,trms=c("(Intercept)","ti(z)"),lst=grids,type="response",keep="z");
  ggplot(prd2mb,aes(x=z,y=emp_sel,ymin=lci,ymax=uci)) + geom_line() + geom_ribbon(alpha=0.3);
  prd2m = prdMod(mdl2m,trms=NULL,lst=grids,type="response");
  ggplot(prd2m,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrESsp,mapping=aes(x=z,y=emp_sel,colour=y,size=n_BSFRF),inherit.aes=FALSE) +  
    scale_size_area() + 
    scale_y_continuous(limits=c(0,1.5),oob=scales::oob_squish) +
    geom_line(data=prd2mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    geom_ribbon(data=prd2mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 
  #--on link scale
  prd2mb = prdMod(mdl2m,trms=c("(Intercept)","ti(z)"),lst=grids,type="link",keep="z");
  prd2m = prdMod(mdl2m,trms=NULL,lst=grids,type="link");
  ggplot(prd2m,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y)) + 
    geom_line() + geom_ribbon(colour=NA,alpha=0.3) + 
    geom_point(data=dfrESsp,mapping=aes(x=z,y=log(emp_sel),colour=y,size=n_BSFRF),inherit.aes=FALSE) +  
    scale_size_area() + 
    scale_y_continuous(limits=c(-5,5),oob=scales::oob_squish) +
    geom_line(data=prd2mb,mapping=aes(x=z,y=emp_sel),inherit.aes=FALSE) + 
    geom_ribbon(data=prd2mb,mapping=aes(x=z,y=emp_sel,ymin=lci,ymax=uci),inherit.aes=FALSE,alpha=0.3); 
  
  lstPrd2m = list();
  for (y_ in levels(dfrESs$y))
    lstPrd2m[[y_]] = prdMod(mdl2m,trms=paste0("ti(z):y",y_),lst=grids,type="link") |> 
                       dplyr::filter(y==as.numeric(y_));
  dfrPrd2m = dplyr::bind_rows(lstPrd2m);
  dfrPrd2m |> dplyr::group_by(z) |> dplyr::summarize(sum=sum(emp_sel)); #--check sum to zero
  dfrPrd2m |> dplyr::group_by(y) |> dplyr::summarize(sum=sum(emp_sel)); #--check sum to zero
  ggplot(dfrPrd2m,aes(x=z,y=emp_sel,colour=y,fill=y)) + geom_line()
  
  #--add in intercept and main smooth terms
  lstPrd2m[["I"]] = prdMod(mdl2m,trms=paste0("(Intercept)"),lst=grids,type="link");
  lstPrd2m[["ti(z)"]] = prdMod(mdl2m,trms=paste0("ti(z)"),lst=grids,type="link");
  ggplot(dfrPrd2m,aes(x=z,y=emp_sel,colour=terms,fill=terms)) + geom_line();
  
  dfrPrd2m = dplyr::bind_rows(lstPrd2m);
  dfrPrd2mp = dfrPrd2m |> dplyr::group_by(z,y) |> dplyr::summarize(link=sum(emp_sel)) |> dplyr::ungroup();
  dfrPrd2ma = prdMod(mdl2m,trms=NULL,lst=grids,type="link");
  ggplot(dfrPrd2ma,aes(x=z,y=emp_sel,colour=term,fill=term)) + geom_line()
   