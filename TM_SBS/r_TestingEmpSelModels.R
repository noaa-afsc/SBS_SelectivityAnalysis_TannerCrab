fitMod<-function(dfr,dfrp,frmla,fam,p=0.05){
  mdl = mgcv::gam(frmla,family=fam,data=dfr,weights=n/mean(n));
  prd = dplyr::bind_cols(
            dfrp,
            tibble::as_tibble(
              mgcv::predict.gam(mdl,dfrp,type="response",se.fit=TRUE),
            ) |> 
            dplyr::mutate(type="fit",
                          lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                          uci=qnorm(p,fit,se.fit,lower.tail=FALSE))
        );
  prd = dplyr::bind_rows(dfr,prd |> dplyr::rename(emp_sel=fit));
  return(list(mdl=mdl,prd=prd));
}
prdMod<-function(mdl,trms,lst,type="response",keep=NULL,p=0.05){
  dfr = wtsMGCV::createGridTbl(lst);
  prd = dplyr::bind_cols(
            dfr,
            tibble::as_tibble(
              mgcv::predict.gam(mdl,dfr,type=type,terms=trms,se.fit=TRUE),
            ) |> 
            dplyr::mutate(type="fit",
                          lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                          uci=qnorm(p,fit,se.fit,lower.tail=FALSE))
        ) |> 
          dplyr::rename(emp_sel=fit);
  if (!is.null(keep)){
    prd = prd |> 
            dplyr::distinct(pick(tidyselect::any_of(keep),emp_sel,se.fit,lci,uci));
    drop = names(lst)[!(names(lst) %in% keep)];
    for (drp in drop) prd[[drp]] = NA;
  }
  return(prd);
}
mdl = mdlESMsp1;
prd0 = prdMod(mdl,trms=c("(Intercept)","ti(z)"),lst=grids,type="response",keep="z"); plotMod(prd0,c(0,1.5))
prd1 = prdMod(mdl,trms=NULL,lst=grids,keep=NULL); plotMod(prd1)
prd2 = prdMod(mdl,trms=c("ti(z):y2015"),lst=grids,type="iterms",keep=NULL); plotMod(prd2)
plotMod(dplyr::bind_rows(prd0,prd1))
plotMod<-function(tmp,ylims=c(0,1.5)){
  if (all(is.na(tmp$y))) tmp$y = "all";
  p = ggplot(tmp,aes(x=z,y=emp_sel,ymin=lci,ymax=uci,colour=y,fill=y));
  if ("n" %in% names(tmp)){
    p = p + geom_point(aes(size=n)) + scale_size_area() + 
            geom_line();
  }
  p = p + 
#         geom_line(data=tmp |> dplyr::filter(type=="fit")) + 
         geom_ribbon(alpha=0.3) + 
#         geom_ribbon(data=tmp |> dplyr::filter(type=="fit"),alpha=0.3) + 
         geom_hline(yintercept=0.5,linetype=3) + 
         scale_y_continuous(limits=ylims,oob=scales::squish) + 
         labs(x="size (mm CW)",y="empirical\nselectivity",
              colour="study\nyear",fill="study\nyear",size="crab\nsampled") + 
         # wtsPlots::getStdTheme() + 
         theme(legend.position=c(0.99,0.99),
               legend.justification=c(1,1),
               legend.byrow=TRUE,
               legend.box="horizontal");
  return(p);
}
p = 0.05;#--tail probability for qnorm: confidence interval = 1-2*p
dfrESsp = dfrESs |> dplyr::mutate(y=factor(y));
grids = list(z=(dfrESsp |> dplyr::distinct(z))[[1]],
             y=(dfrESsp |> dplyr::distinct(y))[[1]]);
grdTbl = wtsMGCV::createGridTbl(grids)

lstESM0g = fitMod(dfr=dfrESsp |> dplyr::filter(x=="male"),
                  dfrp=dfrESsp |> dplyr::distinct(z) |> dplyr::mutate(y=NA,x="male"),
                  frmla = emp_sel ~ te(z,bs="tp",k=10),
                  fam=gaussian);
plotMod(lstESM0g$prd);
prdSMs = wtsMGCV::predSmoothTerms(lstESM0g$mdl,grids);
wtsMGCV::plotSmoothTerms(prdSMs,labs=list(z="z"),ci=0.90);
lstESM0t = fitMod(dfr=dfrESsp |> dplyr::filter(x=="male"),
                  dfrp=dfrESsp |> dplyr::distinct(z) |> dplyr::mutate(y=NA,x="male"),
                  frmla = emp_sel ~ te(z,bs="tp",k=10),
                  fam=tw);
plotMod(lstESM0t$prd)
prdSMs = wtsMGCV::predSmoothTerms(lstESM0t$mdl,grids);
wtsMGCV::plotSmoothTerms(prdSMs,labs=list(z="z"),ci=0.90);
lstESM1g = fitMod(dfr=dfrESsp |> dplyr::filter(x=="male"),
                  dfrp=dfrESsp |> dplyr::distinct(z) |> dplyr::mutate(y=NA,x="male"),
                  frmla=emp_sel ~ s(z,bs="re",k=20),
                  fam=gaussian);
plotMod(lstESM1g$prd)
wtsMGCV::gam.check.plots(lstESM1g$mdl);
lstESM1t = fitMod(dfr=dfrESsp |> dplyr::filter(x=="male"),
                  dfrp=dfrESsp |> dplyr::distinct(z) |> dplyr::mutate(y=NA,x="male"),
                  frmla=emp_sel ~ s(z,bs="re",k=20),
                  fam=tw);
plotMod(lstESM1t$prd);
wtsMGCV::gam.check.plots(lstESM1t$mdl);
wtsMGCV::predSmoothTerms(z=dfrESsp |> distinct(z),y=dfrESsp |> dplyr::distinct(y))
gratia::appraise(lstESM1t$mdl);
  mdlESMsp1 = mgcv::gam(emp_sel ~ 1+ti(z,bs="tp",k=20) + ti(z,by=y,bs="tp",k=20),
                       family=tw,
                       data=dfrESs |> dplyr::filter(x=="male") |> dplyr::mutate(y=factor(y)),
                       weights=n/mean(n),method="ML");
  mdlESMsp1a = mgcv::gam(emp_sel ~ y + s(z,y,bs="fs",k=20),
                       family=tw,
                       data=dfrESs |> dplyr::filter(x=="male") |> dplyr::mutate(y=factor(y)),
                       weights=n/mean(n),method="REML");
prd0 = prdMod(mdlESMsp1a,trms=c("s(z,y)1"),lst=grids,type="link"); plotMod(prd0,c(-1.5,1.5))


  mdlESMsp1b = mgcv::gam(emp_sel ~ s(z,y,bs="fs",k=20),
                       family=tw,
                       data=dfrESs |> dplyr::filter(x=="male") |> dplyr::mutate(y=factor(y)),
                       weights=n/mean(n),method="ML");
  mdlESMsp1c = mgcv::gam(emp_sel ~ s(z,k=20) + ti(z,y,bs="fs",k=20),
                       family=tw,
                       data=dfrESs |> dplyr::filter(x=="male") |> dplyr::mutate(y=factor(y)),
                       weights=n/mean(n),method="ML");
  prd = mgcv::predict.gam(mdlESMsp1c,newdata=grdTbl,type="response",se.fit=TRUE,
                          terms="s(z)",unconditional=TRUE) |> 
        tibble::tibble(f)

  mdlESMsp2 = mgcv::gam(emp_sel ~ s(z,y,bs="fs",k=10),
                       family=gaussian,
                       data=dfrESs |> dplyr::filter(x=="male") |> dplyr::mutate(y=factor(y)),
                       weights=n/mean(n),method="ML");
AIC(mdlESMsp0,mdlESMsp0a,mdlESMsp1,mdlESMsp1a,mdlESMsp2)
wtsMGCV::gam.check.plots(mdlESMsp0at)
wtsMGCV::
plot(mdlESMsp0a,se=TRUE,seWithMean=TRUE,residuals=TRUE,unconditional=TRUE)  
gam.check(mdlESMsp)
k.check(mdlESMsp)
anova(mdlESMsp)
summary(mdlESMsp)
pred<-function(mdl){
  dfrPrd = dplyr::bind_cols(
              tibble::as_tibble(
                mgcv::predict.gam(mdl,dfrESs |> dplyr::distinct(z,y) |> dplyr::mutate(y=factor(y)),
                                  type="response",se.fit=TRUE),
              ), 
              dfrESs |> dplyr::distinct(z,y) |> dplyr::mutate(y=factor(y))) |> 
              dplyr::mutate(type="fit",
                            lci=qnorm(p,fit,se.fit,lower.tail=TRUE),
                            uci=qnorm(p,fit,se.fit,lower.tail=FALSE));
  return(dfrPrd)
}
ggplot(pred(mdlESMsp0b),aes(z,fit,ymin=lci,ymax=uci,colour=y,fill=y)) + 
    geom_ribbon(alpha=0.3) + geom_line()

