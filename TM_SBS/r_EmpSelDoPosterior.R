require(dplyr)

getEmpSelFits<-function(m,
                        ds=NULL,
                        scale_="response",#--could be "response" or "link"
                        terms_="all",     #--smooth term names
                        ci_level_=0.90,
                        seed_=100,
                        n_samples_=50000,
                        post_method="mh",
                        n_cores_=4,
                        burnin_=1000,
                        thin_=1,
                        t_df_=40,
                        rw_scale_=0.25){
  if (is.null(ds)) 
    ds <- gratia::data_slice(m, 
                             z=gratia::evenly(z,lower=27.5,upper=182.5,by=5),
                             y=levels(y)) |>
            mutate(.row = row_number()); #--.row can be used as a join by variable
  
  if (any(terms_=="all")) terms_ = NULL;
  
  #--calculate fitted values
  fv <- gratia::fitted_values(m, 
                              terms=terms_,
                              data=ds,
                              scale=scale_,
                              ci_level=ci_level_,
                              unconditional=TRUE);
  
  #--draw samples from full posterior distribution
  ps <- gratia::posterior_samples(m, 
                                  terms=terms_,
                                  n=n_samples_, 
                                  data=ds, 
                                  seed=seed_,
                                  unconditional=TRUE,
                                  method=post_method,
                                  n_cores=n_cores_,
                                  burnin=burnin_,
                                  thin=thin_,
                                  t_df=t_df_,
                                  rw_scale=rw_scale_) |>
    left_join(ds, by = join_by(.row == .row));
  
  q_fun <- function(x,ci_level, ...) {
    lpi = (1-ci_level)/2;
    upi = 1-lpi;
    probs = c(.lower_pi=lpi,.med=0.5,.upper_pi=upi)
    tibble::tibble(
      .value = quantile(x, probs = probs, ...),
      .q = probs * 100,
      .qty = names(probs)
    )
  }

  #--calc quantiles on posterior samples
  qs = ps |>  dplyr::group_by(.row) |>
              dplyr::reframe(q_fun(.response,ci_level_)) |>
              tidyr::pivot_wider(id_cols = .row, 
                                 names_from = .qty, 
                                 values_from = .value);
  
  dfrMR = fv |> dplyr::left_join(qs,by=dplyr::join_by(.row==.row));
  return(list(dfrMR=dfrMR,dfrPS=ps));
}

plotP1<-function(resLst,n_curves,limits_=c(0,2)){
  n_samples_  = length(unique(resLst$dfrPS$`.draw`));
  p1 = ggplot(resLst$dfrMR,mapping=aes(x=z)) + 
        geom_line(data=resLst$dfrPS |> dplyr::filter(.draw %in% sample(1:n_samples_,n_curves)),
                  aes(y=.response,group=.draw),linewidth=0.1) + 
          # add the lower and upper prediction intervals
          geom_line(aes(y=.lower_pi), colour="#56B4E9",linewidth=1.5) +
          geom_line(aes(y=.upper_pi), colour="#56B4E9",linewidth=1.5) +
          geom_line(aes(y=.med),      colour="#56B4E9",linewidth=1.5) +
          # add the lower and upper credible intervals
          geom_line(aes(y=.lower_ci),colour="green",linewidth=0.5,linetype=2) +
          geom_line(aes(y=.upper_ci),colour="green",linewidth=0.5,linetype=2) +
          geom_line(aes(y=.fitted),colour="green",linewidth=0.5,linetype=1) +
        facet_wrap(~y) + 
        scale_y_continuous(limits=limits_,oob=scales::squish) +
        wtsPlots::getStdTheme();
  return(p1);
}
  
plotP2<-function(resLst,dfr,limits_=c(0,2)){
  p2 = ggplot(resLst$dfrMR,aes(x=z,y=.fitted)) +
      # summarise the posterior samples
      geom_bin_2d(data=resLst$dfrPS, aes(x=z, y=.response, fill=after_stat(ncount), group=y),
                  binwidth=c(5,0.05), alpha=0.7) +
      # add the lower and upper prediction intervals
      geom_line(aes(y=.lower_pi), colour="#56B4E9",linewidth=1.5) +
      geom_line(aes(y=.upper_pi), colour="#56B4E9",linewidth=1.5) +
      geom_line(aes(y=.med),      colour="#56B4E9",linewidth=1.5) +
      # add the lower and upper credible intervals
      geom_line(aes(y=.lower_ci), linewidth=0.5,linetype=2) +
      geom_line(aes(y=.upper_ci), linewidth=0.5,linetype=2) +
      # add the fitted model
      geom_line() +
      # add the observed data
      geom_point(data=dfr, aes(x=z, y=emp_sel, size=n_BSFRF/mean(n_BSFRF)), 
                 colour="green",shape=21,stroke=1) +
      scale_size_area(name="weight") + 
      scale_fill_viridis_c(name="posterior\ndensity",option="plasma",limits=c(0,0.25),oob=scales::squish) +
      scale_y_continuous(limits=limits_,oob=scales::squish) + 
      facet_wrap(~y) + 
      wtsPlots::getStdTheme() + 
      #theme(legend.position="none") +
      labs(x="size (mm CW)", y="Response",colour="year");
  return(p2);
}

#--get "raw" empirical selectivity functions
dfrESs = wtsUtilities::getObj(file.path(child_path$peek(),"rda_EmpSelFcnsRaw.RData"));

#--fit male data
dfrESsp = dfrESs |> dplyr::filter(x=="male",!is.nan(emp_sel));
mdl1m = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10),
                   family=tw,data=dfrESsp,method="REML",
                   weights=n_BSFRF/mean(n_BSFRF));
mdl2m = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,by=y,bs="sz"),
                   family=tw,data=dfrESsp,method="REML",select=TRUE,
                   weights=n_BSFRF/mean(n_BSFRF));
mdl3m = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,y,bs="fs"),
                   family=tw,data=dfrESsp,method="REML",select=TRUE,
                   weights=n_BSFRF/mean(n_BSFRF));

ds = tidyr::expand_grid(z=seq(from=27.5,to=min(182.5,max(dfrESsp$z)+5),by=5.0),
                        y=factor(levels(dfrESsp$y))) |> 
       dplyr::mutate(.row=dplyr::row_number());
resLstM1 = getEmpSelFits(mdl1m,ds); plotP1(resLstM1,100,c(0,2)); plotP2(resLstM1,dfrESsp,c(0,2));
resLstM2 = getEmpSelFits(mdl2m,ds); plotP1(resLstM2,100,c(0,2)); plotP2(resLstM2,dfrESsp,c(0,2));
resLstM2a = getEmpSelFits(mdl2m,ds,terms_=c("(Intercept)","ti(z)")); 
plotP1(resLstM2a,100,c(0,2)); plotP2(resLstM2a,dfrESsp,c(0,2));
resLstM3 = getEmpSelFits(mdl3m,ds); plotP1(resLstM3,100,c(0,2)); plotP2(resLstM3,dfrESsp,c(0,2));
resLstM3a = getEmpSelFits(mdl3m,ds,terms_=c("(Intercept)","ti(z)")); 
plotP1(resLstM3a,100,c(0,2)); plotP2(resLstM3a,dfrESsp,c(0,2));

#--fit female data
dfrESsp = dfrESs |> dplyr::filter(x=="female",!is.nan(emp_sel));
mdl1f = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10),
                   family=tw,data=dfrESsp,method="REML",
                   weights=n_BSFRF/mean(n_BSFRF));
mdl2f = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,by=y,bs="sz"),
                   family=tw,data=dfrESsp,method="REML",select=TRUE,
                   weights=n_BSFRF/mean(n_BSFRF));
mdl3f = mgcv::gam(emp_sel ~ ti(z,bs="tp",k=10) + ti(z,y,bs="fs"),
                   family=tw,data=dfrESsp,method="REML",select=TRUE,
                   weights=n_BSFRF/mean(n_BSFRF));

ds = tidyr::expand_grid(z=seq(from=27.5,to=min(182.5,max(dfrESsp$z)+5),by=5.0),
                        y=factor(levels(dfrESsp$y))) |> 
       dplyr::mutate(.row=dplyr::row_number());
resLstF1 = getEmpSelFits(mdl1f,ds); plotP1(resLstF1,100,c(0,1.0)); plotP2(resLstF1,dfrESsp,c(0,1.0));
resLstF2 = getEmpSelFits(mdl2f,ds); plotP1(resLstF2,100,c(0,1.0)); plotP2(resLstF2,dfrESsp,c(0,1.0));
resLstF2a = getEmpSelFits(mdl2f,ds,terms_=c("(Intercept)","ti(z)")); 
plotP1(resLstF2a,100,c(0,1.0)); plotP2(resLstF2a,dfrESsp,c(0,1.0));
resLstF3 = getEmpSelFits(mdl3f,ds); plotP1(resLstF3,100,c(0,1.0)); plotP2(resLstF3,dfrESsp,c(0,1.0));
resLstF3a = getEmpSelFits(mdl3f,ds,terms_=c("(Intercept)","ti(z)")); 
plotP1(resLstF3a,100,c(0,1.0)); plotP2(resLstF3a,dfrESsp,c(0,1.0));

wtsMGCV::gam.check.plots(mdl1f);
