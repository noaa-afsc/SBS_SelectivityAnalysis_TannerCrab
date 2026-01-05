#--plot model results------------------------------------------------------
#' 
#' @title Plot 1D model smooth
#' @description Function to plot 1D model smooth.
#' @param sme - smooth estimates from [gratia::smooth_estimates()], with confidence intervals
#' @param prs - dataframe with partial residuals from [gratia::add_partial_residuals()]
#' @param x_ - name of x-axis column in `sme`
#' @param y_ - name of value (y-axis) column in `sme`
#' @return a ggplot2 plot object 
#' @import ggplot2
#' @import gratia 
#' @import rlang
#' @export
#'  
plt1D<-function(sme,prs,x_,y_){
  require(rlang);
  vx=sym(x_); vy=sym(y_);
  gratia::draw(sme) + 
        geom_point(aes(x=!!vx,y=!!vy),prs,alpha=0.2);
}
#' 
#' @title Plot 2D model smooth
#' @description Function to plot 2D model smooth.
#' @param sme - smooth estimates from [gratia::smooth_estimates()], with confidence intervals
#' @param prs - dataframe with partial residuals from [gratia::add_partial_residuals()]
#' @param x_ - name of x-axis column in `sme`
#' @param y_ - name of y-axis column in `sme`
#' @param z_ - name of value (z-axis) column in `sme`
#' @return a ggplot2 plot object 
#' @import ggplot2
#' @import gratia 
#' @import rlang
#' @export
#'  
plt2D<-function(sme,prs,x_,y_,z_){
  require(rlang);
  vx=sym(x_); vy=sym(y_); vz=sym(z_);
  gratia::draw(sme) + 
        geom_point(aes(x=!!vx,y=!!vy,size=!!vz),prs,alpha=0.2) + 
        scale_size_area()
}
#' 
#' @title Plot model smooths
#' @description Function to plot model smooths.
#' @param mdl - model to plot 
#' @param subs - named character vector mapping covariates (names) to labels (values) 
#' @return nested list of a ggplot object and caption for each smooth 
#' @import gratia
#' @import mgcv
#' @export
#' 
plotModelSmooths<-function(mdl,
                           subs = c("z"="size (mm CW)",
                                    "d"="depth (m)",
                                    "t"="temperature (deg C)",
                                    "f"="phi",
                                    "s"="sorting coefficient")){
  resp = as.character(formula(mdl))[2];     #--model response
  fmly = family(mdl);
  nsms  = n_smooths(mdl);
  smths = smooths(mdl);
  smtsl = smooth_terms(mdl);
  dfrDatp = mdl$model;
  allPlts = list();
  for (ism in 1:nsms){
    #--testing: ism = 1;
    smth = smths[ism];
    smts = smtsl[[ism]];
    sme = gratia::smooth_estimates(mdl,select=smth,
                                   unconditional=TRUE,
                                   overall_uncertainty=TRUE);
    sme= add_confint(sme,0.80);
    prs = add_partial_residuals(dfrDatp,mdl,select=smth);
    plts = list();
    if (length(smts)==1){
      p = plt1D(sme,prs,smts[1],smth) + 
                    labs(x=subs[smts[1]],y="Partial effect on ln-scale catchability");
      cap = paste0("Partial effect of ",subs[smts[1]]," on ln-scale catachability. ",
                   "Partial residuals for the response '",
                   resp,"' are shown on the link scale of the model, based on its ",
                   fmly$family," error distribution and ",fmly$link," link.")
      plts[[1]] = list(p=p,cap=cap);
      if (smth %in% c("ti(z)","s(z)")) {
        p = gratia::draw(transform_fun(sme,fun=exp,constant=model_constant(mdl))) + 
                      labs(x=subs[smts[1]],y="Base Selectivity");
        cap = paste0("Estimated base selectivity (excluding covariate effects).");
        plts[[2]] = list(p=p,cap=cap);
      } else {
        # plts[[2]] = gratia::draw(transform_fun(sme,fun=inv_link(mdl))) + 
        #               labs(x=subs[smts[1]],y="Response Scale Effect");
        # cap = paste0("Estimated selectivity (independent of covariate effects.");
        # plts[[2]] = list(p=p,cap=cap);
      }
    } else if (length(smts==2)){
      p = plt2D(sme,prs,smts[1],smts[2],smth) + 
                    labs(x=subs[smts[1]],y=subs[smts[2]],size="Partial\nresiduals");
      cap = paste0("Partial effects of ",subs[smts[1]]," and ",subs[smths[2]]," on ln-scale catchbility. ",
                   "Partial residuals for the response '",
                   resp,"' are shown on the link scale of the model, based on its ",
                   fmly$family," error distribution and ",fmly$link," link.");
      plts[[1]] = list(p=p,cap=cap);
    }
    allPlts[[smth]] = list(plts);
    rm(plts);
  }
  return(allPlts)
}
