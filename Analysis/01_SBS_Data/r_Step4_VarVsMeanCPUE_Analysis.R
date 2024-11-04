#--analyze variance to mean relationships for CPUE

calcStep4<-function(){
  source("r_Functions-Figures-Abundance.R");
  
  #--output list
  out = list();
  
  #--load crab abundance/CPUE results
  lst = wtsUtilities::getObj("rda_Step3_SBS_CrabAbundance.RData");
  dfrStatsCPUE = lst$dfrStatsCPUE;
  
#| label: fig-VMRsCPUEMales
  cap = "VMRs for male CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment. The horizontal line indicates a VMR of 1."
  p = plotVMRsCPUE(dplyr::filter(dfrStatsCPUE,x=='male'),mn,vr);
  out = c(out,list(`fig-VMRsCPUEMales`=list(p=p,cap=cap)));
  
#| label: fig-VMRsCPUEFemales
  cap = "VMRs for female CPUE in the SBS studies, by 5-mm size bin and gear type. The dotted vertical line marks the lower limit of the size bins used in the assessment. The horizontal line indicates a VMR of 1."
  p = plotVMRsCPUE(dplyr::filter(dfrStatsCPUE,x=='female'),mn,vr);
  out = c(out,list(`fig-VMRsCPUEFemales`=list(p=p,cap=cap)));

#| label: fig-VarVsMeanCPUEBySex
  cap = "Variance in CPUE by 5 mm CW size bin plotted against mean CPUE, on the natural logarithm scale. Symbols: observed values (circles: females; triangles: males); colored shading/line: GAM smooth fits to observed values; dotted lines: linear regressions; solid black line: 1:1 line."
  p = plotVarVsMean(dplyr::filter(dfrStatsCPUE,z>22.5),label="");
  out = c(out,list(`fig-VarVsMeanCPUEBySex`=list(p=p,cap=cap)));

#| label: fig-VarVsMeanCPUEByFleet
  cap = "Variance in CPUE by 5 mm CW size bin plotted against mean CPUE, on the natural logarithm scale. Symbols: observed values (circles: females; triangles: males); colored shading/line: GAM smooth fits to observed values; dotted lines: linear regressions; solid black line: 1:1 line."
  p = plotVarVsMean(dplyr::filter(dfrStatsCPUE,z>22.5),facet=x,factor=fleet,label="");
  out = c(out,list(`fig-VarVsMeanCPUEByFleet`=list(p=p,cap=cap)));

  #--all
  mdlTWPf = mgcv::gam(log(vr)~fleet:x + lnmn + fleet:x:lnmn,
                      data=dplyr::filter(dfrStatsCPUE |> dplyr::mutate(lnmn=log(mn)),
                                         z>=22.5),
                      family=gaussian(),method="REML");
  mdlTWPr = mgcv::gam(log(vr)~fleet:x + lnmn,
                      data=dplyr::filter(dfrStatsCPUE |> dplyr::mutate(lnmn=log(mn)),
                                         z>=22.5),
                      family=gaussian(),method="REML");
  bicTWP = BIC(mdlTWPf,mdlTWPr) |> 
             dplyr::mutate(model=c("full","reduced"),
                           `delta(BIC)`=BIC-min(BIC,na.rm=TRUE)) |>
             dplyr::arrange(`delta(BIC)`) |> 
             dplyr::select(model,BIC,`delta(BIC)`)
  out = c(out,list(mdlsTWP=list(full=mdlTWPf,reduced=mdlTWPr,dfrBIC=bicTWP)));

  mdlSmryTWP = modelsummary::modelsummary(list(reduced=mdlTWPr,full=mdlTWPf),
                                          estimate="{estimate}",
                                          statistic=c("CI"="[{conf.low}, {conf.high}]"),
                                          shape=term~statistic,
                                          coef_map=c("lnmn"="TW exponent"),
                                          output="modelsummary_list")
  out = c(out,list(mdlSmryTWP=mdlSmryTWP));

  #| label: fig-QQ-mdlTWP
  cap = "Quantile-quantile plot of scaled residuals [@DHARMa] from the reduced model used to estimate the Tweedie power coefficient for the observed CPUE data."
  p = ggplotify::as.ggplot(~DHARMa::plotQQunif(mdlTWPr));
  out = c(out,list(`fig-QQ-mdlTWP`=list(p=p,cap=cap)));
  
  #| label: fig-DRs-mdlTWP
  cap = "Scaled residuals [@DHARMa] from the reduced model used to estimate the Tweedie power coefficient for the observed CPUE data."
  srs = DHARMa::simulateResiduals(mdlTWPr,plot=FALSE);
  p = ggplotify::as.ggplot(~DHARMa::plotResiduals(srs));
  out = c(out,list(`fig-DRs-mdlTWP`=list(p=p,cap=cap)));
  
  #--save objects
  wtsUtilities::saveObj(out,"rda_Step4_SBS_VarVsMeanCPUE_Analysis.RData");
  return(out);
}
out = calcStep4();
#rm(out);
