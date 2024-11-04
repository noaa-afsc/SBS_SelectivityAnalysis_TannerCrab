#--estimating Tweedie power coefficients for CPUE data


The variance in the "observed" CPUE was examined as a function of the mean by size bin for both gear types and sexes ([@fig-VarVsMeanCUE]) to estimate the corresponding variance-to-mean ratio (VMR, also known as the coefficient of dispersion) and determine whether it was a constant or it varied with the mean. Values for sizes less than 25 mm CW (the minimum size included in the stock assessment model) were excluded. Distributions consistent with the Tweedie exponential family (@Tweedie1984) exhibit variance-mean relationships of the form 

$$V =\sigma^2 \mu^p_{tw}$$ {#eq-TweedieVMRel}

where $V$ represents the variance of some random variable, $\mu$ represents its mean, $\sigma$ is a positive scaling coefficient, and $p_{tw}$ is the so-called Tweedie power parameter. Tweedie distributions are a subset of exponential dispersion models and include the normal ($p_{tw}=0$), Poisson  ($p_{tw}=1$), and gamma ($p_{tw}=2$) distributions, among others, as special cases (@ref). For $1 < p_{tw} < 2$, the distributions are compound Poisson-gamma distributions, with positive mass at 0 but otherwise positive and continuous. [@eq-TweedieVMRel] implies a linear relationship between $ln(V)$ and $ln(\mu)$:

$$ln(V) = \alpha + p_{tw} \cdot ln(\mu) $$ {#eq-TweedieLnVMRel}

[@fig-VarVsMeanCUE] suggests the relationship between the ln-scale variance and mean may be slightly convex when treating results across years and size bins as equivalent. We fit GAMs of the form $ln(V) \sim \alpha + \beta \cdot ln(\mu) + s(ln(\mu))$, where $\alpha$ and $\beta$ are fixed effects representing the ln-scale intercept and slope of the linear relationship between $ln(V)$ and the covariate $ln(\mu)$ expected for a Tweedie-distributed process and $s(ln(\mu))$ represents an additional smooth term, to the "observations" in [@fig-VarVsMeanCUE] by gear type and sex.

```{r}
#| label: fig-VarVsMeanCPUE
#| fig-cap: "Variance in CPUE by 5 mm CW size bin plotted against mean CPUE, on the natural logarithm scale. Symbols: observed values (circles: females; triangles: males); colored shading/line: GAM smooth fits to observed values; dotted lines: linear regressions; solid black line: 1:1 line."
  p = plotVarVsMean(dplyr::filter(dfrStatsCPUE,z>22.5));
  lstFigs = c(lstFigs,wtsQMD::printGGplot(p));
```

```{r SBSDataResults_AnalyzeCPUE}
fitModels<-function(full,dfr,keep="log(mn)"){
  fs = wtsUtilities::createFormulas(full,keep=keep);
  mdls = list();
  for (i in 1:length(fs)){
    mdl = lm(fs[[i]],data=dfr);
    if (mdl$df.residual>0) mdls[[as.character(fs[[i]][3])]] = mdl;
  }
  return(mdls)
}

#--formula to fit log-scale mean to log-scale variance for CPUE
full = log(vr)~y*z*log(mn);

#--males, BSFRF
mdlL = mgcv::gam(log(vr)~logmn,
                data=dplyr::filter(dfrStatsCPUE |> dplyr::mutate(logmn=log(mn)),
                                   fleet=="BSFRF",x=="male",z>=22.5),
                family=gaussian(),method="REML")
mdlS = mgcv::gam(log(vr)~logmn + ti(logmn,y,bs="fs"),
                data=dplyr::filter(dfrStatsCPUE |> dplyr::mutate(logmn=log(mn)),
                                   fleet=="BSFRF",x=="male",z>=22.5),
                family=gaussian(),method="REML")
AIC(mdlL,mdlS)
BIC(mdlL,mdlS)
DHARMa::plotQQunif(mdlS)
DHARMa::plotResiduals(mdlS)
summary(mdlS)
# dfr = dfrStatsCPUE |> dplyr::filter(fleet=="BSFRF",x=="male") |>
#                       dplyr::mutate(z=factor(z),y=factor(y));
# mdls = fitModels(full,dfr,keep=NULL)
# aicctbl = AICcmodavg::aictab(mdls)
# summary(mdls[[aicctbl$Modnames[1]]])
# bictbl = AICcmodavg::bictab(mdls)
# summary(mdls[[bictbl$Modnames[1]]])
# mdlstr = "y:log(mn)";
# DHARMa::plotQQunif(mdls[[mdlstr]])
# DHARMa::plotResiduals(mdls[[mdlstr]])

#--females, BSFRF
dfr = dfrStatsCPUE |> dplyr::filter(fleet=="BSFRF",x=="female") |>
                      dplyr::mutate(z=factor(z),y=factor(y));
mdls = fitModels(full,dfr,keep=NULL)
aicctbl = AICcmodavg::aictab(mdls)
summary(mdls[[aicctbl$Modnames[1]]])
bictbl = AICcmodavg::bictab(mdls)
summary(mdls[[bictbl$Modnames[1]]])
mdlstr = "log(mn)";
DHARMa::plotQQunif(mdls[[mdlstr]])
DHARMa::plotResiduals(mdls[[mdlstr]])

#--males, NMFS
dfr = dfrStatsCPUE |> dplyr::filter(fleet=="NMFS",x=="male") |>
                      dplyr::mutate(z=factor(z),y=factor(y));
mdls = fitModels(full,dfr,keep=NULL)
aicctbl = AICcmodavg::aictab(mdls);
summary(mdls[[aicctbl$Modnames[1]]])
bictbl = AICcmodavg::bictab(mdls);
summary(mdls[[bictbl$Modnames[1]]])
mdlstr = "log(mn)";
DHARMa::plotQQunif(mdls[[mdlstr]])
DHARMa::plotResiduals(mdls[[mdlstr]])

#--females, NMFS
dfr = dfrStatsCPUE |> dplyr::filter(fleet=="NMFS",x=="female") |>
                      dplyr::mutate(z=factor(z),y=factor(y));
mdls = fitModels(full,dfr,keep=NULL)
aicctbl = AICcmodavg::aictab(mdls);
summary(mdls[[aicctbl$Modnames[1]]]);
bictbl = AICcmodavg::bictab(mdls);
summary(mdls[[bictbl$Modnames[1]]]);
mdlstr = "log(mn)";
DHARMa::plotQQunif(mdls[[mdlstr]]);
DHARMa::plotResiduals(mdls[[mdlstr]]);

```
