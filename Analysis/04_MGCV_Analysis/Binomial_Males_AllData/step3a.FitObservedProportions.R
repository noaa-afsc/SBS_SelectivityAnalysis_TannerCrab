#--fit various models for ln(r) using mgcv to fit GAMs for a specific sex
require(DHARMa);
require(mgcv);

dfrProps_AYs =  wtsUtilities::getObj("../rda_dfrPropsAllSBSYears.RData");

#--the logit-scale observed proportions "o" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_a is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

#--extract data for a specific sex
x = "MALE";
dfrDat = dfrProps_AYs |> 
         subset(SEX==x) |> 
         dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting,z=SIZE,
                          p=propNMFS,n=numTot,lnq=log(q),lgtp=log(p/(1-p))-lnq);
grd_z = seq(5,180,5)+2.5; med_z = 77.5;
med_d = median(dfrDat$d,na.rm=TRUE); rng_d = range(dfrDat$d,na.rm=TRUE); grd_d = seq(from=rng_d[1],rng_d[2],length.out=50);
med_t = median(dfrDat$t,na.rm=TRUE); rng_t = range(dfrDat$t,na.rm=TRUE); grd_t = seq(from=rng_t[1],rng_t[2],length.out=50);
med_f = median(dfrDat$f,na.rm=TRUE); rng_f = range(dfrDat$f,na.rm=TRUE); grd_f = seq(from=rng_f[1],rng_f[2],length.out=50);
med_s = median(dfrDat$s,na.rm=TRUE); rng_s = range(dfrDat$s,na.rm=TRUE); grd_s = seq(from=rng_s[1],rng_s[2],length.out=50);

#--check for outliers
dfrQs = tibble::tibble(q  =          c(  3,     6,      12),
                       y  =   1000 + c(  0,     0,       0),
                       lbl=as.factor(c("min","nominal","max")));
p = ggplot2::ggplot()+
      ggplot2::geom_histogram(data=dfrDat,mapping=ggplot2::aes(x=exp(lnq),fill=y)) +
      ggplot2::geom_vline(data=dfrQs,mapping=ggplot2::aes(xintercept=q),linetype=2) +
      ggplot2::geom_text(data=dfrQs,mapping=ggplot2::aes(x=q,y=y,label=lbl)) +
      labs(x="q = (A_NMSF/A_BSFRF)*(S_NMFS/S_BSFRF)",y="count",fill="year");
print(p);
nrw_orig = nrow(dfrDat);
p = ggplot2::ggplot()+
      ggplot2::geom_point(data=dfrDat,
                          mapping=ggplot2::aes(x=z,y=lgtp,colour=y)) +
      labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
print(p);
#--edit data
dfrDatp  = dfrDat |> subset(dplyr::between(exp(lnq),dfrQs$q[1],dfrQs$q[3]));
dfrDrop  = dfrDat |> anti_join(dfrDatp,by=NULL) |> dplyr::select(h) |> dplyr::distinct();
nrw_edit = nrow(dfrDat);
p = ggplot2::ggplot()+
      ggplot2::geom_point(data=dfrDatp,mapping=ggplot2::aes(x=z,y=lgtp,colour=y)) +
      ggplot2::geom_point(data=dfrDat |> anti_join(dfrDatp,by=NULL),mapping=ggplot2::aes(x=z,y=lgtp),colour="black") +
      labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
print(p);

#set up BINOMIAL model
fam = stats::binomial(link="logit");
mdlsB = list();

dfrNew = newDFR(zgrd,dmed,tmed,fmed,smed); 

#--ln(r) = ti(z)
frmla = p~ti(z,bs="ts");
mdl = estSelRatio(dfrDat,formula=frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["z"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+ti(f,bs="ts")+ti(s,bs="ts")+
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
           ti(d,t,bs="ts")+ti(d,f,bs="ts")+ti(d,s,bs="ts")+
           ti(t,f,bs="ts")+ti(t,s,bs="ts")+
           ti(f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["all 2-way interactions"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

#--ln(r) = ti(z) + ti(d) + ti(t) + MISSING[ti(f))] + ti(s) + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + ti(t,f) + ti(t,s) + ti(f,s)
frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+ti(s,bs="ts")+
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
           ti(d,t,bs="ts")+ti(d,f,bs="ts")+ti(d,s,bs="ts")+
           ti(t,f,bs="ts")+ti(t,s,bs="ts")+
           ti(f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=NULL);
mdlsB[["signif 2-way interactions"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) +ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s)
frmla  = p~ti(z,bs="ts")   + ti(d,bs="ts")   + ti(f,bs="ts")   + ti(s,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + ti(z,s,bs="ts") +
           ti(d,t,bs="ts") + ti(d,f,bs="ts") + ti(d,s,bs="ts") +
           ti(t,f,bs="ts") + ti(t,s,bs="ts")+
           ti(f,s,bs="ts") +
           ti(z,d,t,bs="ts") + ti(z,d,f,bs="ts") + ti(z,d,s,bs="ts") +
           ti(z,t,f,bs="ts") + ti(z,t,s,bs="ts") +
           ti(z,f,s,bs="ts") +
           ti(d,t,f,bs="ts") + ti(d,t,s,bs="ts") + ti(d,f,s,bs="ts") + ti(t,f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["all 3-way interactions"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + ti(s) +
#--        ti(z,d) + ti(z,t) + ti(z,f) + MISSING[ti(z,s)] + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + MISSING[ti(t,s)] + 
#--        MISSING[ti(f,s)]
frmla  = p~ti(z,bs="ts")   + ti(d,bs="ts")   + ti(f,bs="ts")   + ti(s,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + # ti(z,s,bs="ts") +
           ti(d,t,bs="ts") + ti(d,f,bs="ts") + ti(d,s,bs="ts") +
           ti(t,f,bs="ts") + # ti(t,s,bs="ts") +
           #ti(f,s,bs="ts") +
           ti(z,d,t,bs="ts") + ti(z,d,f,bs="ts") + ti(z,d,s,bs="ts") +
           ti(z,t,f,bs="ts") + ti(z,t,s,bs="ts") +
           ti(z,f,s,bs="ts") +
           ti(d,t,f,bs="ts") + ti(d,t,s,bs="ts") + ti(d,f,s,bs="ts") + ti(t,f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["signif 3-way interactions"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f) + MISSING[ti(s)] +
#--        ti(z,d) + ti(z,t) + ti(z,f) + MISSING[ti(z,s)] + 
#--        ti(d,t) + MISSING[ti(d,f)] + MISSING[ti(d,s)] + 
#--        ti(t,f) + MISSING[ti(t,s)] + 
#--        MISSING[ti(f,s)]
frmla  = p~ti(z,bs="ts")   + ti(d,bs="ts")   + ti(f,bs="ts")   + # ti(s,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + # ti(z,s,bs="ts") +
           ti(d,t,bs="ts") + # ti(d,f,bs="ts") + ti(d,s,bs="ts") +
           ti(t,f,bs="ts") + # ti(t,s,bs="ts") +
           #ti(f,s,bs="ts") +
           ti(z,d,t,bs="ts") + ti(z,d,f,bs="ts") + ti(z,d,s,bs="ts") +
           ti(z,t,f,bs="ts") + ti(z,t,s,bs="ts") +
           ti(z,f,s,bs="ts") +
           ti(d,t,f,bs="ts") + ti(d,t,s,bs="ts") + ti(d,f,s,bs="ts") + ti(t,f,s,bs="ts");
mdl  = estSelRatio(dfrDat,frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
mdlsB[["signif 3-way interactions: 2nd pass"]] = mdl;
summary(mdl$model);
plot(mdl$model);
rm(mdl);

wtsUtilities::saveObj(mdlsB,     paste0("selModels.Binom.All.",tolower(x),"s.RData"));
wtsUtilities::saveObj(mdlsB[[3]],paste0("selModels.Binom.Best.",tolower(x),"s.RData"));

stats::AIC(mdlsB[[1]]$model,mdlsB[[2]]$model,mdlsB[[3]]$model,mdlsB[[4]]$model,mdlsB[[5]]$model,mdlsB[[6]]$model);
stats::BIC(mdlsB[[1]]$model,mdlsB[[2]]$model,mdlsB[[3]]$model,mdlsB[[4]]$model,mdlsB[[5]]$model,mdlsB[[6]]$model);

# ################################################################################
# #--set up NEGATIVE BINOMIAL model
# fam = mgcv::nb(link="log");
# dfrDatNB = dfrDat |> dplyr::mutate(lgtP = log(p/(1-p)));
# mdlsNB = list();
# 
# #--ln(r) = ti(z)
# frmla = p~ti(z,bs="ts");
# mdl = estSelRatio(dfrDatNB,formula=frmla,family=fam,dfrNew=dfrNew,showPlots=TRUE);
# mdlsNB[["z"]] = mdl;
# summary(mdl$model);
# plot(mdl$model);
# rm(mdl);


