#--fit various models for ln(r) using mgcv to fit GAMs for a specific sex

#--load required packages
require(magrittr);
require(mgcv);
#--other required packages
#----dplyr
#----stats
#----stringr
#----tibble
#----wtsUtilities


#--the logit-scale observed proportions "o" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_a is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

#--extract data for a specific sex
x_ = "MALE";
dfrDatRE = wtsUtilities::getObj("../dfrTrimmedData.RData") %>% subset(x==x_);
dfrDatRE$h = as.factor(dfrDatRE$h);
wtsUtilities::saveObj(dfrDatRE,"dfrDatRE.RData");

#set up BINOMIAL model
fam = stats::binomial(link="logit");
mdlsB = list();

#----SIMPLE RE(h)----------------------------
#--ln(r) = ti(z) + RE(h)
frmla = p~ti(z,bs="ts") + s(h,bs="re");
mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=FALSE,scale=0,offset=lnq,method="REML");
summary(mdl);
#plot(mdl);
mdlsB[["z + RE(h)"]] = mdl;
rm(mdl);

#---------ALL 2-WAY INTERACTIONS + RE(h) ----------------------------------
#--ln(r) = ti(z) + ti(d) + ti(t) + ti(f)) + ti(s) + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        ti(d,t) + ti(d,f) + ti(d,s) + 
#--        ti(t,f) + ti(t,s) + 
#--        ti(f,s) +
#--        RE(h)
frmla  = p~ti(z,bs="ts")  +ti(d,bs="ts")  +ti(t,bs="ts")  +ti(f,bs="ts")  +ti(s,bs="ts")  +
           ti(z,d,bs="ts")+ti(z,t,bs="ts")+ti(z,f,bs="ts")+ti(z,s,bs="ts")+
           ti(d,t,bs="ts")+ti(d,f,bs="ts")+ti(d,s,bs="ts")+
           ti(t,f,bs="ts")+ti(t,s,bs="ts")+
           ti(f,s,bs="ts") +
           s(h,bs="re");
mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="REML");
summary(mdl);
#plot(mdl);
mdlsB[["all 2-way interactions + RE"]] = mdl;
rm(mdl);

#---------1st PASS SIGNIFICANT 2-WAY INTERACTIONS + RE(h) ----------------------------------
#--ln(r) = ti(z) + ti(d) + MISSING[ti(t)] + ti(f)) + MISSING[ti(s)] + 
#--        ti(z,d) + ti(z,t) + ti(z,f) + ti(z,s) + 
#--        MISSING[ti(d,t)] + MISSING[ti(d,f)] + ti(d,s) + 
#--        MISSING[ti(t,f)] + ti(t,s) + 
#--        ti(f,s) +
#--        RE(h)
frmla  = p~ti(z,bs="ts") + ti(d,bs="ts") + ti(f,bs="ts") +
           ti(z,d,bs="ts") + ti(z,t,bs="ts") + ti(z,f,bs="ts") + ti(z,s,bs="ts")+
           ti(d,s,bs="ts")+
           ti(t,s,bs="ts")+
           ti(f,s,bs="ts") +
           s(h,bs="re");
mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="REML");
summary(mdl);
#plot(mdl);
mdlsB[["1st pass signif 2-way interactions + RE"]] = mdl;
rm(mdl);

# #---------2nd PASS SIGNIFICANT 2-WAY INTERACTIONS + RE(h) ----------------------------------
# #--ln(r) = ti(z) + ti(d) + ti(t) + MISSING[ti(f))] + MISSING[ti(s)] + 
# #--        ti(z,d) + MISSING[ti(z,t)] + MISSING[ti(z,f)] + ti(z,s) + 
# #--        MISSING[ti(d,t)] + MISSING[ti(d,f)] + ti(d,s) + 
# ##--       MISSING[ti(t,f)] + MISSING[ti(t,s)] + MISSING[ti(f,s)] +
# #--        RE(h)
# frmla  = p~ti(z,bs="ts")+ti(d,bs="ts")+#ti(f,bs="ts")+ti(s,bs="ts")+
#            ti(z,d,bs="ts")+#ti(z,t,bs="ts")+ti(z,f,bs="ts")+
#            ti(z,s,bs="ts")+
#         #   ti(d,t,bs="ts")+ti(d,f,bs="ts")+
#            ti(d,s,bs="ts")+
#         #   ti(t,f,bs="ts")+ti(t,s,bs="ts")+
#         #   ti(f,s,bs="ts") +
#            s(h,bs="re");
# mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=TRUE,scale=0,offset=lnq,method="REML");
# summary(mdl);
# plot(mdl);
# mdlsB[["2nd pass signif 2-way interactions + RE"]] = mdl;
# rm(mdl);

#--compare models using information criteria
str_mdls = paste0("mdlsB[['",names(mdlsB),"']]",collapse=",");
eval(parse(text=paste0("aic = stats::AIC(",str_mdls,")")))
eval(parse(text=paste0("bic = stats::BIC(",str_mdls,")")))
dfrICs = rbind(aic %>% tibble::rownames_to_column(var="model") %>% dplyr::transmute(type="AIC",model=model,df=df,value=AIC),
               bic %>% tibble::rownames_to_column(var="model") %>% dplyr::transmute(type="BIC",model=model,df=df,value=BIC)) %>%
         tidyr::pivot_wider(names_from=type,values_from=value);
dfrICs$model %<>% stringr::str_sub(9,-4);
wtsUtilities::saveObj(dfrICs,"dfrICs.RData");
write.csv(dfrICs,file="dfrICs.csv");

wtsUtilities::saveObj(mdlsB,     paste0("selModels.BinomRE.All.", tolower(x),"s.RData"));
wtsUtilities::saveObj(mdlsB[[3]],paste0("selModels.BinomRE.Best.",tolower(x),"s.RData"));


