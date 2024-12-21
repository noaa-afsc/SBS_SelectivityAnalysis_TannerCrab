#--create list with censored data and prediction grids

dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

#--get trimmed data----
lst = wtsUtilities::getObj(file.path(dirThs,"../rda_Step2_TrimmedDataList.RData"));
dfrDat =  lst$lstTrimmedFinal$dfrDat |> 
            dplyr::mutate(obsR=exp(lnR),
                          nN=round(p*n), #--number caught in NMFS gear
                          nB=n-nN);      #--number caught in BSFRF gear

#--extract characteristics from size and environmental data----
##--size----
grd_z = seq(5,180,5)+2.5; med_z = 100.0; #--for males
##--depth----
med_d = median(dfrDat$d,na.rm=TRUE); rng_d = range(dfrDat$d,na.rm=TRUE);
grd_d = seq(from=rng_d[1],rng_d[2],length.out=50);
##--temperature----
med_t = median(dfrDat$t,na.rm=TRUE); rng_t = range(dfrDat$t,na.rm=TRUE);
grd_t = seq(from=rng_t[1],rng_t[2],length.out=50);
##--phi----
med_f = median(dfrDat$f,na.rm=TRUE); rng_f = range(dfrDat$f,na.rm=TRUE);
grd_f = seq(from=rng_f[1],rng_f[2],length.out=50);
##--sorting----
med_s = median(dfrDat$s,na.rm=TRUE); rng_s = range(dfrDat$s,na.rm=TRUE);
grd_s = seq(from=rng_s[1],rng_s[2],length.out=50);

meds  = list(z=med_z,d=med_d,t=med_t,f=med_f,s=med_s);
grids = list(z=grd_z,d=grd_d,t=grd_t,f=grd_f,s=grd_s);
lbls  = list(z="size (mm CW)",d="depth (m)",t="temperature (deg C)",f="-log2(phi)",s="sorting coeff.")

#--the logit-scale observed proportions "o=nN/(nN+nB)" are related to the 
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q) 
# where q = expF_b/expF_n is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf), 
# and As = area swept and Sf is the sampling fraction (note that 
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).

#--set sex
x_ = "MALE";
dfrDat   = dfrDat |> dplyr::filter(x==x_) |> dplyr::mutate(obsR=exp(lnR));

lst = list(x=x_,dfrDat=dfrDat,meds=meds,grids=grids,lbls=lbls);
wtsUtilities::saveObj(lst,file.path(dirThs,"rda_Step3a.CensoredDataAndGridsList.Males.RData"));
