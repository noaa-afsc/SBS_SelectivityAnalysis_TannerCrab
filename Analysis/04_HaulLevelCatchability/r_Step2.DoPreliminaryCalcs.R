#--subset haul-level data----

#--read project setup info----
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);
fn = file.path(dirPrj,"rda_ProjectSetup.RData");
s  = wtsUtilities::getObj(fn);

#--get proportions dataframe
dfrPropsAllYears   = wtsUtilities::getObj(file.path(dirThs,"rda_Step1b_dfrPropsAllSBSYears.RData"));

#--the logit-scale observed proportions "o" are related to the
#--NMFS-to-BSFRF selectivity ratio "r" (S_nmfs/S_bsfrf) by
#     logit(o) = ln(r) + ln(q)
# where q = expF_b/expF_a is the ratio of expansion factors used to convert
# numbers sampled to CPUE: i.e., CPUE = N_s * expF, where expF = 1/(As * Sf),
# and As = area swept and Sf is the sampling fraction (note that
# the SAMPLING_FACTOR in the haul data tables is 1/(sampling fraction),
# so expF = SAMPLING_FACTOR/AREA_SWEPT).
#--extract proportions data with environmental covariates
dfrDat = dfrPropsAllYears |>
         dplyr::transmute(y=YEAR,h=HAULJOIN,d=BOTTOM_DEPTH,t=GEAR_TEMPERATURE,f=phi,s=sorting,
                          z=SIZE,x=SEX,p=propNMFS,n=numTot,asr=ratioAS,q=q,lnq=log(q),lgtp=log(p/(1-p)),lnR=lgtp-lnq);

##--function to check for outliers based on area swept ratios----
#'
#' @title Edit data for outlier areas swept ratios (NMFS/BSFRF ratios)
#'
#' @description Function to edit data for outlier areas swept ratios.
#'
#' @param dfrDat - tibble with proprtion data to edit
#' @param qs - 3-element vector with min ASR, nominal ASR, max ASR
#'
#' @return list
#'
#' @export
#'
check_ASRs<-function(dfrDat,ASRs,maxASR=20){
  #--
  dfrDatp = dfrDat |> dplyr::distinct(y,h,asr) |>
              dplyr::mutate(asr=ifelse(asr>maxASR,maxASR,asr));
  #--create dataframe for q's defining range
  gghst = ggplot2::stat_bin(data=dfrDatp,mapping=ggplot2::aes(x=asr,fill=y),binwidth=1);
  dfrASRs = tibble::tibble(asr=c(ASRs[1]-0.5,ASRs[2],ASRs[3]+0.5),
                           y=1000+0*ASRs[1:3],
                           lbl=as.factor(c("min","nominal","max")));
  qts = quantile(dfrDatp$asr,probs=c(0.98),na.rm=TRUE);
  #--plot histogram of ASR's
  p1 = ggplot2::ggplot()+
        ggplot2::geom_histogram(data=dfrDatp,mapping=ggplot2::aes(x=asr,fill=y),binwidth=1) +
        ggplot2::geom_vline(data=dfrASRs,mapping=ggplot2::aes(xintercept=asr),linetype=2) +
        # ggplot2::geom_text(data=dfrQs,mapping=ggplot2::aes(x=q,y=y,label=lbl)) +
        ggplot2::scale_colour_viridis_d(option="plasma") +
        ggplot2::xlim(c(0,maxASR)) +
        ggplot2::labs(x="area swept ratio (A_NMSF/A_BSFRF)",y="count",fill="year");
  # print(p1);

  #--plot original data on logit scale
  nrw_orig = nrow(dfrDat);
  dfrDatp = dplyr::bind_rows(
                  dfrDat |> dplyr::filter(dplyr::between(asr,ASRs[1],ASRs[3])) |>
                            dplyr::mutate(cat=paste0("in [",ASRs[1],",",ASRs[3],"]")),
                  dfrDat |> dplyr::filter(asr<ASRs[1]) |>
                            dplyr::mutate(cat=paste0("<",ASRs[1])),
                  dfrDat |> dplyr::filter(asr>ASRs[3]) |>
                            dplyr::mutate(cat=paste0(">",ASRs[3])));
  dfrDatp$cat = factor(dfrDatp$cat,levels=c(paste0("<",ASRs[1]),
                                            paste0("in [",ASRs[1],",",ASRs[3],"]"),
                                            paste0(">",ASRs[3])));
  p2 = ggplot2::ggplot(data=dfrDatp,mapping=ggplot2::aes(x=z,y=lnR,fill=cat,size=n))+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=1),alpha=0.8,shape=21) +
        ggplot2::scale_colour_viridis_d(option="plasma") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        ggplot2::labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)",fill="q");
  # print(p2);

  #--edit data based on ASR limits
  dfrDatp  = dfrDat |> dplyr::filter(dplyr::between(asr,ASRs[1],ASRs[3]));
  dfrDrop  = dfrDat |> dplyr::anti_join(dfrDatp,by=NULL) |> dplyr::select(h) |> dplyr::distinct();
  nrw_edit = nrow(dfrDatp);
  nrw_drop = nrow(dfrDrop);
  p3 = ggplot2::ggplot()+
        ggplot2::geom_point(data=dfrDatp,
                            mapping=ggplot2::aes(x=z,y=lgtp,colour=y),
                            position=ggplot2::position_jitter(width=1)) +
        ggplot2::geom_point(data=dfrDat |> dplyr::anti_join(dfrDatp,by=NULL),
                            mapping=ggplot2::aes(x=z,y=lnR),colour="black") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  # print(p3);

  return(list(dfrDat=dfrDatp,dfrDrop=dfrDrop,ps=list(p1=p1,p2=p2,p3=p3),
              nrw_orig=nrw_orig,nrw_edit=nrw_edit,nrw_drop=nrw_drop));
}

##--function to select data based on number of individuals caught----
#'
#' @title Select data based on minimium number of individuals caught
#'
#' @description Function to select data based on minimium number of individuals caught.
#'
#' @param dfrDat - tibble with proportion data to edit
#' @param n_min - minimum number of individuals required to keep
#'
#' @return list
#'
#' @export
#'
selectSizeData<-function(dfrDat,n_min){
  #--subset the data based on minimum number of individuals caught
  dfrDatpp = dfrDat |> dplyr::filter(n>=n_min);

  #--plot the data that have been subsetted
  p1 = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=lnR,colour=y,size=n))+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=1)) +
        ggplot2::geom_point(data=dfrDat |> dplyr::anti_join(dfrDatpp,by=NULL),mapping=ggplot2::aes(x=z,y=lnR,size=n),
                            position=position_jitter(width=1),colour="black") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  # print(p1);

  #--plot data that will be fit
  p2 = ggplot2::ggplot(data=dfrDatpp,mapping=ggplot2::aes(x=z,y=lnR,colour=y,size=n))+
        ggplot2::geom_point(position=ggplot2::position_jitter(width=1)) +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  # print(p2);
  return(list(dfrDat=dfrDatpp,ps=list(p1=p1,p2=p2)));
}

#--check for outliers in area swept ratios----
#----(should be ~6 with 5 min. BSFRF tow, 30 min. NMFS tow)
ASRs = c(3,6,12);
lst1 = check_ASRs(dfrDat,ASRs);
lst1$ASRs = ASRs;

#----drop cells with < n_min individuals
n_min = 5;
lst2 = selectSizeData(lst1$dfrDat,n_min=n_min);
lst2$n_min = n_min;

#--save results----
lst3 = list(dfrInput=dfrDat,
            lstTrimmedASRs=lst1,
            lstTrimmedFinal=lst2
            );
wtsUtilities::saveObj(lst3,file.path(dirThs,"rda_Step2_TrimmedDataList.RData"));

