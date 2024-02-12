#--calculate numbers caught and numbers sampled between survey types
dirPrj = rstudioapi::getActiveProject();
dirThs = dirname(rstudioapi::getActiveDocumentContext()$path);

#--
cutpts = seq(-0.5,204.5,5);

#--calculate numbers sampled by sex/5-mm size bin
#----NOT NECESSARILY numbers caught, due to subsampling
#----estimate actual numbers caught by rescaling numbers sampled by the (assumed) sampling factor
lst = wtsUtilities::getObj(file.path(dirThs,"rda_Step1_SBS_RawData.RData"));
dfrID_SBS = dplyr::bind_rows(
              lst$dfrID_BSFRF_SBS |> 
                dplyr::inner_join(lst$dfrHD_BSFRF_SBS,by="HAULJOIN") |> 
                dplyr::mutate(bin=round(cutpts[cut(SIZE,cutpts,labels=FALSE)])+2.5) |> 
                dplyr::select(YEAR,GIS_STATION,SEX,MATURITY,SHELL_CONDITION,SIZE,bin,
                              numIndivs,SAMPLING_FACTOR) |> 
                dplyr::mutate(survey="BSFRF"),
              lst$dfrID_NMFS_SBS |> 
                dplyr::inner_join(lst$dfrHD_NMFS_SBS,by="HAULJOIN") |> 
                dplyr::mutate(bin=round(cutpts[cut(SIZE,cutpts,labels=FALSE)])+2.5) |> 
                dplyr::select(YEAR,GIS_STATION,SEX,MATURITY,SHELL_CONDITION,SIZE,bin,
                              numIndivs,SAMPLING_FACTOR) |> 
                dplyr::mutate(survey="NMFS")); 
dfrNumByZB = dfrID_SBS |> 
      dplyr::group_by(survey,YEAR,SEX,bin) |> 
      dplyr::summarize(sampled=wtsUtilities::Sum(numIndivs),
                       caught=wtsUtilities::Sum(SAMPLING_FACTOR)) |> 
      dplyr::ungroup() |> 
      tidyr::pivot_longer(c("sampled","caught"),names_to="type",values_to="num");

wtsUtilities::saveObj(dfrNumByZB,file.path(dirThs,"rda_Step2_SBS_Nums.RData"));

#--plot results
require(ggplot2);
require(rlang);
plotNums<-function(sex,col,label){
  p = ggplot(dfrNumByZB |> dplyr::filter(SEX==sex),
             aes(x=bin,y={{col}},colour=survey,linetype=type)) + 
        geom_line() + 
        geom_vline(xintercept=25,linetype=3) + 
        facet_wrap(~YEAR,nrow=2,scales="free_y") + 
        labs(x="size (mm CW)",y=label) + 
        theme(legend.position=c(0.99,0.48),
              legend.justification=c(1,1)) +
        wtsPlots::getStdTheme();
  return(p);
}
plotNums("MALE",num,"males");
plotNums("FEMALE",num,"females");

