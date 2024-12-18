#' 
#' @param dfr - dataframe with k-folding model statistics diff'd to base model
#' @param dfrd - dataframe with means of statistics in dfr
#' @param types - vector of types of statistics to plot 
#' @return ggplot object
plotStats_BestModels<-function(dfr,dfrd,
                               types=c("impr_mspe_prd",
                                       "impr_aic")){
  # dfr  = dfrCrsVald;
  # dfrd = dfrCrsValdp1;
  fcn<-function(x){
    xp = x |> 
             stringr::str_replace("impr_","") |> 
             stringr::str_replace("_prd","");
    return(xp);
  }
  typesp = toupper(types |> fcn());
  idx  = which(dfrd$i==1);
  dfrd = dfrd |> dplyr::slice_head(n=idx);
  dfrp = dfr |> 
           dplyr::filter(i %in% dfrd$i) |> 
           tidyr::pivot_longer(tidyselect::starts_with("impr"),
                               names_to="type",values_to="improvement") |> 
           dplyr::filter(type %in% types) |>
           dplyr::mutate(smths=factor(smths,levels=dfrd$smths),
                         type=factor(toupper(fcn(type)),levels=typesp));
  dfrdp = dfrd |> 
            dplyr::select(!tidyr::ends_with("impr_mase_prd")) |> 
            tidyr::pivot_longer(tidyselect::starts_with("mn_impr"),
                                names_to="type",values_to="improvement") |> 
            dplyr::filter(type %in% paste0("mn_",types)) |>
            dplyr::mutate(type=gsub("mn_","",type),
                          smths=factor(smths,levels=dfrd$smths),
                          type=factor(toupper(fcn(type)),levels=typesp));
  p = ggplot(dfrp,aes(x=smths,y=improvement,colour=type)) + 
        geom_boxplot() + 
        geom_point()+
        geom_point(data=dfrdp,shape=22,colour="black")+
        geom_hline(yintercept=0,linetype=3) + 
        facet_wrap(~type,nrow=1,scales="free_y") + 
        labs(x="model",y="improvement") + 
        wtsPlots::getStdTheme() + 
        theme(legend.position="none",
              axis.text.x = element_text(angle=45,hjust=1,size=12),
              axis.title.x=element_blank());
      return(p);
}
#--plotStats_BestModels(dfrCrsVald,dfrCrsValdp1);
