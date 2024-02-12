#--function to check for outliers based on q
#'
#' @title Edit data for outlier q's (area swept and sampling ratios)
#' 
#' @description Function to edit data for outlier q's (area swept and sampling ratios).
#' 
#' @param dfrDat - tibble with proprtion data to edit
#' @param qs - 3-element vector with min q, nominal q, max q
#' 
#' @return list
#' 
#' @export
#' 
check_qs<-function(dfrDat,qs){
  #--create dataframe for q's defining range
  dfrQs = tibble::tibble(q=qs[1:3],
                         y=1000+0*qs[1:3],
                         lbl=as.factor(c("min","nominal","max")));
  #--plot histogram of q's
  p1 = ggplot2::ggplot()+
        ggplot2::geom_histogram(data=dfrDat,mapping=ggplot2::aes(x=q,fill=y),binwidth=1) +
        ggplot2::geom_vline(data=dfrQs,mapping=ggplot2::aes(xintercept=q),linetype=2) +
        ggplot2::geom_text(data=dfrQs,mapping=ggplot2::aes(x=q,y=y,label=lbl)) +
        labs(x="q = (A_NMSF/A_BSFRF)*(S_NMFS/S_BSFRF)",y="count",fill="year");
  print(p1);
  
  #--plot original data on logit scale
  nrw_orig = nrow(dfrDat);
  dfrDatp = rbind(dfrDat %>% subset(dplyr::between(q,dfrQs$q[1],dfrQs$q[3])) %>% mutate(cat=paste0("in [",dfrQs$q[1],",",dfrQs$q[3],"]")),
                  dfrDat %>% subset(q<dfrQs$q[1])                            %>% mutate(cat=paste0("<",dfrQs$q[1])),
                  dfrDat %>% subset(q>dfrQs$q[3])                            %>% mutate(cat=paste0(">",dfrQs$q[3])));
  dfrDatp$cat = factor(dfrDatp$cat,levels=c(paste0("<",dfrQs$q[1]),
                                            paste0("in [",dfrQs$q[1],",",dfrQs$q[3],"]"),
                                            paste0(">",dfrQs$q[3])));
  p2 = ggplot2::ggplot(data=dfrDatp,mapping=ggplot2::aes(x=z,y=lgtp,fill=cat,size=n))+
        ggplot2::geom_point(position=position_jitter(width=1),alpha=0.8,shape=21) +
        ggplot2::scale_colour_viridis_d(option="plasma") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        ggplot2::labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)",fill="q");
  print(p2);
  
  #--edit data based on q limits
  dfrDatp  = dfrDat %>% subset(dplyr::between(q,dfrQs$q[1],dfrQs$q[3]));
  dfrDrop  = dfrDat %>% anti_join(dfrDatp,by=NULL) %>% dplyr::select(h) %>% dplyr::distinct();
  nrw_edit = nrow(dfrDatp);
  p3 = ggplot2::ggplot()+
        ggplot2::geom_point(data=dfrDatp,mapping=ggplot2::aes(x=z,y=lgtp,colour=y),position=position_jitter(width=1)) +
        ggplot2::geom_point(data=dfrDat %>% anti_join(dfrDatp,by=NULL),mapping=ggplot2::aes(x=z,y=lgtp),colour="black") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  print(p3);
  
  return(list(dfrDat=dfrDatp,dfrDrop=dfrDrop,ps=list(p1=p1,p2=p2,p3=p3),nrw_orig=nrw_orig,nrw_edit=nrw_edit));
}

selectSizeData<-function(dfrDat,n_min){
  #--subset the data based on minimum number of individuals caught
  dfrDatpp = dfrDat %>% subset(n>=n_min);
  
  #--plot the data that have been subsetted
  p1 = ggplot2::ggplot(data=dfrDat,mapping=ggplot2::aes(x=z,y=lgtp,colour=y))+
        ggplot2::geom_point(position=position_jitter(width=1)) +
        ggplot2::geom_point(data=dfrDat %>% anti_join(dfrDatpp,by=NULL),mapping=ggplot2::aes(x=z,y=lgtp),
                            position=position_jitter(width=1),colour="black") +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  print(p1);
  
  #--plot data that will be fit
  p2 = ggplot2::ggplot(data=dfrDatpp,mapping=ggplot2::aes(x=z,y=lgtp,colour=y,size=n))+
        ggplot2::geom_point(position=position_jitter(width=1)) +
        ggplot2::facet_grid(x~.,scales="free_y")+
        labs(x="size (mm CW)",y="log(R) = logit(phi)-ln(q)");
  print(p2);
  return(list(dfrDat=dfrDatpp,ps=list(p1=p1,p2=p2)));
}

