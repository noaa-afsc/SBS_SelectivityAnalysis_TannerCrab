#--calculate empirical selectivity from data using selfisher package

#require(tcsamSurveyData);
require(magrittr);
require(selfisher); require(bbmle);
source("CalcFitPredictAndPlot.R")

  #--retrieve proportions
  dfrPropsAllYears = wtsUtilities::getObj("dfrPropsAllYears.RData");

  #--apply models to proportions
  models<-list(); #--list, by sex, with all models and "best" model 
  # prdDFRs<-list();#--list, by sex, with "best" selectivity function
  pdf(file="EstimatedSelectivity_SelfisherModels.pdf",width=8,height=4.5);
  for (x in uXs){
    tmp<-dfrPropsAllYears[dfrPropsAllYears$SEX==x,c("SIZE","YEAR","HAULJOIN","BOTTOM_DEPTH","GEAR_TEMPERATURE","q","propNMFS","numTot")];
    #zs<-seq(min(tmp$SIZE),max(tmp$SIZE),delta/2);
    zs<-seq(minSize[x],maxSize[x],delta/2);
    prdDFR_ZY<-expand.grid(SIZE=zs,
                           YEAR=as.character(uYs),
                           BOTTOM_DEPTH=median(tmp$BOTTOM_DEPTH),
                           GEAR_TEMPERATURE=median(tmp$GEAR_TEMPERATURE),
                           HAULJOIN="1",q=1,numTot=1,plot="ZY",stringsAsFactors=FALSE);
    prdDFR_Z <-expand.grid(SIZE=zs,
                           YEAR="--",
                           BOTTOM_DEPTH=median(tmp$BOTTOM_DEPTH),
                           GEAR_TEMPERATURE=median(tmp$GEAR_TEMPERATURE),
                           HAULJOIN="1",q=1,numTot=1,plot="Z",stringsAsFactors=FALSE);
    prdDFR_ZY_RE<-rbind(prdDFR_ZY,prdDFR_Z); prdDFR_ZY_RE$plot<-"ZY_RE";
    prdDFR_D    <-expand.grid(SIZE=center[x],
                           YEAR="--",
                           BOTTOM_DEPTH=seq(min(tmp$BOTTOM_DEPTH),max(tmp$BOTTOM_DEPTH),2.5),
                           GEAR_TEMPERATURE=median(tmp$GEAR_TEMPERATURE),
                           HAULJOIN="1",q=1,numTot=1,plot="D",stringsAsFactors=FALSE);
    prdDFR_ZD<-rbind(prdDFR_Z,prdDFR_D);
    #tmp$SCALED_SIZE   <-as.numeric((tmp$SIZE   -minSize[x])/width[x]);
    #prdDFR$SCALED_SIZE<-as.numeric((prdDFR$SIZE-minSize[x])/width[x]);
    #--find "best" model collapsing data over all years
    mdls<-list();
    #--fit a single "global" selectivity with fixed factors for slope and intercept
    str<-paste0("Global size-based logistic function for ",x,"s");
    cat("\n\n",str,"\n"); writeTextToPDF(str);
    res<-doAllSelfisher(rformula=propNMFS~SIZE,
                        input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR_Z,type_predict="selection",plotFit=FALSE);
    mdls[[paste0("mLF0_",0)]]<-res;

    # # --fit annual logistic selectivity curves with different slopes and intercepts as fixed factors
    # res<-doAllSelfisher(rformula=propNMFS~0+YEAR*SIZE,
    #                     input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR_ZY,type_predict="selection",plotFit=FALSE);
    # mdls[[paste0("mLF1_",1)]]<-res;

    # #--fit a single logistic selectivity function with random effects for annual variability in intercept and slope
    # res<-doAllSelfisher(rformula=propNMFS~1+SIZE+(1+SIZE|YEAR),
    #                     input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR_ZY_RE,type_predict="selection",plotFit=FALSE);
    # mdls[[paste0("mLR0_",0)]]<-res;

    #--fit a single, "global" size-based spline function, with different degrees of freedom
    for (df in 3:10){
      str<-paste0("Global size-based spline function for ",x,"s, df =",df);
      cat("\n\n",str,"\n"); writeTextToPDF(str);
      res<-doAllSelfisher(rformula=propNMFS~0+splines::ns(SIZE,df=(df),intercept=TRUE),
                          input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR_Z,type_predict="selection",plotFit=FALSE);
      mdls[[paste0("mFSZ_",df)]]<-res;
      # for (dfd in 3:5){
      #   res<-doAllSelfisher(rformula=propNMFS~0+splines::ns(SIZE,df=3,intercept=TRUE)+splines::ns(BOTTOM_DEPTH,df=3,intercept=TRUE),
      #                       input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR_ZD,type_predict="selection",plotFit=FALSE);
      #   mdls[[paste0("mFSZD_",df,"_",dfd)]]<-res;
      # }
    }

   #  #--fit annual size-based spline functions, with different degrees of freedom
   #  for (df in 3:10){
   #    res<-doAllSelfisher(rformula=propNMFS~0+splines::ns(SIZE,df=(df),intercept=TRUE):YEAR,
   #                        input_data=tmp,Lp="basic",printSummary=TRUE,predict=TRUE,new_data=prdDFR,type_predict="selection",plotFit=FALSE);
   #    mdls[[paste0("mFSZY_",df)]]<-res;
   #    for (dfd in 3:5){
   #      res<-doAllSelfisher(rformula=propNMFS~0+splines::ns(SIZE,df=(df),intercept=TRUE):YEAR+splines::ns(BOTTOM_DEPTH,df=(dfd),intercept=FALSE),
   #                          input_data=tmp,Lp="basic",printSummary=TRUE,predict=FALSE,new_data=prdDFR,type_predict="selection",plotFit=FALSE);
   #      mdls[[paste0("mFSZYD_",df,"_",dfd)]]<-res;
   #    }
   # }

    # for (df in 2:10){
    #   # mdls[[paste0("ma0_",df)]]<-selfisher(data=tmp,rformula=propNMFS~-1+offset(log(q))+ns(SIZE,df)                                            ,total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("ma1_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+ns(SIZE,df)+ns(BOTTOM_DEPTH,df)                        ,total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("ma2_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+ns(SIZE,df)                    +ns(GEAR_TEMPERATURE,df),total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("ma3_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+ns(SIZE,df)+ns(BOTTOM_DEPTH,df)+ns(GEAR_TEMPERATURE,df),total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("m10_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+YEAR:ns(SIZE,df)                                            ,total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("m11_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+YEAR:ns(SIZE,df)+ns(BOTTOM_DEPTH,df)                        ,total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("m12_",df)]]<-selfisher(data=tmp,rformula=propNMFS~offset(log(q))+YEAR:ns(SIZE,df)                    +ns(GEAR_TEMPERATURE,df),total=numTot,haul=HAULJOIN);
    #   # mdls[[paste0("m13_",df)]]<-selfisher(data=tmp,rformula=propNMFS~-1+offset(log(q))+YEAR:ns(SIZE,df)+ns(BOTTOM_DEPTH,df)+ns(GEAR_TEMPERATURE,df),total=numTot,haul=HAULJOIN);
    # }
    fits<-list(); for (mdl in names(mdls)) fits[[mdl]]<-mdls[[mdl]]$fit;
    str<-paste0("Best model for ",x,"s");
    cat("\n\n",str,"\n"); writeTextToPDF(str);
    tbl<-BICtab(fits); txt<-print(tbl); writeTextToPDF(txt,,cex=1);
    best_name<-attr(tbl,"row.names",exact=TRUE)[1];
    str<-paste0("Best model for ",x,"s is ",best_name);
    cat("\n\n",str,"\n"); writeTextToPDF(str);
    df<-as.numeric(stringr::str_sub(best_name,6,-1));
    best_model<-mdls[[best_name]];
    mdls[["best"]]<-best_model;
    print(best_model$predictedPlot);
    #   # prdDFRA<-expand.grid(SIZE=unique(tmp$SIZE),q=1,numTot=1,HAULJOIN=1,stringsAsFactors=FALSE);
    #   # prdDFRY<-expand.grid(SIZE=unique(tmp$SIZE),YEAR=as.character(uYs),q=1,numTot=1,HAULJOIN=1,stringsAsFactors=FALSE);
    # for (i in 1:145){
    #   cat(i,"\n")
    #   prdSel<-predict(best_model,newdata=prdDFR,se.fit=TRUE,type="selection",na.action=na.omit);
    # }
    #   # prdSel<-predict(bestA,newdata=prdDFRA,se.fit=TRUE,type="selection");
    #   # prdSel<-predict(bestY,newdata=prdDFRY,se.fit=TRUE,type="selection");
    # prdDFR$est_sel<-prdSel$fit;
    # prdDFR$se_sel<-prdSel$se.fit
    # prdRat<-predict(best_model,newdata=prdDFR,se.fit=TRUE,type="ratio");
    # prdDFR$est_ratio<-prdRat$fit;
    # prdDFR$se_ratio <-prdRat$se.fit
    # prdRsp<-predict(best_model,newdata=prdDFR,se.fit=TRUE,type="response");
    # prdDFR$est_resp<-prdRsp$fit;
    # prdDFR$se_resp <-prdRsp$se.fit
    # p <- ggplot(data=prdDFR,mapping=aes(SIZE,y=est_sel,colour=YEAR,group=YEAR)) + geom_line();
    # p <- p + geom_ribbon(aes(ymin=est_sel-se_sel,ymax=est_sel+se_sel),alpha=0.2);
    # p <- p + labs(x="size (mm CW)",y="estimated selectivity",colour="study year",
    #               subtitle="by year")
    # print(p);
    # p <- ggplot(data=prdDFR,mapping=aes(SIZE,y=est_ratio,colour=YEAR,group=YEAR)) + geom_line();
    # p <- p + geom_ribbon(aes(ymin=est_ratio-se_ratio,ymax=est_ratio+se_ratio),alpha=0.2)
    # p <- p + labs(x="size (mm CW)",y="estimated ratio",colour="study year",
    #               subtitle="by year")
    # print(p);
    # p <- ggplot(data=prdDFR,mapping=aes(SIZE,y=est_resp,colour=YEAR,group=YEAR)) + geom_line();
    # p <- p + geom_ribbon(aes(ymin=est_resp-se_resp,ymax=est_resp+se_resp),alpha=0.2)
    # p <- p + labs(x="size (mm CW)",y="estimated response",colour="study year",
    #               subtitle="by year")
    # print(p);
    models[[x]]<-mdls;
    # prdDFRs[[x]]<-prdDFR;
  }#--x
#}#--iB
dev.off();

for (x in uXs){
  tmp<-models[[x]]$best$predicted %>% subset(27.5<=SIZE);
  nr<-nrow(tmp);
  irs<-seq(1,nr,2);
  tmp <- tmp[irs,];
  readr::write_csv(tmp,path=paste0("PredictedSelectivity_",x,".csv"));
}

wtsUtilities.saveObj(models,file="EstimatedSelectivityFromData_SelfisherResults.RData");
