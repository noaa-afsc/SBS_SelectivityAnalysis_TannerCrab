doAllSelfisher<-function(rformula,
                         input_data,
                         start=NULL,
                         Lp="basic",
                         printSummary=TRUE,
                         predict=TRUE,
                         new_data=NULL,
                         type_predict=c("response","selection","prob", "ratio", "link"),
                         plotFit=FALSE,
                         showFitPlot=TRUE,
                         type_residuals=c("deviance","pearson","response"),
                         plotResiduals=TRUE,
                         showResidualsPlot=TRUE,
                         plotQQ=TRUE,
                         showQQPlot=TRUE,
                         plotPredicted=TRUE,
                         showPredictedPlot=TRUE,
                         debug=FALSE){
    fit<-fitSelfisher(rformula=rformula,data=input_data,start=start,Lp=Lp,debug=debug);
    if (printSummary) print(summary(fit));
    if (debug) return(fit);
    resPlot<-fitPlot<-qqPlot<-prd<-prdPlot<-NULL;
    #--plot fit
    if (plotFit) fitPlot<-plotFit(fit,input_data,showPlot=showFitPlot);
    #--plot residuals
    if (plotResiduals) resPlot<-plotResiduals(fit,input_data,type=type_residuals[1],showPlot=showResidualsPlot);
    #--plot residuals
    if (plotQQ) qqPlot<-plotQQ(fit,input_data,type=type_residuals[1],showPlot=showQQPlot);
    #--predict results for new data
    if (predict) {
        prd<-predictSelfisher(fit,new_data,type=type_predict[1],showPlot=FALSE);
        #--plot predicted
        if (plotPredicted) prdPlot<-plotSelfisher(prd,type=type_predict[1],showPlot=showPredictedPlot);
    }
    return(list(fit=fit,predicted=prd,fitPlot=fitPlot,resPlot=resPlot,qqPlot=qqPlot,predictedPlot=prdPlot));
}

fitSelfisher<-function(rformula,
                       data,
                       start=NULL,
                       Lp="basic",
                       debug=FALSE){
    res<-selfisher(data=data,rformula=rformula,total=numTot,qratio=q,haul=HAULJOIN,start=start,Lp=Lp,debug=debug);
    return(res)
}

predictSelfisher<-function(model,
                          newdata,
                          type=c("response","selection","prob", "ratio", "link"),
                          showPlot=TRUE){
    type=type[1];
    nr<-nrow(newdata);
    prdSel<-predict(model,newdata=newdata,se.fit=TRUE,type=type,na.action=na.omit);
    newdata$est<-prdSel$fit[1:nr];
    newdata$se <-prdSel$se.fit[1:nr];
    if (showPlot) plotSelfisher(newdata,type=type,showPlot=TRUE);
    return(newdata);
}

plotSelfisher<-function(newdata,
                          type=c("response","selection","prob", "ratio", "link"),
                          plotSEs=TRUE,
                          showPlot=TRUE){
    type<-type[1];
    if (type=="response")  ylab="estimated response";
    if (type=="selection") ylab="estimated selectivity";
    if (type=="prob")      ylab="estimated p";
    if (type=="ratio")     ylab="estimated ratio (r/[1-r])";
    if (type=="link")      ylab="estimated link";
    uPs<-unique(newdata$plot);
    plots<-list();
    for (uP in uPs){
        tmp0<-newdata[newdata$plot==uP,];
        tmp<-tmp0; tmp1<-NULL;
        uYs<-unique(tmp$YEAR);
        if ((length(uYs)>1)&("--" %in% uYs)) {idx<-tmp$YEAR=="--"; tmp1<-tmp[idx,]; tmp<-tmp[!idx,];}
        if (grepl("Z",uP,fixed=TRUE)) {x<-"SIZE"; xlab<-"size (mm CW)";}
        if (grepl("D",uP,fixed=TRUE)) {x<-"BOTTOM_DEPTH"; xlab<-"bottom depth (m)";}
        if (grepl("T",uP,fixed=TRUE)) {x<-"GEAR_TEMPERATURE"; xlab<-"gear temperature (deg C)";}
        p <- ggplot2::ggplot(data=tmp,mapping=ggplot2::aes_string(x=x,y="est",colour="YEAR",group="YEAR"));
        p <- p + ggplot2::geom_line();
        if (plotSEs) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=est-se,ymax=est+se),alpha=0.3);
        if (!is.null(tmp1)) {
            p <- p + ggplot2::geom_line(data=tmp1,colour="green");
            if (plotSEs) p <- p + ggplot2::geom_ribbon(data=tmp1,mapping=ggplot2::aes(x=as.name(x),y=est,ymin=est-se,ymax=est+se),alpha=0.2,colour="green",fill="green");
        }
        p <- p + ggplot2::labs(x=xlab,y=ylab,colour="study\nyear");
        if (type=="selection") p <- p + ggplot2::ylim(c(0,1));
        if (showPlot) print(p);
        plots[[uP]]<-p;
    }
    return(plots)
}

plotFit<-function(model,
                  input_data,
                  showPlot=TRUE){
    nr<-nrow(input_data);
    prd<-predict(model,newdata=input_data,se.fit=TRUE,type="response",na.action=na.omit);
    input_data$logit_resp<-log(input_data$propNMFS/(1-input_data$propNMFS));
    input_data$est<-prd$fit[1:nr];
    input_data$se <-prd$se.fit[1:nr];
    p <- ggplot2::ggplot(data=input_data,mapping=ggplot2::aes(x=SIZE,y=logit_resp,size=numTot));
    p <- p + ggplot2::geom_point(colour="black",alpha=0.2,stat="identity",position=ggplot2::position_jitter(width=0.5)) + ggplot2::scale_size_area();
    p <- p + ggplot2::geom_errorbar(ggplot2::aes(ymin=est-1.96*se,ymax=est+1.96*se),alpha=0.5,colour="red");
    p <- p + ggplot2::geom_point(ggplot2::aes(y=est),alpha=0.7,shape=17,size=1,colour="red");
    p <- p + ggplot2::facet_grid(rows=vars(YEAR))
    p <- p + ggplot2::labs(x="size (mm CW)",y="logit(NMFS/[NMFS+BSFRF])",size="total\nnumber");
    if (showPlot) print(p);
    return(p)
}
plotResiduals<-function(model,
                        input_data,
                        type=c("deviance","pearson","response"),
                        showPlot=TRUE){
    type<-type[1];
    nr<-nrow(input_data);
    resids<-residuals(model,type=type);
    input_data$residuals<-as.numeric(resids)[1:nr];
    p <- ggplot2::ggplot(data=input_data,mapping=ggplot2::aes(x=SIZE,y=residuals,size=numTot));
    p <- p + ggplot2::geom_point(alpha=0.5,stat="identity",position=ggplot2::position_jitter(width=0.5)) + ggplot2::scale_size_area();
    p <- p + ggplot2::facet_grid(rows=ggplot2::vars(YEAR))
    p <- p + ggplot2::labs(x="size (mm CW)",y="deviance reiduals",size="total\nnumber");
    if (showPlot) print(p);
    return(p)
}
plotQQ<-function(model,
                 input_data,
                 type=c("deviance","pearson","response"),
                 showPlot=TRUE){
    type<-type[1];
    nr<-nrow(input_data);
    resids<-residuals(model,type=type);
    input_data$residuals<-as.numeric(resids)[1:nr];
    #--use ggplot2 to calculate the theoretical and sample quantiles
    qq <- ggplot2::ggplot(data=input_data,mapping=ggplot2::aes(sample=residuals));
    qq <- qq + ggplot2::geom_qq(alpha=0.5) + ggplot2::geom_qq_line();
    #--extract dataframe with  the theoretical and sample quantiles
    dfr<-ggplot2::ggplot_build(qq)$data[[1]][,c("sample","theoretical")];
    dfr<-ggplot2::layer_data(qq,i=1)[,c("sample","theoretical")];
    #--add in auxilliary information about the residuals
    idx<-order(input_data$residuals);
    dfr$YEAR<-input_data$YEAR[idx];
    dfr$numTot<-input_data$numTot[idx];
    #--set y-axis label
    ylab<-"observed quantiles\n";
    if (type=="deviance") ylab<-paste0(ylab,"(deviance residuals)");
    if (type=="pearson")  ylab<-paste0(ylab,"(pearson residuals)");
    if (type=="response") ylab<-paste0(ylab,"(response residuals)");
    #--create plot
    p <- ggplot2::ggplot(data=dfr,mapping=ggplot2::aes(x=theoretical,y=sample,colour=YEAR,size=numTot));
    p <- p + ggplot2::geom_abline(slope=1,intercept=0);
    p <- p + ggplot2::geom_point(alpha=0.5) + ggplot2::scale_size_area();
    p <- p + ggplot2::labs(x="theoretical quantiles",y=ylab,size="total\nnumber",colour="year");
    if (showPlot) print(p);
    return(p)
}

getPredictorVariables<-function(model){
    return(attr(terms(model), "predvars"));
}

getSpline<-function(model,
                    predvar="BOTTOM_DEPTH"){
    predvars<-getPredictorVariables(model);
    for (i in 1:length(predvars)) {
        if (any(grepl(predvar,as.character(predvars[[i]]),fixed=TRUE))){
            return(predvars[[i]]);
        }
    }
    return(NULL)
}

predictSplineComponent<-function(model,
                                 predvar="BOTTOM_DEPTH"){
    scall<-getSpline(model,predvar);
    bndry<-scall$Boundary.knots;
    xvals<-seq(from=bndry[1],to=bndry[2],length.out=30);
    if (is.null(scall)) return(NULL);
    dfr<-data.frame(val=xvals); names(dfr)<-predvar;
    s<-eval(scall,envir=dfr);

    #--get estimated coefficient values
    fxdeffs<-selfisher::fixef(model);
    coeffs<-as.numeric(fxdeffs$r[grepl(predvar,names(fxdeffs$r),fixed=TRUE)])
    rss<-s %*% coeffs;
    rss<-1/(1+exp(-rss));#--assumes link function is logit

    #--plot
    plot(xvals,rss)
    return(data.frame(x=xvals,y=rss));
}

writeTextToPDF<-function(txt,cex=1,x=0.5,y=0.5){
    oldpar<-par(oma=c(0,0,0,0));
    on.exit(par(oldpar));
    plot.new();
    txtp<-txt;
    if (inherits(txt,"matrix")|inherits(txt,"data.frame")){
        rwn<-max(nchar(rownames(txt),type="width"));
        txtp<-paste(paste0(" ",paste(character(rwn),collapse=" ")),paste0(colnames(txt),collapse=" "),"\n");
        for (r in 1:nrow(txt)){
            rtxt<-paste(rownames(txt)[r],paste(as.character(txt[r,]),collapse=" "),"\n");
            txtp<-paste0(txtp,rtxt);     
        }#--r
    }#--data.frame
    cat(txtp)
    text(x=x,y=y,labels=txtp,cex=cex,offset=0);
}
