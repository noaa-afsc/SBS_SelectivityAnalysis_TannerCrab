dfrDatRE = dfrDat %>% subset(n>5);
dfrDatRE$h = as.factor(dfrDatRE$h);

frmla = p~ti(z,bs="ts") + s(h,bs="re");
#frmla = p~ti(z,bs="ts");
mdl = gam(data=dfrDatRE,family=fam,formula=frmla,weights=n,select=FALSE,scale=0,offset=lnq,method="ML");
dfrDatRE %<>% dplyr::mutate(fits=fitted(mdl),                         #--fitted values on response scale (i.e., prNMFS)
                            residuals=residuals(mdl,type="deviance"), #--deviance residuals for qqplot
                            lnR=log(fits/(1-fits))-lnq,               #--corresponding values of lnR
                            R=exp(lnR));                              #--corresponding values of R

dfrPrd  = predict(mdl,type="iterms",se.fit=TRUE);
dfrPrdp = tibble::tibble(pred=predict(mdl,type="link"),
                         int=mdl$coefficients[1]) %>%
          dplyr::mutate(`ti(z)`=dfrPrd$fit[,1],
                        `se[ti(z)]`=dfrPrd$se.fit[,1],
                        `lower[ti(z)]`=`ti(z)`-`se[ti(z)]`,
                        `upper[ti(z)]`=`ti(z)`+`se[ti(z)]`,
                        pred1 = `ti(z)`+`int`,
                        `lower(pred1)`=`lower[ti(z)]`+`int`,
                        `upper(pred1)`=`upper[ti(z)]`+`int`) %>%
          dplyr::mutate(`R[ti(z)]`=exp(`ti(z)`+`int`),
                        `lower(R[ti(z)])`=exp(`ti(z)`+`int`-`se[ti(z)]`),
                        `upper(R[ti(z)])`=exp(`ti(z)`+`int`+`se[ti(z)]`));
dfrPred = cbind(dfrDatRE,dfrPrdp);
dfrR = dfrPred %>% dplyr::select(z,pred1,`lower(pred1)`,`upper(pred1)`,
                                 `R[ti(z)]`,`lower(R[ti(z)])`,`upper(R[ti(z)])`) %>% 
                   dplyr::distinct_all();

ggplot2::ggplot(data=dfrPred,mapping=aes(x=z,y=log(p/(1-p))-lnq,size=n)) + 
  ggplot2::geom_point() +scale_size_area() +
  ggplot2::geom_line(data=dfrPred,mapping=ggplot2::aes(x=z,y=pred1,size=NULL),color="blue")+
  ggplot2::geom_ribbon(data=dfrPred,mapping=ggplot2::aes(x=z,ymin=`lower(pred1)`,ymax=`upper(pred1)`,size=NULL),color="blue",alpha=0.2);
ggplot2::ggplot(data=dfrPred,mapping=ggplot2::aes(x=z,y=exp(pred),group=h)) + 
  ggplot2::geom_point(alpha=0.5) +
  ggplot2::geom_line(alpha=0.3) +
  ggplot2::geom_line(data=dfrR,mapping=ggplot2::aes(x=z,y=`R[ti(z)]`,size=NULL,group=NULL),color="blue",size=1) +
  ggplot2::geom_ribbon(data=dfrR,mapping=ggplot2::aes(x=z,y=NULL,ymin=`lower(R[ti(z)])`,ymax=`upper(R[ti(z)])`,size=NULL,group=NULL),color="blue",alpha=0.2)
# ggplot2::ggplot(data=dfrDatRE,mapping=aes(x=p,y=fits,size=n)) + 
#   ggplot2::geom_point() +scale_size_area();

summary(mdl);
summary(mdlsB[[1]]$model)    
