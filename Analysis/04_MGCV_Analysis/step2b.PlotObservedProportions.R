#plot observed proportions from step2a
#--load required packages
require(ggplot2);
#--other required packages
#----wtsUtilities

#--get previously-created results
dfrPropsAllYears   = wtsUtilities::getObj("dfrPropsAllYears.RData");
dfrMnPropsAllYears = wtsUtilities::getObj("dfrMnPropsAllYears.RData");

#--plot results
doPDFs = FALSE;
uXs = c("MALE","FEMALE");
for (x in uXs){
  if (doPDFs) pdf(file=paste0("plots2a1.ObservedProportions_",x,".pdf"),width=9,height=6.5);
  tmp0<-dfrPropsAllYears[dfrPropsAllYears$SEX==x,];
  tmp1<-dfrMnPropsAllYears[dfrMnPropsAllYears$SEX==x,];
  p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=propNMFS,size=numTot,colour=YEAR));
  p <- p + geom_point(alpha=0.2,stat="identity",position=position_jitter(width=2.0)) + scale_size_area();
  p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=propNMFS,colour=YEAR),alpha=0.4,size=2);
  p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=mnPropNMFS,colour=YEAR),alpha=0.4,size=1);
  p <- p + facet_grid(rows=YEAR~.);
  p <- p + labs(x="size (mm CW)",y="proportion (NMFS/[NMFS+BSFRF])",colour="year",size="total\nnumber");
  print(p);
  p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=lnR,size=numTot,colour=YEAR));
  p <- p + geom_point(alpha=0.2,stat="identity",position=position_jitter(width=2.0)) + scale_size_area();
  p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=lnR,colour=YEAR),alpha=0.4,size=2);
#  p <- p + geom_line(data=tmp1,mapping=aes(x=SIZE,y=mnLnR,colour=YEAR),alpha=0.4,size=1);
  p <- p + facet_grid(rows=YEAR~.);
  p <- p + labs(x="size (mm CW)",y="ln(R) = logit(NMFS/[NMFS+BSFRF])-log(q)",colour="year",size="total\nnumber");
  print(p);
  if (doPDFs) dev.off();
}
# if (doPDFs) pdf(file=paste0("plots2a2.LnQ.pdf"),width=9,height=4.5);
# p <- ggplot(data=tmp0,mapping=aes(x=SIZE,y=log(q),size=numTot,colour=YEAR));
# p <- p + geom_point(alpha=0.4,stat="identity",position=position_jitter(width=0.5)) + scale_size_area();
# p <- p + labs(x="size (mm CW)",y="log(q)",colour="sex",size="total/number",shape="year");
# print(p);
# if (doPDFs) dev.off();
rm(p,tmp0,tmp1,x);
