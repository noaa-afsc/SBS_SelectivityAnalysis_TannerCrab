#--NOTE: run McConnaughey_opimodel.R first
#
Survey.weighted.efficiency<-function (sex){
        if(sex==1) dat1<-survey2010.males.phi
        if(sex==2) dat1<-survey2010.females.phi
# get rid of the unmeasurable crabs
        dat1<-dat1[!is.na(dat1$width),]

#
# constrain width to be a max of 125 so that size matches the prediction functions
#
        dat1<-dat1[dat1$width>=20,]
        if(sex==1) dat1$width[dat1$width>=125]<-125
        if(sex==2) dat1$width[dat1$width>=75]<-75
        phi<-dat1$phi
        width<-dat1$width
        width=trunc((width-1)/5)*5 + 3
        area<-dat1$area
        samp<-dat1$sampfac
        station<-as.character(dat1$station)
        dep<-dat1$depth
#
        sed<-phi
        indat<-data.frame(width,sed,dep)
        if(sex==1) gam.object<-predict.gam.males.new     #--Dave didn't provide this!! Assuming it's the 'fit1' from Bob's function
        if(sex==2) gam.object<-predict.gam.females.new
        pcap<-predict(gam.object, newdata=indat)

# calculate the ratio of the station area to the swept area
        stat.area<-rep(400,length(station))
        dat3<-strata.definitions.new           #--missing
        stationx<-as.character(dat3$station)
        areax<-dat3$area
        for (i in 1:length(areax)){
           stat.area[station==stationx[i]]<-areax[i]
        }
        area.ratio<-stat.area/area             #--area swept expansion factor (multiplicative)
#
# calculate the weighted average cap prob weighting by subsampling expansion and area ratio
        wt<-samp*area.ratio
        num<-wt*pcap
        den<-wt
# sum both over each  width interval
        numx<-tapply(num,width,sum);  #--summing "num"erator by width bin
        denx<-tapply(den,width,sum);  #--summing "den"ominator by width bin
# calculate weighted prob capture
        pcaptured<-numx/denx;             #--weighted capture probability
        pcaptured[pcaptured < 0.0]<-0.0
        widths<-as.numeric(names(denx))
        pcaptured=pcaptured[widths > 20]
        widths=widths[widths >20]
        plot(widths, pcaptured, ylab="Proportion captured", ylim=c(0,1),
             xlim=c(20,125),xlab="Carapace width (mm)", type="l",
             lty=2, las=1, cex=3.0, cex.lab=1.2, cex.axis=1.2,lwd=2.0)
# plot the previous version
        if(sex==1) dat1=pcapture.select.males.pooled
        if(sex==2) dat1=pcapture.select.females.pooled
        width<-dat1$width
        sed=dat1$phi
        dep<-dat1$depth
        indat=data.frame(width,sed,dep)
        pcap<-predict(gam.object, newdata=indat)
        pcap.x=tapply(pcap,width,mean)
        widths<-as.numeric(names(pcap.x))

        lines(widths,pcap.x)
        invisible()
}

#--run function
sex<-1;
survey2010.males.phi<-read.csv("Somerton_survey2010.males.phi.csv",stringsAsFactors=FALSE);
pcapture.select.males.pooled<-read.csv("Someron_pcapture.select.males.pooled.csv",stringsAsFactors=FALSE)
Survey.weighted.efficiency(sex);

