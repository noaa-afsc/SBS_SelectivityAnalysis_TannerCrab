Binomial.gam.model4.new<-function (){
    dat <- pcapture.data.selectmodel.males
    width <- dat$width;    #--carapace widths
    meas.b <- dat$meas.b   #--number of measaured crab, BSFRF
    meas.n <- dat$meas.n   #--number of measured crab, NMFS
    sed <- dat$phi         #--sediment size, expressed in units of phi = -ln(grain size in mm)
    dep <- dat$depth       #--trawl depth
    #--fit to binomial model (successes, failures) with mean response s(Z) + s(phi, D))
    #----note: model includes an intercept term
    fit1 <- gam(cbind(meas.n, meas.b) ~ s(width) + s(sed, dep), family = binomial)
    print(summary(fit1))
    cat("AIC=", AIC(fit1), "\n")
    cat.b <- dat$num.b; #--expanded catch numbers, BSFRF
    cat.n <- dat$num.n; #--expanded catch numbers, NMFS
    cat.b[cat.b == 0] <- 1
    cat.n[cat.n == 0] <- 1
    sp.b <- meas.b/cat.b; #--sampling factor, BSFRF
    sp.n <- meas.n/cat.n; #--sampling factor, NMFS
    sp.b[sp.b == 0] <- 1
    sp.n[sp.n == 0] <- 1
    sam.rat <- sp.b/sp.n; #--sampling ratio
    area.b <- dat$area.b;  #--area swept, BSFRF
    area.n <- dat$area.n;  #--area swept, NMFS
    a.rat <- area.b/area.n;  #--relative area swept, BSFRF:NMFS
    #--predict response
    ypred <- predict(fit1, type = "response")
    ypredx <- tapply(ypred, width, mean)
    y <- ypredx
    widths <- as.numeric(names(ypredx))
    plot(widths, ypredx, ylim = c(0, 1), las = 1, type="l", lty=1, lwd=2);
    points(width,meas.n/(meas.n+meas.b));
    title("predicted response");
    cat(widths, "\n")
    #--what's this section doing?
    a.rat.mean <- mean(a.rat); #--was 0.144
    sam.rat <- (sam.rat)^0.25; #--variance stabilization?
    gfit <- gam(sam.rat ~ s(width))
    s.pred <- predict.gam(gfit, newdata = data.frame(width = widths))
    s.pred <- s.pred^4
    plot(widths, s.pred, type="l", lty=1, lwd=2);
    points(width,sam.rat);
    title("sampling ratio (BSFRF/NMFS)");

    #--plot capture probability (similar to Fig. 6 from Somerton et al., 2013)
    p <- (y * s.pred * a.rat.mean)/(1 - y); #--follows from eq. 3 in Somerton et al 2013
    plot(widths, p, ylim = c(0, 1), las = 1, type = "l")
    title("estimated capture probability (r(width))")

    #--can't do the following: "loglike.sol" is unavailable from what Bob sent me
    #----what is it?--assessment model estimate of selectivity?
    if (FALSE){
        points(loglike.sol$width, loglike.sol$pcap, col = 5)
        fit3 <- gam(loglike.sol$pcap ~ s(loglike.sol$width))
        lines(loglike.sol$width, fit3$fitted, col = 5)
    }
    return(fit1)
}

#--run function
require(mgcv);
pcapture.data.selectmodel.males<-read.csv("McConnaughey_data.csv",stringsAsFactors=FALSE,header=TRUE,skip=1);
predict.gam.males.new<-Binomial.gam.model4.new();#--fitted model for gam(cbind(meas.n, meas.b) ~ s(width) + s(sed, dep),family = binomial)

