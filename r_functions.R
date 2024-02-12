
lstPredTermTbls = predSmoothTerms(model,grids);

plotSmoothTerm1D("ti(z)",lstPredTermTbls,dfrDat,xlab="size (mm CW)")

plotSmoothTerm2D("ti(z,d)",lstPredTermTbls,dfrDat,xlab="size (mm CW)",ylab="depth (m)")

ps = plotSmoothTerms(lstPredTermTbls,dfrDat,labs=lbls);

