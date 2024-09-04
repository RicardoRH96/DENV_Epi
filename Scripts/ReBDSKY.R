#Plot BDSKY
library(coda); library(bdskytools); library(beastio)

bdskytrace <- readLog('~/bdsky_D2.log', burnin = 0.1)

Re_sky    <- beastio::getLogFileSubset(bdskytrace, "reproductiveNumber_BDSKY_Serial")
Re_hpd    <- t(beastio::getHPDMedian(Re_sky))
delta_hpd <- beastio::getHPDMedian(bdskytrace[, "becomeUninfectiousRate_BDSKY_Serial"])


#Plot
bdskytools::plotSkyline(1:10, Re_hpd, type='step', ylab="R")

#Smoothing
tmrca_med  <- median(bdskytrace[, "Tree.height"])
gridTimes  <- seq(0, median(tmrca_med), length.out=1040)  

Re_gridded <- mcmc(bdskytools::gridSkyline(Re_sky, bdskytrace[, "origin_BDSKY_Serial"], gridTimes))
Re_gridded_hpd <- t(getHPDMedian(Re_gridded))

#Plot smooth
times <- 2024.306 - gridTimes
plotSkyline(times, Re_gridded_hpd, xlab="Date", ylab="Re", type="smooth")


#Final plot
par(mar=c(5,4,4,4)+0.1)

plotSkylinePretty(range(times), as.matrix(delta_hpd), type='step', axispadding=0.0, 
                  col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), col.axis=pal.dark(cblue), 
                  ylab=expression(delta), side=4, yline=2, ylims=c(0,1), xaxis=FALSE)

plotSkylinePretty(times, Re_gridded_hpd, type='smooth', axispadding=0.0, 
                  col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange), 
                  xlab="Time", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=TRUE, 
                  ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,7), new=TRUE, add=TRUE)
