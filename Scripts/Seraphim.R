#install
install_github("sdellicour/seraphim/unix_OS")
#load libraries
library(seraphim)


treeDirectory <- '~/f12_extTrees/'
allTrees <- scan(file = '~/lineage_F12_aln.trees', what = '', sep = '\n', quiet = TRUE)
burnIn =60
randomSampling = FALSE
numberoftreestosample = 100
mostRecentSamplingDate = 2024.306
coordinateAttributeName = 'location'


treeExtractions(treeDirectory,allTrees,burnIn,randomSampling,numberoftreestosample,
                mostRecentSamplingDate,coordinateAttributeName)

#Get dispersal velocities
spreadStatistics(treeDirectory, nberOfExtractionFiles = 100, timeSlices = 100, onlyTipBranches = FALSE,
                 showingPlots = TRUE, outputName = 'DENV2', nberOfCores = 4, slidingWindow = 1)

#Map
library(diagram)
source('~/mccExtractions.r')


mcc_tre = readAnnotatedNexus('~/lineage_f12.annotated.mcc.tre')

mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDate)

nberOfExtractionFiles = numberoftreestosample
prob = 0.95
precision = 0.025
startDatum = min(mcc_tab[,'startYear'])


polygons = suppressWarnings(spreadGraphic2(treeDirectory,
                                           nberOfExtractionFiles, prob, startDatum, precision))


colour_scale = colorRampPalette(brewer.pal(11,"RdYlGn"))(141)[21:121]
minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons)) {
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"40")
}


geodata = sf::read_sf('~/Downloads/gadm41_COL_shp/gadm41_COL_1.shp')
template_raster = raster(geodata$geometry)

borders = crop(getData("GADM", country="COL", level=1), extent(template_raster))


#plot
dev.new(width=6, height=6.3)
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(geodata$geometry, col="white", box=F, axes=F, colNA="grey90", legend=F)
for (i in 1:length(polygons)) {
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
plot(borders, add=T, lwd=0.1, border="gray10")
for (i in 1:dim(mcc_tab)[1]) {
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]),
              cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA,
              arr.pos=F, curve=0.1, dr=NA, endhead=F)
}
for (i in dim(mcc_tab)[1]:1) {
  if (i == 1) {
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16,
           col=colour_scale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1,
           col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16,
         col=endYears_colours[i], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1,
         col="gray10", cex=0.8)
}
rect(xmin(geodata$geometry), ymin(geodata$geometry), xmax(geodata$geometry),
     ymax(geodata$geometry), xpd=T, lwd=0.2)
axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))),
     pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))),
     pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=1, tck=-0.01, col.axis="gray30")
rast = raster(matrix(nrow=1, ncol=2))
rast[1] = min(mcc_tab[,"startYear"])
rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5,
     legend.shrink=0.3, smallplot=c(0.40,0.80,0.14,0.155), legend.args=list(text="",
                                                                            cex=0.7, line=0.3, col="gray30"), horizontal=T, axis.args=list(cex.axis=0.6,
                                                                                                                                           lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0),
                                                                                                                                             at=seq(2016.4,2017.2,0.2)))
