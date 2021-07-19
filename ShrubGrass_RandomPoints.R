# ShrubGrass_RandomPoints.R

# Run this after generating the NDVI rasters in GEE
# Find difference between the 2 NDVI rasters, mask out selected NDVI classes, classify the differences, then find random points 
# based on the classes

library(rgdal)
library(raster)
require(rgeos) # gbuffer
library(classInt) # classIntervals
setwd("N:\\project\\monitoring\\SoCalShrubs\\NDVI")

# NDVI rasters from GEE
springNDVI <- raster("absmaxNDVI_0401_0515_2018.tif")
fallNDVI <- raster("absminNDVI_0901_0930_2018.tif")

nlcd <- raster("N:\\rasterlib\\NLCD\\nlcd_2006_landcover_2011_edition_2014_10_10\\nlcd_2006_landcover_2011_edition_2014_10_10.img")
nlcdc <- crop(nlcd, springNDVI)
setExtent(nlcdc,extent(springNDVI),keepres=TRUE,snap=TRUE );origin(nlcdc) <- origin(springNDVI)

m = c(10.5,12,NA,
      21.5,24,NA,
      30.5,31,NA,
      40.5,43,NA,
      80.5,82,NA)

rmatrix <- matrix(m, ncol=3, byrow=TRUE)
nlcdr <- reclassify(nlcdc, rmatrix)
writeRaster(nlcdr, "nlcdr.tif",overwrite=TRUE)

values(nlcdr)[values(nlcdr) == 0] <- NA

values(springNDVI)[values(springNDVI) == -9999] = NA
values(springNDVI)[values(springNDVI) == 0] = NA
values(fallNDVI)[values(fallNDVI) == -9999] = NA
values(fallNDVI)[values(fallNDVI) == 0] = NA

diffNDVI <- springNDVI - fallNDVI

diffNDVIc <- crop(diffNDVI, nlcdr)
nlcdrc <- crop(nlcdr, diffNDVIc)
diffNDVI_m <- mask(diffNDVIc, nlcdrc, maskvalue=NA) 

writeRaster(diffNDVI_m, "maskeddiffndvi.tif",overwrite=TRUE)

plot(diffNDVI_m, main="NDVI Difference", axes=FALSE)
hist(diffNDVI_m[!diffNDVI_m == 0],breaks=40,col="springgreen4",main="NDVI Difference",ylab="Num Pixels",xlab="NDVI Diff")

# Generate classes based on values of diffNDVI
GDALinfo("diffndvi.tif")

m = c(-500000,75,2,
      75,130,3,
      130,170,4,
      170,205,5,
      205,243,6,
      243,290,7,
      290,360,8,
      360,470,9,
      470,650,10,
      650,910,11,
      910,1200,12,
      1200,1500,13,
      1500,1800,14,
      1800,2130,15,
      2130,2500,16,
      2500,3000,17,
      3000,3780,18,
      3780,65000,19)

reclassmatrix <- matrix(m, ncol=3, byrow=TRUE)
reclassmatrix
reclassed <- reclassify(diffNDVI_m,reclassmatrix)
#writeRaster(reclassed, "maskedreclassed.tif",overwrite=TRUE)
f <- freq(reclassed)
f

#  ------------ Get random stratified sample & write squares for each point  ------------
# take 500 from each of the 20 classes
samp <- sampleStratified(reclassed,size=555,na.rm=TRUE, xy=TRUE, sp=TRUE)
dsn <- file.path(getwd())
outfile <- '\\samp_pts18class'
writeOGR(samp, dsn,outfile, driver="ESRI Shapefile", overwrite_layer = TRUE, morphToESRI = TRUE)

buff30 = gBuffer(samp,byid = TRUE, width=15.0,capStyle="SQUARE")
writeOGR(buff30, dsn,"\\samp10k_poly18cl", driver = "ESRI Shapefile", overwrite_layer = TRUE, morphToESRI = TRUE)

#   ------------ Select more points in 3 specific classes of NLCD & write those squares ------------
nlcdr <- raster("nlcdr.tif")
values(nlcdr)[values(nlcdr) == 0] <- NA

pts <- which(getValues(nlcdr) %in% c(21,71))
extent(nlcdr)
tmp <- raster(ncols=ncol(nlcdr),nrows=nrow(nlcdr),res=res(nlcdr),crs=crs(nlcdr),ext=extent(nlcdr),vals=NULL)
tmp[] <- NA
tmp[pts] <- nlcdr[pts]
samp2171 <- sampleStratified(tmp,size=1000,na.rm=TRUE, xy=TRUE, sp=TRUE)

tmp[] <- NA
pts <- which(getValues(nlcdr) == 95)
tmp[pts] <- nlcdr[pts]
samp95 <- sampleStratified(tmp,size=17190,na.rm=TRUE, xy=TRUE, sp=TRUE) 

allsamp <- rbind(samp2171,samp95)

dsn <- file.path(getwd())
outfile <- '\\samp_ptsnlcd'

writeOGR(allsamp, dsn,outfile, driver="ESRI Shapefile", overwrite_layer = TRUE, morphToESRI = TRUE)
buff30 = gBuffer(allsamp,byid = TRUE, width=15.0,capStyle="SQUARE")
writeOGR(buff30, dsn,"\\sampnlcd_poly", driver = "ESRI Shapefile", overwrite_layer = TRUE, morphToESRI = TRUE)




