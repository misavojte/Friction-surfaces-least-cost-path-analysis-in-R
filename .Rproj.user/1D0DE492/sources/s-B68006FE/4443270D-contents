#Lib managment
if(!require("pacman")) install.packages("pacman")
pacman::p_load("raster", "sf", "sp","readxl","rgdal")
  
#read table for speed through each of CORINE landcover types
speedThroughTerrainByTypes <- readxl::read_xls("input/clckoodid.xls")

#get area of interest shapefile
areaOfInterest <- sf::st_read("input/Service_Area.shp")

# prepare it for masking purposes
areaOfInterest <- sf::st_union(areaOfInterest$geometry)

#prepare raster template for rasterization of vectors ====
rasterTemplate <- raster()
raster::extent(rasterTemplate) <- sf::st_bbox(areaOfInterest)
raster::crs(rasterTemplate) <- raster::crs(areaOfInterest)
raster::res(rasterTemplate) <- 50 #set raster res to 50 m

#fun for loading, uniting crs and masking new layers by area of interest
loadNewShapefile <- function(pathToShp) {
  newLayer <- sf::st_transform(st_read(pathToShp), raster::crs(areaOfInterest))
  return(sf::st_intersection(newLayer, areaOfInterest))
}

corine <- loadNewShapefile("input/clc00ee_L-EST97.shp")

#keep only relevant codes for getting speeds
corine <- corine[,c("CODE_00")]

#assign speeds to CORINE
for (i in 1:nrow(speedThroughTerrainByTypes)) {
  currentCode <- as.numeric(speedThroughTerrainByTypes[i,"code"])
  currentSpeedForCode <- as.numeric(speedThroughTerrainByTypes[i,"speed_km_h"])
  corine[corine$CODE_00==currentCode,"speed"] <- currentSpeedForCode
}

#plot(corine[,c("speed")])

#rasterize corine
#keep highest possible speed for each cell (fun="max")
#(this can take a while)
corineRaster <- raster::rasterize(corine,rasterTemplate,field="speed",fun="max")

#get rivers (barriers)
rivers <- loadNewShapefile("input/Rivers.shp")
rivers <- rivers[,c("geometry")]
riversRaster <- raster::rasterize(rivers,rasterTemplate)

#prepare rivers raster for multiplying
riversRaster$layer[riversRaster$layer>0] <- 0 #barrier - for multiplying
riversRaster$layer[is.na(riversRaster$layer)] <- 1 #not speed! just for keeping values from corineRaster

#get roads and set their speed
mainRoads <- loadNewShapefile("input/MainRoads.shp")
mainRoads <- mainRoads[,c("geometry")]
mainRoadsRaster <- raster::rasterize(mainRoads,rasterTemplate)
mainRoadsRaster$layer[mainRoadsRaster$layer>0] <- 90 #speed in km/h

localRoads <- loadNewShapefile("input/LocalRoads.shp")
localRoads <- localRoads[,c("geometry")]
localRoadsRaster <- raster::rasterize(localRoads,rasterTemplate)
localRoadsRaster$layer[localRoadsRaster$layer>0] <- 80 #speed in km/h

otherRoads <- loadNewShapefile("input/OtherRoads.shp")
otherRoads <- otherRoads[,c("geometry")]
otherRoadsRaster <- raster::rasterize(otherRoads,rasterTemplate)
otherRoadsRaster$layer[otherRoadsRaster$layer>0] <- 60 #speed in km/h

#start generating speed raster
speedRaster <- corineRaster * riversRaster
speedRaster <- raster::stack(mainRoadsRaster,localRoadsRaster,otherRoadsRaster,speedRaster)
speedRaster <- raster::stackApply(speedRaster,indices=c(1,1,1,1),fun="max")

#from km/h to m/min (prepare for accumulation cost analysis)
sR2 <- speedRaster * 16.66667
sR2[sR2 == 0] <- NA
sR3 <- 50/sR2

# create a transition object from the raster for accumulation cost analysis
tr <- gdistance::transition(sR3,function(x) 1/mean(x),8)

#load emergency stations
stations <- loadNewShapefile("input/EmrgStat.shp")

# create accumulation cost raster from stations
ac <- gdistance::accCost(tr,sf::as_Spatial(sf::st_as_sf(stations)))

raster::writeRaster(ac,"output/accCost.tiff")
