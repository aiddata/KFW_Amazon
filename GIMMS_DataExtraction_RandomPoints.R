#R script to drop 100,000 random points in KFW communities,
#weighted by area.

library(maptools)
library(sp)

get_area <- function (x)
{
  #prj <- CRS("+proj=longlat +ellps=clrk66")
  xx <- x
  library(rgdal)
  #Lambert Equal Area
  xxx <- spTransform(xx, CRS("+init=epsg:3571"))

  sapply(slot(xxx, "polygons"), function(x) length(slot(x, "Polygons")))
  
  sapply(slot(xxx, "polygons"), function(x) {
    xi <- slot(x, "Polygons")
    any(sapply(xi, slot, "hole"))
  })
  # no declared holes
  sapply(slot(xxx, "polygons"), slot, "area")
  # areas in square metres
}


KFW <- readShapePoly("/home/dan/Desktop/GitRepo/KFW_Amazon/processed_data/kfw_analysis_inputs.shp", proj4string=CRS("+proj=longlat +datum=WGS84"))

KFW_noNA <- KFW[which(!is.na(KFW@data$demend_y)),]

KFW_post2001 <- KFW_noNA[KFW_noNA@data$demend_y >= 2001,]

KFW_post2001$area <- get_area(KFW_post2001)



tot_area <- sum(KFW_post2001$area)
samples = 100000

all_pts <- matrix(as.numeric(NA), nlevels(KFW_post2001$reu_id), 1)
all_pts <- as.list(all_pts)

reu_ids <- matrix(as.numeric(NA), nlevels(KFW_post2001$reu_id), 1)
reu_ids <- as.list(reu_ids)

for(i in 1:length(KFW_post2001))
{
  
  num_points = round(KFW_post2001@data$area[[i]]/tot_area, 5) * samples
  
  cur_id = KFW_post2001@data$reu_id[[i]]
  KFW_SinglePoly <- KFW_post2001[which(KFW_post2001@data$reu_id == cur_id),]

  rnd_pts <- spsample(KFW_SinglePoly, n=num_points, "stratified")
  all_pts[i] <- rnd_pts
  reu_ids[i] <- cur_id
}

final_list <- matrix(as.numeric(NA), nlevels(KFW_post2001$reu_id), 1)
final_list <- as.list(final_list)

for(i in 1:length(all_pts))
{
  reu_id <- as.data.frame(matrix(reu_ids[[i]],length(rownames(all_pts[[i]]@coords)),1))
  pts_dataframe = SpatialPointsDataFrame(all_pts[[i]]@coords, reu_id, match.ID=FALSE)
  final_list[i] = pts_dataframe
}

#merge(all_pts,reu_ids)
KFW_allPts_post2001 <- do.call("rbind", final_list)

## Plot 
plot(KFW_post2001)
plot(KFW_allPts_post2001 , add=TRUE)

writePointsShape(KFW_allPts_post2001, "/home/dan/Desktop/GitRepo/KFW_Amazon/interimData/rand_pts_r0.shp")

