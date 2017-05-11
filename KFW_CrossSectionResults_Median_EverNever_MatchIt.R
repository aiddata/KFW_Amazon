
#-----------------------
#KFW 1 Cross-Sectional Model
#Treatment: Ever Demarcated through PPTAL
#Outcome: Median NDVI change in level from 1995-2010
#Using MatchIt package instead of SCI
#-----------------------

library(devtools)
devtools::install_github("itpir/SCI@master")
library(SCI)
library(stargazer)
loadLibs()
library(MatchIt)
library(rgeos)
library(maptools)
library(rgdal)
library(sp)


shpfile = "processed_data/kfw_analysis_inputs.shp"
dta_Shp = readShapePoly(shpfile)

#-------------------------------------------------
#-------------------------------------------------
#Pre-processing to create cross-sectional variable summaries
#-------------------------------------------------
#-------------------------------------------------
#Calculate NDVI Trends
dta_Shp$pre_trend_NDVI_mean <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_NDVI_max <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_NDVI_med <- timeRangeTrend(dta_Shp,"MedL_[0-9][0-9][0-9][0-9]",1982,1995,"id")

dta_Shp$pre_trend_NDVI_mean <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_NDVI_max <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$NDVIslope_95_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp@data["NDVILevelChange_95_10"] <- dta_Shp$MaxL_2010 - dta_Shp$MaxL_1995
dta_Shp@data["NDVILevelChange_95_10_Med"] <- dta_Shp$MedL_2010 - dta_Shp$MedL_1995

#Calculate Temp and Precip Pre and Post Trends
dta_Shp$pre_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1982,1995,"id")

dta_Shp$post_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp$post_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp$post_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1995,2010,"id")

dta_Shp$pre_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1982,1995,"id")

dta_Shp$post_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp$post_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp$post_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1995,2010,"id")

#-------------------------------------------------
#-------------------------------------------------
#Define the Treatment Variable and Population
#-------------------------------------------------
#-------------------------------------------------
#Eliminate non-PPTAL indigenous lands
dta_Shp@data$proj_check <- 0
dta_Shp@data$proj_check[is.na(dta_Shp@data$reu_id)] <- 1
proj_Shp <- dta_Shp[dta_Shp@data$proj_check !=1,]
dta_Shp <- proj_Shp

projtable <- table(proj_Shp@data$proj_check)
View(projtable)

#Make a binary for ever treated vs. never treated
dta_Shp@data["TrtBin"] <- 0
dta_Shp@data$NA_check <- 0
dta_Shp@data$NA_check[is.na(dta_Shp@data$demend_y)] <- 1
dta_Shp@data$TrtBin[dta_Shp@data$NA_check != 1] <- 1

demtable <- table(dta_Shp@data$TrtBin)
View(demtable)

#--------------------------
#Matching, with replacement
#--------------------------

aVars <- c("reu_id","UF","TrtBin", "terrai_are","Pop_1990", "MeanT_1995", "pre_trend_temp_mean",
           "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995", "pre_trend_precip_min", 
           "pre_trend_NDVI_mean", "pre_trend_NDVI_max","pre_trend_NDVI_med",
           "Slope","Elevation","MaxL_1995","MedL_1995","Riv_Dist","Road_dist",
           "pre_trend_precip_mean", "pre_trend_precip_max",
           "NDVILevelChange_95_10","NDVILevelChange_95_10_Med",
           "post_trend_temp_mean","post_trend_temp_min","post_trend_temp_max",
           "post_trend_precip_mean","post_trend_precip_min","post_trend_precip_max")


psmModel <- matchit(TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                      pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
                      pre_trend_NDVI_mean + pre_trend_NDVI_max + pre_trend_NDVI_med+
                      Slope + Elevation + MedL_1995 + Riv_Dist + Road_dist +
                      pre_trend_precip_mean + pre_trend_precip_max,
                    data=dta_Shp@data[aVars],
                    method="nearest",replace=TRUE, exact="UF",discard="both")

print(summary(psmModel))
#create new dataset with matches
model_data<-match.data(psmModel)

#check states that were dropped out
summary(model_data$UF)

#--------------
#Analytic Models
#--------------

##Matched Data Models

model2 <- lm(NDVILevelChange_95_10_Med ~ TrtBin, data=model_data, weights=(weights))

model3<-lm(NDVILevelChange_95_10_Med ~ TrtBin + pre_trend_NDVI_med + MedL_1995 + terrai_are + Pop_1990 + MeanT_1995 + 
             post_trend_temp_mean + 
             post_trend_precip_mean + 
             MeanP_1995 + Slope + Elevation  + Riv_Dist + Road_dist, data=model_data, weights=(weights))


#-------------
#Stargazer
#-------------

stargazer(model2, model3,
          keep=c("TrtBin", "pre_trend_NDVI_med","MedL_1995", "terrai_are","Pop_1990","MeanT_1995","post_trend_temp","MeanP_1995",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Precipitation Trends","Baseline Precipitation", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Median NDVI 1995-2010"),
          title="Regression Results", type="html", 
          star.cutoffs = c(0.05, 0.01, 0.001),
          omit.stat=c("f","ser"), align=TRUE)




