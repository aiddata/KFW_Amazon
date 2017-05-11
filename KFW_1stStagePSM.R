#-------------------------------------------------
#-------------------------------------------------
#Cross Sectional Models - KFW
#Testing in Cross Section the impact of being treated AFTER April 2001
#On the Max Level of NDVI, measured as a change in the level of NDVI between start and end year (1995-2001, 2001-2010)
#-------------------------------------------------
#-------------------------------------------------
library(devtools)
devtools::install_github("itpir/SAT@master")
library(SAT)
library(stargazer)
loadLibs()
#-------------------------------------------------
#-------------------------------------------------
#Load in Processed Data - produced from script KFW_dataMerge.r
#-------------------------------------------------
#-------------------------------------------------

shpfile = "processed_data/kfw_analysis_inputs.shp"
dta_Shp = readShapePoly(shpfile)

#-------------------------------------------------
#-------------------------------------------------
#Pre-processing to create cross-sectional variable summaries
#-------------------------------------------------
#-------------------------------------------------
#Calculate NDVI Trends
dta_Shp$pre_trend_NDVI_mean <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$pre_trend_NDVI_max <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$NDVIslope_95_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp@data["NDVILevelChange_95_10"] <- dta_Shp$MaxL_2010 - dta_Shp$MaxL_1995

#NDVI Trends for 1995-2001
dta_Shp$post_trend_NDVI_95_01 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
dta_Shp@data["NDVILevelChange_95_01"] <- dta_Shp$MaxL_2001 - dta_Shp$MaxL_1995
#dta_Shp@data["NDVIslopeChange_95_01"] <- dta_Shp@data["post_trend_NDVI_95_01"] - dta_Shp@data["pre_trend_NDVI_max"]

#NDVI Trends for 2001-2010
dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")
dta_Shp@data["NDVILevelChange_01_10"] <- dta_Shp$MaxL_2010 - dta_Shp$MaxL_2001
#dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp@data["post_trend_NDVI_01_10"] - dta_Shp@data["pre_trend_NDVI_max"]

#Calculate Temp and Precip Pre and Post Trends
dta_Shp$pre_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$pre_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$pre_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")

dta_Shp$post_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_temp_95_01 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
dta_Shp$post_trend_temp_01_10 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")

dta_Shp$pre_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$pre_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")
dta_Shp$pre_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1982,1995,"SP_ID")

dta_Shp$post_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp$post_trend_precip_95_01 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
dta_Shp$post_trend_precip_01_10 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")

#-------------------------------------------------
#-------------------------------------------------
#EVER VS. NEVER: Define the Treatment Variable and Population
#-------------------------------------------------
#-------------------------------------------------
#Make a binary to test treatment..
dta_Shp_ever <- dta_Shp

#Eliminate non-PPTAL indigenous lands
dta_Shp_ever@data$proj_check <- 0
dta_Shp_ever@data$proj_check[is.na(dta_Shp_ever@data$reu_id)] <- 1
proj_Shp_ever <- dta_Shp_ever[dta_Shp_ever@data$proj_check !=1,]
dta_Shp_ever <- proj_Shp_ever

projtable_ever <- table(proj_Shp_ever@data$proj_check)
View(projtable_ever)


#Make a binary for ever treated vs. never treated
dta_Shp_ever@data["TrtBin"] <- 0
dta_Shp_ever@data$NA_check <- 0
dta_Shp_ever@data$NA_check[is.na(dta_Shp_ever@data$demend_y)] <- 1
dta_Shp_ever@data$TrtBin[dta_Shp_ever@data$NA_check != 1] <- 1

demtable_ever <- table(dta_Shp_ever@data$TrtBin)
View(demtable_ever)

#-------------------------------------------------
#-------------------------------------------------
#EARLY VS. LATE: Define the Treatment Variable and Population
#-------------------------------------------------
#-------------------------------------------------
#Make a binary to test treatment..
dta_Shp_early <- dta_Shp
dta_Shp_early@data["TrtBin"] <- 0
dta_Shp_early@data$TrtBin[dta_Shp_early@data$demend_y <= 2001] <- 1
dta_Shp_early@data$TrtBin[(dta_Shp_early@data$demend_m > 4) & (dta_Shp_early@data$demend_y==2001)] <- 0

#Remove units that did not ever receive any treatment (within-sample test)
dta_Shp_early@data$NA_check <- 0
dta_Shp_early@data$NA_check[is.na(dta_Shp_early@data$demend_y)] <- 1
int_Shp_early <- dta_Shp_early[dta_Shp_early@data$NA_check != 1,]
dta_Shp_early <- int_Shp_early
table(dta_Shp_early@data$TrtBin)

#-------------------------------------------------
#-------------------------------------------------
#Define and run the first-stage of the PSM, calculating propensity scores
#-------------------------------------------------
#-------------------------------------------------
psmModel <-  "TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
pre_trend_NDVI_mean + pre_trend_NDVI_max + Slope + Elevation + MaxL_1995 + Riv_Dist + Road_dist +
pre_trend_precip_mean + pre_trend_precip_max"
#MeanL_1995

psmRes_ever <- SAT::SpatialCausalPSM(dta_Shp_ever,mtd="logit",psmModel,drop="support",visual=TRUE)

psmRes_early <- SAT::SpatialCausalPSM(dta_Shp_early,mtd="logit",psmModel,drop="support",visual=TRUE)


##-----------------------------------------------
#Stargazer
##-----------------------------------------------


stargazer(psmRes_ever$model, psmRes_early$model,
          type="html", align=TRUE,
          covariate.labels=c("Pre Trend Min Precipitation","Pre Trend Mean Precipitation","Pre Trend Max Precipitation",
          "Pre Trend Mean Temperature","Pre Trend Min Temperature","Pre Trend Max Temperature",
          "Pre Trend NDVI Mean","Pre Trend NDVI Max","Max NDVI Baseline",
          "Population Baseline","Mean Temperature Baseline","Mean Precipitation Baseline", "Area (hectares)",
          "Slope","Elevation","Distance to River","Distance to Road"),
          omit.stat=c("f","ser"),
          order=c("pre_trend_precip","pre_trend_temp","pre_trend_NDVI","MaxL_","Pop","Mean"),
          title="PSM First Stage Results",
          star.cutoffs = c(0.05, 0.01, 0.001),
          dep.var.labels=c("Ever Demarcated   Early vs. Late"))


