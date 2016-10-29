
#-----------------------
#KFW 1 Cross-Sectional Model
#Treatment: Early Demarcated through PPTAL by April 2001 (vs. demarcated after April 2001)
#Outcome: Median NDVI change in level from 1995-2001, Median NDVI change in level from 2001-2010
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

#NDVI Max Trends for 1995-2001
dta_Shp$post_trend_NDVI_95_01 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2001,"id")
dta_Shp@data["NDVILevelChange_95_01"] <- dta_Shp$MaxL_2001 - dta_Shp$MaxL_1995
#dta_Shp@data["NDVIslopeChange_95_01"] <- dta_Shp@data["post_trend_NDVI_95_01"] - dta_Shp@data["pre_trend_NDVI_max"]

#NDVI Max Trends for 2001-2010
dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",2001,2010,"id")
dta_Shp@data["NDVILevelChange_01_10"] <- dta_Shp$MaxL_2010 - dta_Shp$MaxL_2001
#dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp@data["post_trend_NDVI_01_10"] - dta_Shp@data["pre_trend_NDVI_max"]

#NDVI Med Trends for 1995-2001
dta_Shp$post_trend_NDVI_95_01_Med <- timeRangeTrend(dta_Shp,"MedL_[0-9][0-9][0-9][0-9]",1995,2001,"id")
dta_Shp@data["NDVILevelChange_95_01_Med"] <- dta_Shp$MedL_2001 - dta_Shp$MedL_1995

#NDVI Med Trends for 2001-2010
dta_Shp$post_trend_NDVI_01_10_Med <- timeRangeTrend(dta_Shp,"MedL_[0-9][0-9][0-9][0-9]",2001,2010,"id")
dta_Shp@data["NDVILevelChange_01_10_Med"] <- dta_Shp$MedL_2010 - dta_Shp$MedL_2001

#Calculate Temp and Precip Pre and Post Trends
dta_Shp$pre_trend_temp_mean <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_temp_max <- timeRangeTrend(dta_Shp,"MaxT_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_temp_min <- timeRangeTrend(dta_Shp,"MinT_[0-9][0-9][0-9][0-9]",1982,1995,"id")

dta_Shp$post_trend_temp_95_01 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",1995,2001,"id")
dta_Shp$post_trend_temp_01_10 <- timeRangeTrend(dta_Shp,"MeanT_[0-9][0-9][0-9][0-9]",2001,2010,"id")

dta_Shp$pre_trend_precip_mean <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_precip_max <- timeRangeTrend(dta_Shp,"MaxP_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_precip_min <- timeRangeTrend(dta_Shp,"MinP_[0-9][0-9][0-9][0-9]",1982,1995,"id")

dta_Shp$post_trend_precip_95_01 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",1995,2001,"id")
dta_Shp$post_trend_precip_01_10 <- timeRangeTrend(dta_Shp,"MeanP_[0-9][0-9][0-9][0-9]",2001,2010,"id")

#-------------------------------------------------
#-------------------------------------------------
#Define the Treatment Variable and Population
#-------------------------------------------------
#-------------------------------------------------
#Make a binary to test treatment..
dta_Shp@data["TrtBin"] <- 0
dta_Shp@data$TrtBin[dta_Shp@data$demend_y <= 2001] <- 1
dta_Shp@data$TrtBin[(dta_Shp@data$demend_m > 4) & (dta_Shp@data$demend_y==2001)] <- 0

#Remove units that did not ever receive any treatment (within-sample test)
dta_Shp@data$NA_check <- 0
dta_Shp@data$NA_check[is.na(dta_Shp@data$demend_y)] <- 1
int_Shp <- dta_Shp[dta_Shp@data$NA_check != 1,]
dta_Shp <- int_Shp
table(dta_Shp@data$TrtBin)


#--------------------------
#Matching, with replacement
#--------------------------

#identify vars needed for psm model and analytic models (only these will be included in new matched dataset)
aVars <- c("reu_id","UF","TrtBin", "terrai_are","Pop_1990","Pop_2000", "MeanT_1995","MeanT_2001", "pre_trend_temp_mean",
           "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995", "MeanP_2001","pre_trend_precip_min", 
           "pre_trend_NDVI_mean", "pre_trend_NDVI_max","pre_trend_NDVI_med",
           "MaxL_1995","MedL_1995","NDVILevelChange_95_01","NDVILevelChange_01_10",
           "NDVILevelChange_95_01_Med","NDVILevelChange_01_10_Med",
           "Slope","Elevation","Riv_Dist","Road_dist",
           "pre_trend_precip_mean", "pre_trend_precip_max",
           "post_trend_temp_95_01","post_trend_temp_01_10",
           "post_trend_precip_95_01","post_trend_precip_01_10")

#propensity score model
#replace=TRUE to match with replacement
#exact="UF" restricts matches to same Brazilian state and discard="both" must accompany it
#use resulting weights in models to account for matching with replacement

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

##create standardized dataset to produce standardized coefficients in models that are easy to output

#identify vars for inclusion in standardized dataset
#include all numeric variables from psm equation and that will be included in models
#exclude any id fields and weights created from matchit
# stvars <- c("TrtBin", "terrai_are","Pop_1990","Pop_2000" ,"MeanT_1995","MeanT_2001", "pre_trend_temp_mean",
#             "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995","MeanP_2001", "pre_trend_precip_min", 
#             "pre_trend_NDVI_mean", "pre_trend_NDVI_max","pre_trend_NDVI_med",
#             "NDVILevelChange_95_01_Med","NDVILevelChange_01_10_Med",
#             "MaxL_1995","MedL_1995",
#             "Slope","Elevation","Riv_Dist","Road_dist",
#             "pre_trend_precip_mean", "pre_trend_precip_max",
#             "NDVILevelChange_95_01","NDVILevelChange_01_10","post_trend_temp_95_01","post_trend_temp_01_10",
#             "post_trend_precip_95_01","post_trend_precip_01_10")
# 
# model_data_st<- model_data
# model_data_st[stvars]<-lapply(model_data_st[stvars],scale)

#--------------
#Analytic Models
#--------------

##Early Models, Outcome: 1995-2001 Median

#Create dataset with some common names for stargazer
model_data_early <- model_data
colnames(model_data_early)[(colnames(model_data_early)=="Pop_1990")] <- "Pop_B"
colnames(model_data_early)[(colnames(model_data_early)=="MeanT_1995")] <- "MeanT_B"
colnames(model_data_early)[(colnames(model_data_early)=="MeanP_1995")] <- "MeanP_B"
colnames(model_data_early)[(colnames(model_data_early)=="post_trend_temp_95_01")] <- "post_trend_temp"
colnames(model_data_early)[(colnames(model_data_early)=="post_trend_precip_95_01")] <- "post_trend_precip"

#ModelEarly2, treatment effect + weights, 1995-2001 Median
ModelEarly2_Med <- lm(NDVILevelChange_95_01_Med ~ TrtBin, data=model_data, weights=(weights))

#ModelEarly3, treatment effect + weights + covars, 1995-2001 Median
ModelEarly3_Med<-lm(NDVILevelChange_95_01_Med~TrtBin +pre_trend_NDVI_med + MedL_1995 + terrai_are+Pop_B+
                      MeanT_B + post_trend_temp+
                      MeanP_B + post_trend_precip+
                      Slope+Elevation+Riv_Dist+Road_dist,
                    data=model_data_early,
                    weights=(weights))

##Late Models

#Create dataset with some common names for stargazer
model_data_late<-model_data
colnames(model_data_late)[(colnames(model_data_late)=="Pop_2000")] <- "Pop_B"
colnames(model_data_late)[(colnames(model_data_late)=="MeanT_2001")] <- "MeanT_B"
colnames(model_data_late)[(colnames(model_data_late)=="MeanP_2001")] <- "MeanP_B"
colnames(model_data_late)[(colnames(model_data_late)=="post_trend_temp_01_10")] <- "post_trend_temp"
colnames(model_data_late)[(colnames(model_data_late)=="post_trend_precip_01_10")] <- "post_trend_precip"

#ModelLate, treatment effect + weights + covars, 2001-2010 Median
ModelLate_Med<-lm(NDVILevelChange_01_10_Med~TrtBin+ pre_trend_NDVI_med + MedL_1995+terrai_are+Pop_B+
                    MeanT_B+post_trend_temp+
                    MeanP_B + post_trend_precip+
                    Slope + Elevation + Riv_Dist + Road_dist,
                  data=model_data_late,
                  weights=(weights))



#-------------
#Stargazer
#-------------

stargazer(ModelEarly2_Med, ModelEarly3_Med,ModelLate_Med,
          keep=c("TrtBin", "pre_trend_NDVI_med","MedL_1995", "terrai_are","Pop_B","MeanT_B","post_trend_temp","MeanP_B",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment (Early Demarcation)", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends","Baseline Precipitation","Precipitation Trends", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Median NDVI 1995-2010"," Median NDVI 2001-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)




