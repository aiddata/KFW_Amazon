#-----------------------
#KFW 1 Cross-Sectional Model
#Treatment: Demarcated through PPTAL (vs. never demarcated)
#Outcome: Max NDVI change in level from 1995-2010
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
dta_Shp$pre_trend_NDVI_mean <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$pre_trend_NDVI_max <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1982,1995,"id")
dta_Shp$NDVIslope_95_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",1995,2010,"id")
dta_Shp@data["NDVILevelChange_95_10"] <- dta_Shp$MaxL_2010 - dta_Shp$MaxL_1995

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

#plot demarcation year and NDVI pre-trend
# reg<-lm(dta_Shp@data$pre_trend_NDVI_max~dta_Shp@data$demend_y)
# summary(reg)
# plot(dta_Shp@data$demend_y,dta_Shp@data$pre_trend_NDVI_max,
#      xlab="Year of Demarcation",ylab="NDVI Max Pre-Trend")
# 
# abline(lm(dta_Shp@data$pre_trend_NDVI_max~dta_Shp@data$demend_y))
# abline(h=0)

ggplot(dta_Shp@data, aes(x = demend_y, y = pre_trend_NDVI_max)) + 
  geom_point() +
  stat_smooth(method = "lm",se=TRUE,col="black")+
  geom_hline(yintercept=.00249109,linetype="dashed")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="Year of Demarcation",y="NDVI Pre-Trend")+
  theme_bw()

#--------------------------
#Matching, with replacement
#Uses MatchIt
#--------------------------

aVars <- c("reu_id","UF","TrtBin", "terrai_are","Pop_1990", "MeanT_1995", "pre_trend_temp_mean",
           "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995", "pre_trend_precip_min", 
           "pre_trend_NDVI_mean", "pre_trend_NDVI_max","Slope","Elevation","MaxL_1995","MeanL_1995","Riv_Dist","Road_dist",
           "pre_trend_precip_mean", "pre_trend_precip_max",
           "NDVILevelChange_95_10","post_trend_temp_mean","post_trend_temp_min","post_trend_temp_max",
           "post_trend_precip_mean","post_trend_precip_min","post_trend_precip_max")


psmModel <- matchit(TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
                      pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
                      pre_trend_NDVI_mean + pre_trend_NDVI_max + Slope + Elevation + MaxL_1995 + Riv_Dist + Road_dist +
                      pre_trend_precip_mean + pre_trend_precip_max,
                    data=dta_Shp@data[aVars],
                    method="nearest",replace=TRUE, exact="UF",discard="both")

print(summary(psmModel))

model_data<-match.data(psmModel)

#create variable to use for weighting in models (1/pscore * weight)
model_data$psmweight<-(1/model_data$distance)*(model_data$weights)

#-------------------------
#Matching without Replacement
#Using old SCI package, not MatchIt (to preserve pair ids, and results from initial journal submission)
#-------------------------

psmModel <-  "TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
pre_trend_NDVI_mean + pre_trend_NDVI_max + Slope + Elevation + MaxL_1995 + Riv_Dist + Road_dist +
pre_trend_precip_mean + pre_trend_precip_max"
#MeanL_1995

psmRes <- SAT::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel,drop="support",visual=FALSE)

#-------------------------------------------------
#-------------------------------------------------
#Based on the Propensity Score Matches, pair comprable treatment and control units.
#-------------------------------------------------
#-------------------------------------------------
drop_set<- c(drop_unmatched=TRUE,drop_method="None",drop_thresh=0.5)
psm_Pairs <- SAT(dta = psmRes$data, mtd = "fastNN",constraints=c(groups="UF"),psm_eq = psmModel, ids = "id", drop_opts = drop_set, visual="TRUE", TrtBinColName="TrtBin")
#c(groups=c("UF"),distance=NULL)
trttable <- table (psm_Pairs@data$TrtBin)
View(trttable)


##create standardized dataset to produce standardized coefficients in models that are easy to output
#or to create normalized difference of means statistics for summary statistics

stvars <- c("TrtBin", "terrai_are","Pop_1990", "MeanT_1995", "pre_trend_temp_mean",
           "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995", "pre_trend_precip_min",
           "pre_trend_NDVI_mean", "pre_trend_NDVI_max","Slope","Elevation","MaxL_1995","MeanL_1995","Riv_Dist","Road_dist",
           "pre_trend_precip_mean", "pre_trend_precip_max",
           "NDVILevelChange_95_10","post_trend_temp_mean","post_trend_temp_min","post_trend_temp_max",
           "post_trend_precip_mean","post_trend_precip_min","post_trend_precip_max")
#
model_data_st<- model_data
model_data_st[stvars]<-lapply(model_data_st[stvars],scale)

psm_Pairs_st<-psm_Pairs@data
psm_Pairs_st[stvars]<-lapply(psm_Pairs_st[stvars],scale)
#
dta_Shp_st<-dta_Shp@data
dta_Shp_st[stvars]<-lapply(dta_Shp_st[stvars],scale)

#-------------
#MODELS
#-------------

##Unmatched Data Models
#no propensity score matching (for comparison)
model2u <- lm(NDVILevelChange_95_10 ~ TrtBin, data=dta_Shp@data)

model3u<-lm(NDVILevelChange_95_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_1990 + MeanT_1995 + 
              post_trend_temp_mean + 
              post_trend_precip_mean + 
              MeanP_1995 + Slope + Elevation  + Riv_Dist + Road_dist, data=dta_Shp@data)

## Matching without Replacement, Data Models (also uses old SCI package, not MatchIt)
#analyticModelEver2, pair FEs, no covars, 1995-2010
analyticModelEver2 <- "NDVILevelChange_95_10 ~ TrtBin + factor(PSM_match_ID)"

OutputEver2=Stage2PSM(analyticModelEver2,psm_Pairs,type="lm",table_out=TRUE)

#analyticModelEver3, pair FEs, covars, 1995-2010
analyticModelEver3 <- "NDVILevelChange_95_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_1990 + MeanT_1995  +
MeanP_1995 + post_trend_temp_mean + post_trend_precip_mean + Slope + Elevation  + Riv_Dist + Road_dist + 
factor(PSM_match_ID)"

OutputEver3=Stage2PSM(analyticModelEver3,psm_Pairs,type="lm",table_out=TRUE)

##Matching with Replacement, Data Models
#when use model_data_st, coefficients are standardized

model2 <- lm(NDVILevelChange_95_10 ~ TrtBin, data=model_data, weights=(psmweight))

model3<-lm(NDVILevelChange_95_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_1990 + MeanT_1995 + 
             post_trend_temp_mean + 
             post_trend_precip_mean + 
             MeanP_1995 + Slope + Elevation  + Riv_Dist + Road_dist, data=model_data, weights=(psmweight))


#-------------
#Stargazer
#-------------

stargazer(model2u, model3u,
          keep=c("TrtBin", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_1990","MeanT_1995","post_trend_temp","MeanP_1995",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Precipitation Trends","Baseline Precipitation", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)


stargazer(model2, model3,
          keep=c("TrtBin", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_1990","MeanT_1995","post_trend_temp","MeanP_1995",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Precipitation Trends","Baseline Precipitation", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)

stargazer(model2u, model3u,model2, model3,
          keep=c("TrtBin", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_1990","MeanT_1995","post_trend_temp","MeanP_1995",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Precipitation Trends","Baseline Precipitation", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)

## FINAL TABLE FOR JEEM REVISION

stargazer(model2u, model3u,OutputEver2$unstandardized,OutputEver3$unstandardized,model2, model3,
          omit=c("factor"),
          keep=c("TrtBin", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_1990","MeanT_1995","post_trend_temp","MeanP_1995",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Precipitation Trends","Baseline Precipitation", 
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Change in Max NDVI: 2010 Level - 1995 Level"),
          column.labels=c("Unmatched","Matched w/o Replacement","Matched w/ Replacement"),
          column.separate=c(2,2,2),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)


#---------------
#Summary Stats, by groups
#this provides values and then enter in Excel to format it as needed
#---------------


##Creating normalized difference in means, by groups
#to create statistic that allows comparison across unmatched, matched w/o replacement, matched w/replacement

#Unmatched
describeBy(dta_Shp_st$terrai_are, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$Pop_1990, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$MeanL_1995, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$MaxL_1995, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$MeanT_1995, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$MeanP_1995, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$Slope, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$Elevation, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$Riv_Dist, dta_Shp_st$TrtBin)
describeBy(dta_Shp_st$Road_dist, dta_Shp_st$TrtBin)

#Matched w/out replacement
describeBy(psm_Pairs_st$terrai_are, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$Pop_1990, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$MeanL_1995, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$MaxL_1995, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$MeanT_1995, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$MeanP_1995, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$Slope, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$Elevation, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$Riv_Dist, psm_Pairs_st$TrtBin)
describeBy(psm_Pairs_st$Road_dist, psm_Pairs_st$TrtBin)

#Matched w/replacement
describeBy(model_data_st$terrai_are, model_data_st$TrtBin)
describeBy(model_data_st$Pop_1990, model_data_st$TrtBin)
describeBy(model_data_st$MeanL_1995, model_data_st$TrtBin)
describeBy(model_data_st$MaxL_1995, model_data_st$TrtBin)
describeBy(model_data_st$MeanT_1995, model_data_st$TrtBin)
describeBy(model_data_st$MeanP_1995, model_data_st$TrtBin)
describeBy(model_data_st$Slope, model_data_st$TrtBin)
describeBy(model_data_st$Elevation, model_data_st$TrtBin)
describeBy(model_data_st$Riv_Dist, model_data_st$TrtBin)
describeBy(model_data_st$Road_dist, model_data_st$TrtBin)







