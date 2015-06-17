#-------------------------------------------------
#-------------------------------------------------
#Cross Sectional Models - KFW
#Testing in Cross Section the impact of being treated EVER vs. not being treated
#On the Mean Level of NDVI, measured as a change in the level of NDVI between start and end year (1995-2010)
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
dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MaxL_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")
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

#-------------------------------------------------
#-------------------------------------------------
#Define and run the first-stage of the PSM, calculating propensity scores
#-------------------------------------------------
#-------------------------------------------------
psmModel <-  "TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
pre_trend_NDVI_mean + pre_trend_NDVI_max + Slope + Elevation +  MeanL_1995 + MaxL_1995 + Riv_Dist + Road_dist +
pre_trend_precip_mean + pre_trend_precip_max"

psmRes <- SAT::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel,drop="support",visual=TRUE)

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



#-------------------------------------------------
#-------------------------------------------------
#Cross-section Models
#-------------------------------------------------
#-------------------------------------------------

#Scale all of the data to get standardized coefficients, create psm_PairsB
psm_PairsB <- psm_Pairs
ind <- sapply(psm_PairsB@data, is.numeric)
psm_PairsB@data[ind] <- lapply(psm_PairsB@data[ind],scale)


#Ever vs. Never

#OLS, no pair FEs, no covars, 1995-2010

summary(analyticModelEver1 <- lm(NDVILevelChange_95_10 ~ TrtBin, data=psm_Pairs))
summary(analyticModelEver1B <- lm(NDVILevelChange_95_10 ~ TrtBin, data=psm_PairsB))

#analyticModelEver2, pair FEs, no covars, 1995-2010

analyticModelEver2 <- "NDVILevelChange_95_10 ~ TrtBin + factor(PSM_match_ID)"

OutputEver2=Stage2PSM(analyticModelEver2,psm_Pairs,type="lm",table_out=TRUE)

#analyticModelEver3, pair FEs, covars, 1995-2010

#create new dataset and rename column names in new dataset to enable multiple columns in stargazer
Data_Ever3 <- psm_Pairs
colnames(Data_Ever3@data)[(colnames(Data_Ever3@data)=="Pop_1990")] <- "Pop_B"
colnames(Data_Ever3@data)[(colnames(Data_Ever3@data)=="MeanT_1995")] <- "MeanT_B"
colnames(Data_Ever3@data)[(colnames(Data_Ever3@data)=="MeanP_1995")] <- "MeanP_B"
colnames(Data_Ever3@data)[(colnames(Data_Ever3@data)=="post_trend_temp_mean")] <- "post_trend_temp"
colnames(Data_Ever3@data)[(colnames(Data_Ever3@data)=="post_trend_precip_mean")] <- "post_trend_precip"
#colnames(Data_Ever3@data)

analyticModelEver3 <- "NDVILevelChange_95_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID)"

OutputEver3=Stage2PSM(analyticModelEver3,Data_Ever3,type="lm",table_out=TRUE)

stargazer(OutputEver2$standardized, OutputEver3$standardized,
          keep=c("TrtBin","pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_B","MeanT_B","post_trend_temp","MeanP_B",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI","Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Baseline Precipitation", "Precipitation Trends",
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)