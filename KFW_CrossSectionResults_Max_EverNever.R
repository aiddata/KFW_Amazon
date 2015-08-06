#-------------------------------------------------
#-------------------------------------------------
#Cross Sectional Models - KFW
#Testing in Cross Section the impact of being treated EVER vs. not being treated
#On the Mean Level of NDVI, measured as a change in the level of NDVI between start and end year (1995-2010)
#-------------------------------------------------
#-------------------------------------------------
library(devtools)
devtools::install_github("itpir/SCI@master")
library(SCI)
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

#analyticalModelEver4, pair FEs, include enforcement years total as covar, 1995-2010
Data_Ever4 <- psm_Pairs
colnames(Data_Ever4@data)[(colnames(Data_Ever4@data)=="Pop_1990")] <- "Pop_B"
colnames(Data_Ever4@data)[(colnames(Data_Ever4@data)=="MeanT_1995")] <- "MeanT_B"
colnames(Data_Ever4@data)[(colnames(Data_Ever4@data)=="MeanP_1995")] <- "MeanP_B"
colnames(Data_Ever4@data)[(colnames(Data_Ever4@data)=="post_trend_temp_mean")] <- "post_trend_temp"
colnames(Data_Ever4@data)[(colnames(Data_Ever4@data)=="post_trend_precip_mean")] <- "post_trend_precip"

analyticModelEver4 <- "NDVILevelChange_95_10 ~ TrtBin + enforce_to + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID)"

OutputEver4=Stage2PSM(analyticModelEver4,Data_Ever4,type="lm",table_out=TRUE)

#________________________________________________________________________

#Create high pressure variable from pre_trend_NDVI
temp <- fivenum(Data_Ever3@data$pre_trend_NDVI_max)

#Categorical
temp_median <- fivenum(Data_Ever3@data$pre_trend_NDVI_max)[3]
temp_categorical <- ifelse(Data_Ever3@data$pre_trend_NDVI_max > temp_median, 1, 0)
length(which(temp_categorical<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Pretrend NDVI max continuous and treatment]
Data_Ever3@data$pre_trend_NDVI_max_int <- Data_Ever3@data$pre_trend_NDVI_max*Data_Ever3@data$TrtBin
#Interaction term [Pretrend NDVI max categorical and treatment]
pre_trend_NDVI_max_cat_int <- temp_categorical*Data_Ever3@data$TrtBin

#_________________________________________________________________________________________

analyticModelEver5 <- "NDVILevelChange_95_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + pre_trend_NDVI_max_int"

OutputEver5=Stage2PSM(analyticModelEver5,Data_Ever3,type="lm",table_out=TRUE)

analyticModelEver6 <- "NDVILevelChange_95_10 ~ TrtBin + temp_categorical + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + pre_trend_NDVI_max_cat_int +factor(PSM_match_ID) "

OutputEver6=Stage2PSM(analyticModelEver6,Data_Ever3,type="lm",table_out=TRUE)

#_________________________________________________________________________________________
#Predicted pre trends
pressure_model_1_max <- "pre_trend_NDVI_max ~ terrai_are + Pop_B + MeanT_B + MeanP_B + pre_trend_temp_mean + 
pre_trend_temp_min + pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_max + pre_trend_precip_mean + 
Slope + Elevation + Riv_Dist + Road_dist + AvgD_FedCU + AvgD_StaCU + AvgD_Log + AvgD_Rail + 
AvgD_Mine + AvgD_City" #+ AvgDistanceToMajorCities"

Output_1=Stage2PSM(pressure_model_1_max,Data_Ever3,type="lm",table_out=TRUE)

pre_coefs_max <- coef(Output_1$standardized)
pre_coefs_max_1 <- coef(Output_1$unstandardized)

pre_coefs_max_1

#Create unstandardized predicted variables based on coefficients

model_int_1 <- pre_coefs_max_1[1]
model_int_2 <- pre_coefs_max_1[2] * Data_Ever3@data$terrai_are
model_int_3 <- pre_coefs_max_1[3] * Data_Ever3@data$Pop_B
model_int_4 <- pre_coefs_max_1[4] * Data_Ever3@data$MeanT_B
model_int_5 <- pre_coefs_max_1[5] * Data_Ever3@data$MeanP_B
model_int_6 <- pre_coefs_max_1[6] * Data_Ever3@data$pre_trend_temp_mean
model_int_7 <- pre_coefs_max_1[7] * Data_Ever3@data$pre_trend_temp_min
model_int_8 <- pre_coefs_max_1[8] * Data_Ever3@data$pre_trend_temp_max
model_int_9 <- pre_coefs_max_1[9] * Data_Ever3@data$pre_trend_precip_min
model_int_10 <- pre_coefs_max_1[10] * Data_Ever3@data$pre_trend_precip_max
model_int_11 <- pre_coefs_max_1[11] * Data_Ever3@data$pre_trend_precip_mean
model_int_12 <- pre_coefs_max_1[12] * Data_Ever3@data$Slope
model_int_13 <- pre_coefs_max_1[13] * Data_Ever3@data$Elevation
model_int_14 <- pre_coefs_max_1[14] * Data_Ever3@data$Riv_Dist
model_int_15 <- pre_coefs_max_1[15] * Data_Ever3@data$Road_dist
model_int_16 <- pre_coefs_max_1[16] * Data_Ever3@data$AvgD_FedCU
model_int_17 <- pre_coefs_max_1[17] * Data_Ever3@data$AvgD_StaCU
model_int_18 <- pre_coefs_max_1[18] * Data_Ever3@data$AvgD_Log
model_int_19 <- pre_coefs_max_1[19] * Data_Ever3@data$AvgD_Rail
model_int_20 <- pre_coefs_max_1[20] * Data_Ever3@data$AvgD_Mine
model_int_21 <- pre_coefs_max_1[21] * Data_Ever3@data$AvgD_City
#model_int_22 <- pre_trend_coefs_1[22] * Data_Ever3@data$AvgDistanceToMajorCities

predict_NDVI_max <- model_int_1+model_int_2+model_int_3+model_int_4+model_int_5+model_int_6+
  model_int_7+model_int_8+model_int_9+model_int_10+model_int_11+model_int_12+
  model_int_13+model_int_14+model_int_15+model_int_16+model_int_17+model_int_18+model_int_19+
  model_int_20+model_int_21#+model_int_22

predict_NDVI_max

Data_Ever3@data["predict_NDVI_max"] <- predict_NDVI_max

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#Run models with predicted NDVI as high pressure & categorical high pressure
#------------------------------------------------------------------------------------

#Create high pressure variable from predicted pre-trend NDVI
temp <- fivenum(Data_Ever3@data$predict_NDVI_max)

#Categorical
temp_median <- fivenum(Data_Ever3@data$predict_NDVI_max)[3]
predict_categorical <- ifelse(Data_Ever3@data$predict_NDVI_max > temp_median, 1, 0)
length(which(temp_categorical<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Predicted pretrend NDVI max continuous and treatment]
Data_Ever3@data$predict_NDVI_max_int <- Data_Ever3@data$predict_NDVI_max*Data_Ever3@data$TrtBin
#Interaction term [Predicted pretrend NDVI max categorical and treatment]
predict_NDVI_max_cat_int <- predict_categorical*Data_Ever3@data$TrtBin

#_________________________________________________________________________________________

analyticModelEver7 <- "NDVILevelChange_95_10 ~ TrtBin + predict_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + predict_NDVI_max_int"

OutputEver7=Stage2PSM(analyticModelEver7,Data_Ever3,type="lm",table_out=TRUE)

analyticModelEver8 <- "NDVILevelChange_95_10 ~ TrtBin + predict_categorical + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + predict_NDVI_max_cat_int +factor(PSM_match_ID) "

OutputEver8=Stage2PSM(analyticModelEver8,Data_Ever3,type="lm",table_out=TRUE)



#-------------------------------------------------------------------------------------



stargazer(OutputEver2$standardized, OutputEver3$standardized, OutputEver4$standardized,
          keep=c("TrtBin", "enforce_to", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_B","MeanT_B","post_trend_temp","MeanP_B",
                 "post_trend_precip","Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Enforcement Years", "Pre-Trend NDVI", "Baseline NDVI", "Area (hectares)","Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Baseline Precipitation", "Precipitation Trends",
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)

#---------------------------------
#---------------------------------
# Tabulating Descriptive Statistics for Treatment and Control Groups, pre- and post-balance
#---------------------------------
#---------------------------------

#Using dta_Shp for the full dataset, no TrtBin

mean(dta_Shp$terrai_are)
summary(dta_Shp$Pop_1990)
summary(dta_Shp$MeanL_1995)
summary(dta_Shp$MaxL_1995)
summary(dta_Shp$MeanT_1995)
summary(dta_Shp$MeanP_1995)
summary(dta_Shp$Slope)
summary(dta_Shp$Elevation)
summary(dta_Shp$Riv_Dist)
summary(dta_Shp$Road_dist)

#Using dta_Shp to get pre-matching, pre-paired data
describeBy(dta_Shp$terrai_are, dta_Shp$TrtBin)
describeBy(dta_Shp$Pop_1990, dta_Shp$TrtBin)
describeBy(dta_Shp$MeanL_1995, dta_Shp$TrtBin)
describeBy(dta_Shp$MaxL_1995, dta_Shp$TrtBin)
describeBy(dta_Shp$MeanT_1995, dta_Shp$TrtBin)
describeBy(dta_Shp$MeanP_1995, dta_Shp$TrtBin)
describeBy(dta_Shp$Slope, dta_Shp$TrtBin)
describeBy(dta_Shp$Elevation, dta_Shp$TrtBin)
describeBy(dta_Shp$Riv_Dist, dta_Shp$TrtBin)
describeBy(dta_Shp$Road_dist, dta_Shp$TrtBin)

#Using psm_pairs to get post-matching, post-paired data
describeBy(psm_Pairs$terrai_are, psm_Pairs$TrtBin)
describeBy(psm_Pairs$Pop_1990, psm_Pairs$TrtBin)
describeBy(psm_Pairs$MeanL_1995, psm_Pairs$TrtBin)
describeBy(psm_Pairs$MaxL_1995, psm_Pairs$TrtBin)
describeBy(psm_Pairs$MeanT_1995, psm_Pairs$TrtBin)
describeBy(psm_Pairs$MeanP_1995, psm_Pairs$TrtBin)
describeBy(psm_Pairs$Slope, psm_Pairs$TrtBin)
describeBy(psm_Pairs$Elevation, psm_Pairs$TrtBin)
describeBy(psm_Pairs$Riv_Dist, psm_Pairs$TrtBin)
describeBy(psm_Pairs$Road_dist, psm_Pairs$TrtBin)

