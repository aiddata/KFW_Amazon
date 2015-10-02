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

## Early vs. Late

#analyticModelEarly1, no pair FE, no covars, 1995-2001
summary(analyticModelEarly1 <- lm(NDVILevelChange_95_01 ~ TrtBin, data=psm_Pairs))
#Standardized Betas
summary(analyticModelEarly1B <- lm(NDVILevelChange_95_01 ~ TrtBin, data=psm_PairsB))

#analyticModelEarly2, treatment effect + pair fixed effects, 1995-2001
analyticModelEarly2 <- "NDVILevelChange_95_01 ~ TrtBin + factor(PSM_match_ID)"

OutputEarly2=Stage2PSM(analyticModelEarly2,psm_Pairs,type="lm",table_out=TRUE)

#analyticModelEarly3, treatment effect + pair fixed effects + covars 1995-2001

#create new dataset and rename column names in new dataset to enable multiple columns in stargazer
Data_Early3 <- psm_Pairs
colnames(Data_Early3@data)[(colnames(Data_Early3@data)=="Pop_1990")] <- "Pop_B"
colnames(Data_Early3@data)[(colnames(Data_Early3@data)=="MeanT_1995")] <- "MeanT_B"
colnames(Data_Early3@data)[(colnames(Data_Early3@data)=="MeanP_1995")] <- "MeanP_B"
colnames(Data_Early3@data)[(colnames(Data_Early3@data)=="post_trend_temp_95_01")] <- "post_trend_temp"
colnames(Data_Early3@data)[(colnames(Data_Early3@data)=="post_trend_precip_95_01")] <- "post_trend_precip"
#colnames(Data_Early3@data)

analyticModelEarly3 <- "NDVILevelChange_95_01 ~ TrtBin+ pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"
OutputEarly3=Stage2PSM(analyticModelEarly3,Data_Early3,type="lm",table_out=TRUE)


#analyticalModelEver4, pair FEs, include enforcement years total as covar, 1995-2010
Data_Early4 <- psm_Pairs
colnames(Data_Early4@data)[(colnames(Data_Early4@data)=="Pop_1990")] <- "Pop_B"
colnames(Data_Early4@data)[(colnames(Data_Early4@data)=="MeanT_1995")] <- "MeanT_B"
colnames(Data_Early4@data)[(colnames(Data_Early4@data)=="MeanP_1995")] <- "MeanP_B"
colnames(Data_Early4@data)[(colnames(Data_Early4@data)=="post_trend_temp_mean")] <- "post_trend_temp"
colnames(Data_Early@data)[(colnames(Data_Early4@data)=="post_trend_precip_mean")] <- "post_trend_precip"

analyticModelEarly4 <- "NDVILevelChange_95_10 ~ TrtBin + enforce_to + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID)"

OutputEarly4=Stage2PSM(analyticModelEarly4,Data_Early4,type="lm",table_out=TRUE)

#________________________________________________________________________

#Create high pressure variable from pre_trend_NDVI
temp <- fivenum(Data_Early3$pre_trend_NDVI_Max)

#Categorical
temp_median <- fivenum(Data_Early3@data$pre_trend_NDVI_max)[3]
temp_categorical <- ifelse(Data_Early3@data$pre_trend_NDVI_max > temp_median, 1, 0)
length(which(temp_categorical<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Pretrend NDVI max continuous and treatment]
Data_Early3@data$pre_trend_NDVI_max_int <- Data_Early3@data$pre_trend_NDVI_max*Data_Early@data$TrtBin
#Interaction term [Pretrend NDVI max categorical and treatment]
pre_trend_NDVI_max_cat_int <- temp_categorical*Data_Early3@data$TrtBin

#_______________________________________________________________________


analyticModelEarly_5 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + pre_trend_NDVI_max_int"

OutputEarly_5=Stage2PSM(analyticModelEarly_5,Data_Early3,type="lm",table_out=TRUE)

analyticModelEarly_6 <- "NDVILevelChange_95_10 ~ TrtBin + temp_categorical + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + pre_trend_NDVI_max_cat_int"

OutputEarly_6=Stage2PSM(analyticModelEarly_6,Data_Early3,type="lm",table_out=TRUE)

#-----------------------------------------------------------------------
#Predicted pre trends
dta_Shp2 <- subset(dta_Shp, pre_trend_NDVI_max<0)
dta_Shp3 <- dta_Shp
dta_Shp3$BinNDVI=0
dta_Shp3$BinNDVI[dta_Shp3$pre_trend_NDVI_max<0]=1

summary(pressure_model_early_max <- lm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + MeanT_1995 + MeanP_1995 + pre_trend_temp_mean + 
pre_trend_temp_min + pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_max + pre_trend_precip_mean + 
Slope + Elevation + Riv_Dist + Road_dist + AvgD_FedCU + AvgD_StaCU + AvgD_Log + AvgD_Rail + 
AvgD_Mine + AvgD_City, data=dta_Shp)) 
summary(BinNDVI <- glm(pre_trend_NDVI_max ~ terrai_are + Pop_1995 + MeanT_1995 + MeanP_1995 + pre_trend_temp_mean + 
                                         pre_trend_temp_min + pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_max + pre_trend_precip_mean + 
                                         Slope + Elevation + Riv_Dist + Road_dist + AvgD_FedCU + AvgD_StaCU + AvgD_Log + AvgD_Rail + 
                                         AvgD_Mine + AvgD_City, data=dta_Shp2)) 

#+ AvgDistanceToMajorCities"



Output_1=Stage2PSM(pressure_model_early_max,Data_Early3,type="lm",table_out=TRUE)

pre_coefs_max_early <- coef(Output_1$standardized)
pre_coefs_max_early_1 <- coef(Output_1$unstandardized)

pre_coefs_max_early_1

#Create unstandardized predicted variables based on coefficients

model_int_early_1 <- pre_coefs_max_early_1[1]
model_int_early_2 <- pre_coefs_max_early_1[2] * Data_Early3@data$terrai_are
model_int_early_3 <- pre_coefs_max_early_1[3] * Data_Early3@data$Pop_B
model_int_early_4 <- pre_coefs_max_early_1[4] * Data_Early3@data$MeanT_B
model_int_early_5 <- pre_coefs_max_early_1[5] * Data_Early3@data$MeanP_B
model_int_early_6 <- pre_coefs_max_early_1[6] * Data_Early3@data$pre_trend_temp_mean
model_int_early_7 <- pre_coefs_max_early_1[7] * Data_Early3@data$pre_trend_temp_min
model_int_early_8 <- pre_coefs_max_early_1[8] * Data_Early3@data$pre_trend_temp_max
model_int_early_9 <- pre_coefs_max_early_1[9] * Data_Early3@data$pre_trend_precip_min
model_int_early_10 <- pre_coefs_max_early_1[10] * Data_Early3@data$pre_trend_precip_max
model_int_early_11 <- pre_coefs_max_early_1[11] * Data_Early3@data$pre_trend_precip_mean
model_int_early_12 <- pre_coefs_max_early_1[12] * Data_Early3@data$Slope
model_int_early_13 <- pre_coefs_max_early_1[13] * Data_Early3@data$Elevation
model_int_early_14 <- pre_coefs_max_early_1[14] * Data_Early3@data$Riv_Dist
model_int_early_15 <- pre_coefs_max_early_1[15] * Data_Early3@data$Road_dist
model_int_early_16 <- pre_coefs_max_early_1[16] * Data_Early3@data$AvgD_FedCU
model_int_early_17 <- pre_coefs_max_early_1[17] * Data_Early3@data$AvgD_StaCU
model_int_early_18 <- pre_coefs_max_early_1[18] * Data_Early3@data$AvgD_Log
model_int_early_19 <- pre_coefs_max_early_1[19] * Data_Early3@data$AvgD_Rail
model_int_early_20 <- pre_coefs_max_early_1[20] * Data_Early3@data$AvgD_Mine
model_int_early_21 <- pre_coefs_max_early_1[21] * Data_Early3@data$AvgD_City
#model_int_early_22 <- pre_coefs_max_early_1[22] * Data_Early3@data$AvgD_MajCi

predict_NDVI_early_max <- model_int_early_1+model_int_early_2+model_int_early_3+model_int_early_4+model_int_early_5+model_int_early_6+
  model_int_early_7+model_int_early_8+model_int_early_9+model_int_early_10+model_int_early_11+model_int_early_12+
  model_int_early_13+model_int_early_14+model_int_early_15+model_int_early_16+model_int_early_17+model_int_early_18+model_int_early_19+
  model_int_early_20+model_int_early_21#+model_int_early_22

predict_NDVI_early_max

Data_Early3@data["predict_NDVI_early_max"] <- predict_NDVI_early_max

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#Run models with predicted NDVI as high pressure & categorical high pressure
#------------------------------------------------------------------------------------

#Create high pressure variable from predicted pre-trend NDVI
temp <- fivenum(Data_Early3@data$predict_NDVI_early_max)

#Categorical
temp_median <- fivenum(Data_Early3@data$predict_NDVI_early_max)[3]
predict_categorical_early <- ifelse(Data_Early3@data$predict_NDVI_early_max > temp_median, 1, 0)
length(which(temp_categorical<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Predicted pretrend NDVI max continuous and treatment]
Data_Early3@data$predict_NDVI_early_max_int <- Data_Early3@data$predict_NDVI_early_max*Data_Early3@data$TrtBin
#Interaction term [Predicted pretrend NDVI max categorical and treatment]
predict_NDVI_early_max_cat_int <- predict_categorical*Data_Early3@data$TrtBin

#_________________________________________________________________________________________

analyticModelEarly_7 <- "NDVILevelChange_95_10 ~ TrtBin + predict_NDVI_early_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + predict_NDVI_early_max_int"

OutputEarly_7=Stage2PSM(analyticModelEarly_7,Data_Early3,type="lm",table_out=TRUE)

analyticModelEarly_8 <- "NDVILevelChange_95_10 ~ TrtBin + predict_categorical_early + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + predict_NDVI_early_max_cat_int +factor(PSM_match_ID) "

OutputEarly_8=Stage2PSM(analyticModelEarly_8,Data_Early3,type="lm",table_out=TRUE)
#_________________________________________________________________________________________

#analyticModelLate, treatment effect + pair fixed effects + covars 2001-2010
#create new dataset and rename column names in new dataset to enable multiple columns in stargazer
Data_Late <- psm_Pairs
colnames(Data_Late@data)[(colnames(Data_Late@data)=="Pop_2000")] <- "Pop_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="MeanT_2001")] <- "MeanT_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="MeanP_2001")] <- "MeanP_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="post_trend_temp_01_10")] <- "post_trend_temp"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="post_trend_precip_01_10")] <- "post_trend_precip"
#colnames(Data_Late@data)

analyticModelLate <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"
OutputLate=Stage2PSM(analyticModelLate,Data_Late,type="lm",table_out=TRUE)

#analyticModelLate_Enf, treatment effect + pair fixed effects + covars + enforcement years covar, 2001-2010
Data_Late_Enf <- psm_Pairs
colnames(Data_Late_Enf@data)[(colnames(Data_Late_Enf@data)=="Pop_2000")] <- "Pop_B"
colnames(Data_Late_Enf@data)[(colnames(Data_Late_Enf@data)=="MeanT_2001")] <- "MeanT_B"
colnames(Data_Late_Enf@data)[(colnames(Data_Late_Enf@data)=="MeanP_2001")] <- "MeanP_B"
colnames(Data_Late_Enf@data)[(colnames(Data_Late_Enf@data)=="post_trend_temp_01_10")] <- "post_trend_temp"
colnames(Data_Late_Enf@data)[(colnames(Data_Late_Enf@data)=="post_trend_precip_01_10")] <- "post_trend_precip"

analyticModelLate_Enf <- "NDVILevelChange_01_10 ~ TrtBin + enforce_to + pre_trend_NDVI_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"

OutputLate_Enf=Stage2PSM(analyticModelLate_Enf,Data_Late_Enf,type="lm",table_out=TRUE)

#___________________________________________________________________________________________

#Create high pressure variable from pre_trend_NDVI
temp <- fivenum(Data_Late$pre_trend_NDVI_Max)

#Categorical
temp_median_late <- fivenum(Data_Late@data$pre_trend_NDVI_max)[3]
temp_categorical_late <- ifelse(Data_Late@data$pre_trend_NDVI_max > temp_median_late, 1, 0)
length(which(temp_categorical_late<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Pretrend NDVI max continuous and treatment]
Data_Late@data$pre_trend_NDVI_max_int <- Data_Late@data$pre_trend_NDVI_max*Data_Late@data$TrtBin
#Interaction term [Pretrend NDVI max categorical and treatment]
pre_trend_NDVI_max_cat_int <- temp_categorical_late*Data_Late@data$TrtBin

#___________________________________________________________________________________________


analyticModelLate1 <- "NDVILevelChange_01_10 ~ TrtBin + enforce_to + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + pre_trend_NDVI_max_int"

OutputLate1=Stage2PSM(analyticModelLate1,Data_Late,type="lm",table_out=TRUE)

analyticModelLate2 <- "NDVILevelChange_95_10 ~ TrtBin + temp_categorical_late + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + pre_trend_NDVI_max_cat_int"

OutputLate2=Stage2PSM(analyticModelLate2,Data_Late,type="lm",table_out=TRUE)

#________________________________________________________________________

#_________________________________________________________________________________________
#Predicted pre trends
pressure_model_late_max <- "pre_trend_NDVI_max ~ terrai_are + Pop_B + MeanT_B + MeanP_B + pre_trend_temp_mean + 
pre_trend_temp_min + pre_trend_temp_max + pre_trend_precip_min + pre_trend_precip_max + pre_trend_precip_mean + 
Slope + Elevation + Riv_Dist + Road_dist + AvgD_FedCU + AvgD_StaCU + AvgD_Log + AvgD_Rail + 
AvgD_Mine + AvgD_City" #+ AvgDistanceToMajorCities"

Output_1=Stage2PSM(pressure_model_late_max,Data_Late,type="lm",table_out=TRUE)

pre_coefs_max_late <- coef(Output_1$standardized)
pre_coefs_max_late_1 <- coef(Output_1$unstandardized)

pre_coefs_max_late_1

#Create unstandardized predicted variables based on coefficients

model_int_late_1 <- pre_coefs_max_late_1[1]
model_int_late_2 <- pre_coefs_max_late_1[2] * Data_Late@data$terrai_are
model_int_late_3 <- pre_coefs_max_late_1[3] * Data_Late@data$Pop_B
model_int_late_4 <- pre_coefs_max_late_1[4] * Data_Late@data$MeanT_B
model_int_late_5 <- pre_coefs_max_late_1[5] * Data_Late@data$MeanP_B
model_int_late_6 <- pre_coefs_max_late_1[6] * Data_Late@data$pre_trend_temp_mean
model_int_late_7 <- pre_coefs_max_late_1[7] * Data_Late@data$pre_trend_temp_min
model_int_late_8 <- pre_coefs_max_late_1[8] * Data_Late@data$pre_trend_temp_max
model_int_late_9 <- pre_coefs_max_late_1[9] * Data_Late@data$pre_trend_precip_min
model_int_late_10 <- pre_coefs_max_late_1[10] * Data_Late@data$pre_trend_precip_max
model_int_late_11 <- pre_coefs_max_late_1[11] * Data_Late@data$pre_trend_precip_mean
model_int_late_12 <- pre_coefs_max_late_1[12] * Data_Late@data$Slope
model_int_late_13 <- pre_coefs_max_late_1[13] * Data_Late@data$Elevation
model_int_late_14 <- pre_coefs_max_late_1[14] * Data_Late@data$Riv_Dist
model_int_late_15 <- pre_coefs_max_late_1[15] * Data_Late@data$Road_dist
model_int_late_16 <- pre_coefs_max_late_1[16] * Data_Late@data$AvgD_FedCU
model_int_late_17 <- pre_coefs_max_late_1[17] * Data_Late@data$AvgD_StaCU
model_int_late_18 <- pre_coefs_max_late_1[18] * Data_Late@data$AvgD_Log
model_int_late_19 <- pre_coefs_max_late_1[19] * Data_Late@data$AvgD_Rail
model_int_late_20 <- pre_coefs_max_late_1[20] * Data_Late@data$AvgD_Mine
model_int_late_21 <- pre_coefs_max_late_1[21] * Data_Late@data$AvgD_City
#model_int_late_22 <- pre_coefs_max_late_1[22] * Data_Late@data$AvgD_MajCi

predict_NDVI_late_max <- model_int_late_1+model_int_late_2+model_int_late_3+model_int_late_4+
  model_int_late_5+model_int_late_6+model_int_late_7+model_int_late_8+model_int_late_9+
  model_int_late_10+model_int_late_11+model_int_late_12+model_int_late_13+model_int_late_14+
  model_int_late_15+model_int_late_16+model_int_late_17+model_int_late_18+model_int_late_19+
  model_int_late_20+model_int_late_21#+model_int_late_22

predict_NDVI_late_max

Data_Late@data["predict_NDVI_late_max"] <- predict_NDVI_late_max

#-------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------
#Run models with predicted NDVI as high pressure & categorical high pressure
#------------------------------------------------------------------------------------

#Create high pressure variable from predicted pre-trend NDVI
temp_late <- fivenum(Data_Late@data$predict_NDVI_late_max)

#Categorical
temp_median_late <- fivenum(Data_Late@data$predict_NDVI_late_max)[3]
predict_categorical_late <- ifelse(Data_Late@data$predict_NDVI_late_max > temp_median_late, 1, 0)
length(which(predict_categorical_late<1))

#variable containing all data in regions with pretrend NDVI above the median
#temp_new <- dta_Shp@data[dta_Shp@data$pre_trend_NDVI_max>temp_median,]

#Interaction term [ Predicted pretrend NDVI max continuous and treatment]
Data_Late@data$predict_NDVI_late_max_int <- Data_Late@data$predict_NDVI_late_max*Data_Late@data$TrtBin
#Interaction term [Predicted pretrend NDVI max categorical and treatment]
predict_NDVI_late_max_cat_int <- predict_categorical_late*Data_Late@data$TrtBin

#_________________________________________________________________________________________

analyticModelLate3 <- "NDVILevelChange_95_10 ~ TrtBin + predict_NDVI_late_max + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + factor(PSM_match_ID) + predict_NDVI_late_max_int"

OutputLate3=Stage2PSM(analyticModelLate3,Data_Late,type="lm",table_out=TRUE)

analyticModelLate4 <- "NDVILevelChange_95_10 ~ TrtBin + predict_categorical_late + MaxL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation  + Riv_Dist + Road_dist + predict_NDVI_late_max_cat_int +factor(PSM_match_ID) "

OutputLate4=Stage2PSM(analyticModelLate4,Data_Late,type="lm",table_out=TRUE)

#-------------------------------------------------------------------------------------
#without enforcement years
stargazer(OutputEarly2$standardized,OutputEarly3$standardized,OutputLate$standardized,
          keep=c("TrtBin", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_B", "MeanT_B","post_trend_temp","MeanP_B",
                 "post_trend_precip", "Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Pre-Trend NDVI", "Baseline NDVI","Area (hectares)", "Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Baseline Precipitation", "Precipitation Trends",
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2001 "," Max NDVI 2001-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)

#with enforcement years
stargazer(OutputEarly2$standardized,OutputEarly3$standardized,OutputLate$standardized,OutputLate_Enf$standardized,
          keep=c("TrtBin", "enforce_to", "pre_trend_NDVI_max","MaxL_1995", "terrai_are","Pop_B", "MeanT_B","post_trend_temp","MeanP_B",
                 "post_trend_precip", "Slope","Elevation","Riv_Dist","Road_dist"),
          covariate.labels=c("Treatment", "Enforcement", "Pre-Trend NDVI", "Baseline NDVI","Area (hectares)", "Baseline Population Density",
                             "Baseline Temperature", "Temperature Trends", "Baseline Precipitation", "Precipitation Trends",
                             "Slope", "Elevation", "Distance to River", "Distance to Road"),
          dep.var.labels=c("Max NDVI 1995-2001 "," Max NDVI 2001-2010"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)

#---------------------------------
#---------------------------------
# Tabulating Descriptive Statistics for Treatment and Control Groups, pre- and post-balance
#---------------------------------
#---------------------------------

#Using dta_Shp for the full dataset, no TrtBin

summary(dta_Shp$terrai_are)
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




