#-------------------------------------------------
#-------------------------------------------------
#Cross Sectional Models - KFW
#Testing in Cross Section the impact of being treated AFTER April 2001
#On the Mean Level of NDVI, measured as a change in the level of NDVI between start and end year (1995-2001, 2001-2010)
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
dta_Shp$NDVIslope_95_10 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1995,2010,"SP_ID")
dta_Shp@data["NDVILevelChange_95_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_1995

#NDVI Trends for 1995-2001
dta_Shp$post_trend_NDVI_95_01 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
dta_Shp@data["NDVILevelChange_95_01"] <- dta_Shp$MeanL_2001 - dta_Shp$MeanL_1995
#dta_Shp@data["NDVIslopeChange_95_01"] <- dta_Shp@data["post_trend_NDVI_95_01"] - dta_Shp@data["pre_trend_NDVI_mean"]

#NDVI Trends for 2001-2010
dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")
dta_Shp@data["NDVILevelChange_01_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_2001
#dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp@data["post_trend_NDVI_01_10"] - dta_Shp@data["pre_trend_NDVI_mean"]

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

analyticModelEarly3 <- "NDVILevelChange_95_01 ~ TrtBin+ pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"
OutputEarly3=Stage2PSM(analyticModelEarly3,Data_Early3,type="lm",table_out=TRUE)


##INTERACTION TERMS EARLY
#Interaction pretrend NDVI & treatment
Data_Early3@data["pretrend_Trt"] <- ((Data_Early3@data$pre_trend_NDVI_mean)*(Data_Early3@data$TrtBin))

#Interaction MeanL1995 & treatment
Data_Early3@data["MeanL1995_Trt"] <- ((Data_Early3@data$MeanL_1995)*(Data_Early3@data$TrtBin))

#Interaction area & treatment
Data_Early3@data["area_Trt"] <- ((Data_Early3@data$terrai_are)*(Data_Early3@data$TrtBin))

#Interaction Population & treatment
Data_Early3@data["Pop_Trt"] <- ((Data_Early3@data$Pop_B)*(Data_Early3@data$TrtBin))

#Interaction MeanT & treatment
Data_Early3@data["MeanT_Trt"] <- ((Data_Early3@data$MeanT_B)*(Data_Early3@data$TrtBin))

#Interaction postT & treatment
Data_Early3@data["postT_Trt"] <- ((Data_Early3@data$post_trend_temp)*(Data_Early3@data$TrtBin))

#Interaction MeanP & Treatment
Data_Early3@data["MeanP_Trt"] <- ((Data_Early3@data$MeanP_B)*(Data_Early3@data$TrtBin))

#Interaction posttrend & Treatment
Data_Early3@data["posttrendP_Trt"] <- ((Data_Early3@data$post_trend_precip)*(Data_Early3@data$TrtBin))

#Interaction Slope & Treatment
Data_Early3@data["Slope_Trt"] <- ((Data_Early3@data$Slope)*(Data_Early3@data$TrtBin))

#Interaction Elevation & treatment
Data_Early3@data["Elevation_Trt"] <- ((Data_Early3@data$Elevation)*(Data_Early3@data$TrtBin))

#Interaction Riv Dist & treatment
Data_Early3@data["RivDist_Trt"] <- ((Data_Early3@data$Riv_Dist)*(Data_Early3@data$TrtBin))

#Interaction Road Distance & treatment
Data_Early3@data["RoadDist_Trt"] <- ((Data_Early3@data$Road_dist)*(Data_Early3@data$TrtBin))


## MODELS WITH INTERACTIONS
analyticModelEarly_1 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + pretrend_Trt"

analyticModelEarly_2 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanL1995_Trt"

analyticModelEarly_3 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + area_Trt"

analyticModelEarly_4 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Pop_Trt"

analyticModelEarly_5 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanT_Trt"

analyticModelEarly_6 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + postT_Trt"

analyticModelEarly_7 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanP_Trt"

analyticModelEarly_8 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + posttrendP_Trt"

analyticModelEarly_9 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Slope_Trt"

analyticModelEarly_10 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Elevation_Trt"

analyticModelEarly_11 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + RivDist_Trt"

analyticModelEarly_12 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + RoadDist_Trt"

OutputEarly3_1=Stage2PSM(analyticModelEarly_1,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_2=Stage2PSM(analyticModelEarly_2,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_3=Stage2PSM(analyticModelEarly_3,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_4=Stage2PSM(analyticModelEarly_4,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_5=Stage2PSM(analyticModelEarly_5,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_6=Stage2PSM(analyticModelEarly_6,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_7=Stage2PSM(analyticModelEarly_7,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_8=Stage2PSM(analyticModelEarly_8,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_9=Stage2PSM(analyticModelEarly_9,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_10=Stage2PSM(analyticModelEarly_10,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_11=Stage2PSM(analyticModelEarly_11,Data_Early3,type="lm",table_out=TRUE)
OutputEarly3_12=Stage2PSM(analyticModelEarly_12,Data_Early3,type="lm",table_out=TRUE)

stargazer(OutputEarly3$standardized, OutputEarly3_1$standardized, 
          OutputEarly3_2$standardized, OutputEarly3_3$standardized, OutputEarly3_4$standardized, 
          OutputEarly3_5$standardized, OutputEarly3_6$standardized, OutputEarly3_7$standardized, 
          OutputEarly3_8$standardized, OutputEarly3_9$standardized, OutputEarly3_10$standardized, 
          OutputEarly3_11$standardized, OutputEarly3_12$standardized, 
          keep=c("TrtBin", "pretrend_Trt", "MeanL1995_Trt", "area_Trt", "Pop_Trt", 
                 "MeanT_Trt", "postT_Trt", "MeanP_Trt", "posttrendP_Trt", "Slope_Trt", 
                 "Elevation_Trt", "RivDist_Trt", "RoadDist_Trt"),
          covariate.labels=c("Treatment", "Pretrend Interaction", "MeanL1995 Interaction", 
                             "Area Interaction", "Pop Interaction", "MeanT Interaction", "PostT Interaction", "MeanP INteraction", 
                             "PosttrendP Interaction", "Slope Interaction", "Elevation Interaction", "RivDist Int", "RoadDist Int"),
          dep.var.labels=c("Mean NDVI 1995-2001"),
          title="Regression Results", type="html", omit.stat=c("f","ser"), align=TRUE)


#analyticModelLate, treatment effect + pair fixed effects + covars 2001-2010
#create new dataset and rename column names in new dataset to enable multiple columns in stargazer
Data_Late <- psm_Pairs
colnames(Data_Late@data)[(colnames(Data_Late@data)=="Pop_2000")] <- "Pop_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="MeanT_2001")] <- "MeanT_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="MeanP_2001")] <- "MeanP_B"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="post_trend_temp_01_10")] <- "post_trend_temp"
colnames(Data_Late@data)[(colnames(Data_Late@data)=="post_trend_precip_01_10")] <- "post_trend_precip"
#colnames(Data_Late@data)

analyticModelLate <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"
OutputLate=Stage2PSM(analyticModelLate,Data_Late,type="lm",table_out=TRUE)

##LATE MODEL INTERACTIONS

##INTERACTION TERMS
#Interaction pretrend NDVI & treatment
Data_Late@data["pretrend_Trt_1"] <- ((Data_Late@data$pre_trend_NDVI_mean)*(Data_Late@data$TrtBin))

#Interaction MeanL1995 & treatment
Data_Late@data["MeanL1995_Trt_1"] <- ((Data_Late@data$MeanL_1995)*(Data_Late@data$TrtBin))

#Interaction area & treatment
Data_Late@data["area_Trt_1"] <- ((Data_Late@data$terrai_are)*(Data_Late@data$TrtBin))

#Interaction Population & treatment
Data_Late@data["Pop_Trt_1"] <- ((Data_Late@data$Pop_B)*(Data_Late@data$TrtBin))

#Interaction MeanT & treatment
Data_Late@data["MeanT_Trt_1"] <- ((Data_Late@data$MeanT_B)*(Data_Late@data$TrtBin))

#Interaction postT & treatment
Data_Late@data["postT_Trt_1"] <- ((Data_Late@data$post_trend_temp)*(Data_Late@data$TrtBin))

#Interaction MeanP & Treatment
Data_Late@data["MeanP_Trt_1"] <- ((Data_Late@data$MeanP_B)*(Data_Late@data$TrtBin))

#Interaction posttrend precipitation & Treatment
Data_Late@data["posttrendP_Trt_1"] <- ((Data_Late@data$post_trend_precip)*(Data_Late@data$TrtBin))

#Interaction Slope & Treatment
Data_Late@data["Slope_Trt_1"] <- ((Data_Late@data$Slope)*(Data_Late@data$TrtBin))

#Interaction Elevation & treatment
Data_Late@data["Elevation_Trt_1"] <- ((Data_Late@data$Elevation)*(Data_Late@data$TrtBin))

#Interaction Riv Dist & treatment
Data_Late@data["RivDist_Trt_1"] <- ((Data_Late@data$Riv_Dist)*(Data_Late@data$TrtBin))

#Interaction Road Distance & treatment
Data_Late@data["RoadDist_Trt_1"] <- ((Data_Late@data$Road_dist)*(Data_Late@data$TrtBin))

##MODELS WITH INTERACTIONS

analyticModelLate_1 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + pretrend_Trt_1"

analyticModelLate_2 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanL1995_Trt_1"

analyticModelLate_3 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + area_Trt_1"

analyticModelLate_4 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Pop_Trt_1"

analyticModelLate_5 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanT_Trt_1"

analyticModelLate_6 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + postT_Trt_1"

analyticModelLate_7 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + MeanP_Trt_1"

analyticModelLate_8 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + posttrendP_Trt_1"

analyticModelLate_9 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Slope_Trt_1"

analyticModelLate_10 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + Elevation_Trt_1"

analyticModelLate_11 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + RivDist_Trt_1"

analyticModelLate_12 <- "NDVILevelChange_01_10 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B + post_trend_temp + 
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID) + RoadDist_Trt_1"

OutputLate_1=Stage2PSM(analyticModelLate_1,Data_Late,type="lm",table_out=TRUE)
OutputLate_2=Stage2PSM(analyticModelLate_2,Data_Late,type="lm",table_out=TRUE)
OutputLate_3=Stage2PSM(analyticModelLate_3,Data_Late,type="lm",table_out=TRUE)
OutputLate_4=Stage2PSM(analyticModelLate_4,Data_Late,type="lm",table_out=TRUE)
OutputLate_5=Stage2PSM(analyticModelLate_5,Data_Late,type="lm",table_out=TRUE)
OutputLate_6=Stage2PSM(analyticModelLate_6,Data_Late,type="lm",table_out=TRUE)
OutputLate_7=Stage2PSM(analyticModelLate_7,Data_Late,type="lm",table_out=TRUE)
OutputLate_8=Stage2PSM(analyticModelLate_8,Data_Late,type="lm",table_out=TRUE)
OutputLate_9=Stage2PSM(analyticModelLate_9,Data_Late,type="lm",table_out=TRUE)
OutputLate_10=Stage2PSM(analyticModelLate_10,Data_Late,type="lm",table_out=TRUE)
OutputLate_11=Stage2PSM(analyticModelLate_11,Data_Late,type="lm",table_out=TRUE)
OutputLate_12=Stage2PSM(analyticModelLate_12,Data_Late,type="lm",table_out=TRUE)

stargazer(OutputLate$standardized,OutputLate_1$standardized,OutputLate_2$standardized,
          OutputLate_3$standardized,OutputLate_4$standardized,OutputLate_5$standardized,
          OutputLate_6$standardized,OutputLate_7$standardized,OutputLate_8$standardized,
          OutputLate_9$standardized,OutputLate_10$standardized,OutputLate_11$standardized,
          OutputLate_12$standardized,
          keep=c("TrtBin", "pretrend_Trt_1","MeanL1995_Trt_1", "area_Trt_1","Pop_Trt_1", "MeanT_Trt_1","postT_Trt_1","MeanP_Trt_1",
                 "posttrendP_Trt_1", "Slope_Trt_1","Elevation_Trt_1","RivDist_Trt_1","RoadDist_Trt_1"),
          covariate.labels=c("Treatment", "PretrendMean Interaction", "MeanL1995 Interaction", 
                             "Area Interaction", "Pop Interaction,", "MeanT Interaction", "PostT Interaction", "MeanP INteraction", 
                             "PosttrendP Interaction", "Slope Interaction", "Elevation Interaction", "RivDist Int", "RoadDist Int"),
          dep.var.labels=c(" Mean NDVI 2001-2010"),
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

