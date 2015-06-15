analyticModelEarly3 <- "NDVILevelChange_95_01 ~ TrtBin + pre_trend_NDVI_mean + MeanL_1995 + terrai_are + Pop_B + MeanT_B  + post_trend_temp +
MeanP_B + post_trend_precip + Slope + Elevation + Riv_Dist + Road_dist + factor(PSM_match_ID)"
OutputEarly3=Stage2PSM(analyticModelEarly3,Data_Early3,type="lm",table_out=TRUE)

##INTERACTION TERMS
#Interaction terrain area & treatment
Data_Early3@data["pretrend_Trt"] <- ((Data_Early3@data$pre_trend_NDVI_mean)*(Data_Early3@data$TrtBin))

#Interaction Pop_ & treatment
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

OutputEarly_2=Stage2PSM(analyticModelEarly_12,Data_Early3,type="lm",table_out=TRUE)


##LATE MODEL INTERACTIONS

##INTERACTION TERMS
#Interaction terrain area & treatment
Data_Late@data["pretrend_Trt_1"] <- ((Data_Late@data$pre_trend_NDVI_mean)*(Data_Late@data$TrtBin))

#Interaction Pop_ & treatment
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

#Interaction posttrend & Treatment
Data_Late@data["posttrendP_Trt_1"] <- ((Data_Late@data$post_trend_precip)*(Data_Late@data$TrtBin))

#Interaction Slope & Treatment
Data_Late@data["Slope_Trt_1"] <- ((Data_Late@data$Slope)*(Data_Late@data$TrtBin))

#Interaction Elevation & treatment
Data_Late@data["Elevation_Trt_1"] <- ((Data_Late@data$Elevation)*(Data_Late@data$TrtBin))

#Interaction Riv Dist & treatment
Data_Late@data["RivDist_Trt_1"] <- ((Data_Late@data$Riv_Dist)*(Data_Late@data$TrtBin))

#Interaction Road Distance & treatment
Data_Late@data["RoadDist_Trt_1"] <- ((Data_Late@data$Road_dist)*(Data_Late@data$TrtBin))

######

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

OutputLate=Stage2PSM(analyticModelLate_12,Data_Late,type="lm",table_out=TRUE)
