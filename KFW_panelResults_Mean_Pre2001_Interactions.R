#-------------------------------------------------
#-------------------------------------------------
#Panel Models - KFW
#Testing in Cross Section the impact of being treated BEFORE March, 2001
#On the Mean Level of NDVI, measured as the yearly mean NDVI value (LTDR)
#-------------------------------------------------
#-------------------------------------------------
library(devtools)
devtools::install_github("itpir/SAT@master")
library(SAT)
library(stargazer)
library(lmtest)
library(multiwayvcov)
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
dta_Shp@data["NDVIslopeChange_95_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_1995

#NDVI Trends for 1995-2001
dta_Shp$post_trend_NDVI_95_01 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",1995,2001,"SP_ID")
dta_Shp@data["NDVIslopeChange_95_01"] <- dta_Shp$MeanL_2001 - dta_Shp$MeanL_1995
#NDVI Trends for 2001-2010
dta_Shp$post_trend_NDVI_01_10 <- timeRangeTrend(dta_Shp,"MeanL_[0-9][0-9][0-9][0-9]",2001,2010,"SP_ID")
dta_Shp@data["NDVIslopeChange_01_10"] <- dta_Shp$MeanL_2010 - dta_Shp$MeanL_2001
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
#View(trttable)


#-------------------------------------------------
#-------------------------------------------------
#Convert from a wide-form dataset for the Cross-sectional 
#to a long-form dataset for the panel model.

varList = c("MeanL_","MaxL_")
psm_Long <- BuildTimeSeries(dta=psm_Pairs,idField="reu_id",varList_pre=varList,1982,2010,colYears=c("demend_y","enforce_st"),interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_","MaxP_","MinP_","MinT_"))
psm_Long$Year <- as.numeric(psm_Long$Year)

#Create post-2004 dummy
psm_Long$Post2004 <- 0
psm_Long$Post2004[psm_Long$Year >= 2004] <- 1

##INTERACTIONS

psm_Long["area_Trt"] <- ((psm_Long$terrai_are)*(psm_Long$TrtMnt_demend_y))
psm_Long["Pop_Trt"] <-((psm_Long$Pop_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MeanT_Trt"] <-((psm_Long$MeanT_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MeanP_Trt"] <- ((psm_Long$MeanP_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MaxT_Trt"] <- ((psm_Long$MaxT_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MaxP_Trt"] <- ((psm_Long$MaxP_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MinT_Trt"] <- ((psm_Long$MinT_)*(psm_Long$TrtMnt_demend_y))
psm_Long["MinP_Trt"] <- ((psm_Long$MinP_)*(psm_Long$TrtMnt_demend_y))
psm_Long["Year_Trt"] <- ((psm_Long$Year)*(psm_Long$TrtMnt_demend_y))
psm_Long["Slope_Trt"] <- ((psm_Long$Slope)*(psm_Long$TrtMnt_demend_y))
psm_Long["Elevation_Trt"] <- ((psm_Long$Elevation)*(psm_Long$TrtMnt_demend_y))
psm_Long["RivDist_Trt"] <- ((psm_Long$Riv_Dist)*(psm_Long$TrtMnt_demend_y))
psm_Long["RoadDist_Trt"] <- ((psm_Long$Road_dist)*(psm_Long$TrtMnt_demend_y))


pModelMax_A <- "MaxL_ ~ TrtMnt_demend_y + TrtMnt_enforce_st + factor(reu_id) "
pModelMax_B <- "MaxL_ ~ TrtMnt_demend_y + TrtMnt_enforce_st + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) "
pModelMax_C <- "MaxL_ ~ TrtMnt_demend_y + TrtMnt_enforce_st + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year"
pModelMax_D <- "MaxL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + factor(reu_id) + Year + Post2004 + Post2004*TrtMnt_demend_y + Post2004*TrtMnt_demend_y*Road_dist + Post2004*Road_dist"
pModelMax_E <- "MaxL_ ~ TrtMnt_demend_y + TrtMnt_enforce_st + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_ + factor(reu_id) + Year + Post2004 + Post2004*TrtMnt_demend_y + Post2004*TrtMnt_demend_y*Road_dist + Post2004*Road_dist"


pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_D_fit <- Stage2PSM(pModelMax_D ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_E_fit <- Stage2PSM(pModelMax_E ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
#-------------------------------------------------

##MODELS WITH INTERACTIONS
pModelMean_Int1 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + area_Trt"
pModelMean_Int2 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + Pop_Trt"
pModelMean_Int3 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MeanT_Trt"
pModelMean_Int4 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MeanP_Trt"
pModelMean_Int5 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MaxT_Trt"
pModelMean_Int6 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MaxP_Trt"
pModelMean_Int7 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MinT_Trt"
pModelMean_Int8 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + MinP_Trt"
pModelMean_Int9 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + Year_Trt"
pModelMean_Int10 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + Slope_Trt"
pModelMean_Int11 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + Elevation_Trt"
pModelMean_Int12 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + RivDist_Trt"
pModelMean_Int13 <- "MeanL_ ~ TrtMnt_demend_y + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year + RoadDist_Trt"

#####################

pModelMean_Int_fit1 <- Stage2PSM(pModelMean_Int1 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit2 <- Stage2PSM(pModelMean_Int2 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit3 <- Stage2PSM(pModelMean_Int3 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit4 <- Stage2PSM(pModelMean_Int4 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit5 <- Stage2PSM(pModelMean_Int5 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit6 <- Stage2PSM(pModelMean_Int6 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit7 <- Stage2PSM(pModelMean_Int7 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit8 <- Stage2PSM(pModelMean_Int8 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit9 <- Stage2PSM(pModelMean_Int9 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit10 <- Stage2PSM(pModelMean_Int10 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit11 <- Stage2PSM(pModelMean_Int11 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit12 <- Stage2PSM(pModelMean_Int12 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_Int_fit13 <- Stage2PSM(pModelMean_Int13 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))




stargazer(pModelMean_Int_fit1 $cmreg,pModelMean_Int_fit2 $cmreg,pModelMean_Int_fit3 $cmreg,
          pModelMean_Int_fit4 $cmreg,pModelMean_Int_fit5 $cmreg,pModelMean_Int_fit6 $cmreg,
          pModelMean_Int_fit7 $cmreg,pModelMean_Int_fit8 $cmreg,pModelMean_Int_fit9 $cmreg,
          pModelMean_Int_fit10 $cmreg,pModelMean_Int_fit11 $cmreg,pModelMean_Int_fit12 $cmreg,
          pModelMean_Int_fit13 $cmreg,type="html",align=TRUE,
          keep=c("TrtMnt_demend_y", "area_Trt", "Pop_Trt", "MeanT_Trt","MeanP_Trt","MaxT_Trt","MaxP_Trt","MinT_Trt","MinP_Trt","Year_Trt", "Slope_Trt", "Elevation_Trt", "RivDist_Trt", "RoadDist_Trt"),
          #covariate.labels=c("TrtMnt","MeanT","MeanP","Pop","MaxT","MaxP","MinT","MinP","Year"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Mean NDVI")
)
