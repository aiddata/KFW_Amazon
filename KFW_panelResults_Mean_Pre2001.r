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
#-------------------------------------------------
#-------------------------------------------------
varList = c("MeanL_","MaxL_")
psm_Long <- BuildTimeSeries(dta=psm_Pairs,idField="reu_id",varList_pre=varList,1982,2010,colYears="demend_y",interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","TrtMnt","MaxT_","MaxP_","MinP_","MinT_"))
psm_Long$Year <- as.numeric(psm_Long$Year)


##INTERACTIONS

#Create interaction term with post2001 & treatment
psm_Long["post2001_Bin"] <- 0
psm_Long$post2001_Bin[psm_Long$Year > 2001] <- 1
psm_Long["post2001_Trt"] <- ((psm_Long$post2001_Bin)*(psm_Long$TrtMnt))

#Interaction terrain area & treatment
psm_Long["area_Trt"] <- ((psm_Long$terrai_are)*(psm_Long$TrtMnt))

#Interaction Pop_ & treatment
psm_Long["Pop_Trt"] <-((psm_Long$Pop_)*(psm_Long$TrtMnt))

#Interaction MeanT & treatment
psm_Long["MeanT_Trt"] <-((psm_Long$MeanT_)*(psm_Long$TrtMnt))

#Interaction MeanP & treatment
psm_Long["MeanP_Trt"] <- ((psm_Long$MeanP_)*(psm_Long$TrtMnt))

#Interaction MaxT & treatment
psm_Long["MaxT_Trt"] <- ((psm_Long$MaxT_)*(psm_Long$TrtMnt))

#Interaction MaxP & treatment
psm_Long["MaxP_Trt"] <- ((psm_Long$MaxP_)*(psm_Long$TrtMnt))

#Interaction MinT & Treatment
psm_Long["MinT_Trt"] <- ((psm_Long$MinT_)*(psm_Long$TrtMnt))

#Interaction MinP & Treatment
psm_Long["MinP_Trt"] <- ((psm_Long$MinP_)*(psm_Long$TrtMnt))

#Interaction Year & Treatment
psm_Long["Year_Trt"] <- ((psm_Long$Year)*(psm_Long$TrtMnt))

#Interaction Slope & treatment
psm_Long["Slope_Trt"] <- ((psm_Long$Slope)*(psm_Long$TrtMnt))

#Interaction Elevation & treatment
psm_Long["Elevation_Trt"] <- ((psm_Long$Elevation)*(psm_Long$TrtMnt))

#Interaction River Distance & treatment
psm_Long["RivDist_Trt"] <- ((psm_Long$Riv_Dist)*(psm_Long$TrtMnt))

#Interaction Road Distance & treatment
psm_Long["RoadDist_Trt"] <- ((psm_Long$Road_dist)*(psm_Long$TrtMnt))

pModelMean_A <- "MeanL_ ~ TrtMnt + factor(reu_id) "
pModelMean_B <- "MeanL_ ~ TrtMnt + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) "
pModelMean_C <- "MeanL_ ~ TrtMnt + MeanT_ + MeanP_ + Pop_ + MaxT_ + MaxP_ + MinT_ + MinP_  + factor(reu_id) + Year"

pModelMean_A_fit <- Stage2PSM(pModelMean_A ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_B_fit <- Stage2PSM(pModelMean_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMean_C_fit <- Stage2PSM(pModelMean_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))

#Model with interactions
pModelMean_Int <- "MeanL_ ~ TrtMnt + MeanT_ + MeanP_ + Pop_ + MaxT_ + \n
MaxP_ + MinT_ + MinP_ + factor(reu_id) + Year + area_Trt+ Pop_Trt + \n
MeanT_Trt + MeanP_Trt + MaxT_Trt + MaxP_Trt + MinT_Trt + MinP_Trt + \n
Year_Trt + Slope_Trt + Elevation_Trt + RivDist_Trt + RoadDist_Trt"

#RoadDist was only significant interaction
pModelMean_Int4 <- "MeanL_ ~ TrtMnt + MeanT_ + MeanP_ + Pop_ + MaxT_ + \n
MaxP_ + MinT_ + MinP_ + factor(reu_id) + Year + RoadDist_Trt"

pModelMean_Int4_fit <- Stage2PSM(pModelMean_Int4 ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


stargazer(pModelMean_A_fit $cmreg,pModelMean_B_fit $cmreg,pModelMean_C_fit $cmreg,pModelMean_Int_fit $cmreg,type="html",align=TRUE,
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Mean NDVI")
)

