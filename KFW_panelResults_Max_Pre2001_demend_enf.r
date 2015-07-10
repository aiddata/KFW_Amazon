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
drop_set<- c(drop_unmatched=TRUE,drop_method="NONE",drop_thresh=0.25)
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
psm_Long <- BuildTimeSeries(dta=psm_Pairs,idField="reu_id",varList_pre=varList,1982,2010,colYears=c("demend_y","enforce_st"),interpYears=c("Slope","Road_dist","Riv_Dist","UF","Elevation","terrai_are","Pop_","MeanT_","MeanP_","MaxT_","MaxP_","MinP_","MinT_"))
psm_Long$Year <- as.numeric(psm_Long$Year)

#Create post-2004 dummy
psm_Long$Post2004 <- 0
psm_Long$Post2004[psm_Long$Year >= 2004] <- 1

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


##texreg to output regression results visualizations
texreg::plotreg(pModelMax_B_fit$cmreg, custom.model.names=c("Panel Results, Max NDVI"), 
                omit.coef="(Pop_)|(Min)|(Mean)|(Max)|(Year)|(match)|(Intercept)|(factor)", 
                custom.note="standard deviation")

texreg::plotreg(pModelMax_C_fit$cmreg,custom.model.names=c("Panel Results, Max NDVI"), 
                omit.coef="(Pop_)|(Min)|(Mean)|(Max)|(Year)|(Intercept)|(factor)", 
                custom.note="standard deviation")

## stargazer output with variable labels
stargazer(pModelMax_A_fit$cmreg,pModelMax_B_fit$cmreg,pModelMax_C_fit$cmreg,pModelMax_D_fit$cmreg,pModelMax_E_fit$cmreg,type="html",align=TRUE,keep=c("TrtMnt_demend_y","TrtMnt_enforce_st","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year","Post2004","TrtMnt_demend_y:Post2004","Post2004:Road_dist","TrtMnt_demend_y:Road_dist","TrtMnt_demend_y:Post2004:Road_dist"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Enforcement Support)","Mean Temperature","Mean Precipitation","Population","Max Temperature","Max Precipitation","Min Temperature","Min Precipitation","Year","Post2004","Post2004*TreatmentDemarcation","Post2004*Road Distance","TreatmentDemarcation*Road Distance","Post2004*TreatmentDemarcation*Road Distance"),
          omit.stat=c("f","ser"),
          title="Regression Results",
          dep.var.labels=c("Max NDVI")
)

## stargazer output with variable names instead of labels
# stargazer(pModelMax_A_fit$cmreg,pModelMax_B_fit$cmreg,pModelMax_C_fit$cmreg,pModelMax_D_fit$cmreg,pModelMax_E_fit$cmreg,type="html",align=TRUE,keep=c("TrtMnt_demend_y","TrtMnt_enforce_st","MeanT_","MeanP_","Pop_","MaxT_","MaxP_","MinT_","MinP_","Year","Post2004","TrtMnt_demend_y:Post2004","Post2004:Road_dist","TrtMnt_demend_y:Road_dist","TrtMnt_demend_y:Post2004:Road_dist"),
#           covariate.labels=c("TrtMntDem","TrtMntEnf","MeanT","MeanP","Pop","MaxT","MaxP","MinT","MinP","Year","Post2004","Post04*TrtMntDem","Post04*RoadDist","TrtMntDem*RoadDist","TrtMntDem*RoadDist*Post2004"),
#           omit.stat=c("f","ser"),
#           title="Regression Results",
#           dep.var.labels=c("Max NDVI")
# )