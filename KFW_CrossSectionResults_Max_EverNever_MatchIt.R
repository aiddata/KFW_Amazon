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
           "post_trend_precip_mean","post_trend_precip_min","post_trend_precip_max",
           "MaxL_2010","MeanL_2010")


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

psmModel_R <-  "TrtBin ~ terrai_are + Pop_1990 + MeanT_1995 + pre_trend_temp_mean + pre_trend_temp_min + 
pre_trend_temp_max + MeanP_1995 + pre_trend_precip_min + 
pre_trend_NDVI_mean + pre_trend_NDVI_max + Slope + Elevation + MaxL_1995 + Riv_Dist + Road_dist +
pre_trend_precip_mean + pre_trend_precip_max"
#MeanL_1995

psmRes <- SCI::SpatialCausalPSM(dta_Shp,mtd="logit",psmModel_R,drop="support",visual=FALSE)

#-------------------------------------------------
#-------------------------------------------------
#Based on the Propensity Score Matches, pair comprable treatment and control units.
#-------------------------------------------------
#-------------------------------------------------
drop_set<- c(drop_unmatched=TRUE,drop_method="None",drop_thresh=0.5)
psm_Pairs <- SAT(dta = psmRes$data, mtd = "fastNN",constraints=c(groups="UF"),psm_eq = psmModel_R, ids = "id", drop_opts = drop_set, visual="TRUE", TrtBinColName="TrtBin")
#c(groups=c("UF"),distance=NULL)
trttable <- table (psm_Pairs@data$TrtBin)
View(trttable)


##create standardized dataset to produce standardized coefficients in models that are easy to output
#or to create normalized difference of means statistics for summary statistics

stvars <- c("terrai_are","Pop_1990", "MeanT_1995", "pre_trend_temp_mean",
           "pre_trend_temp_min", "pre_trend_temp_max", "MeanP_1995", "pre_trend_precip_min",
           "pre_trend_NDVI_mean", "pre_trend_NDVI_max","Slope","Elevation","MaxL_1995","MeanL_1995","Riv_Dist","Road_dist",
           "pre_trend_precip_mean", "pre_trend_precip_max",
           "NDVILevelChange_95_10","post_trend_temp_mean","post_trend_temp_min","post_trend_temp_max",
           "post_trend_precip_mean","post_trend_precip_min","post_trend_precip_max")

#standardize full unmatched dataset, then subset from that
#(rather than standardizing each dataset after it has been subsetted, which won't allow comparisons across all datasets)

#standardize unmatched dataset
dta_Shp_st<-dta_Shp@data
dta_Shp_st[stvars]<-lapply(dta_Shp_st[stvars],scale)

#subset to standardized matched without replacement dataset
pairs_id<-psm_Pairs@data$reu_id
psm_Pairs_st<-dta_Shp_st[dta_Shp_st$reu_id %in% pairs_id,]

#subset to standardized matched with replacement dataset
model_id<-model_data$reu_id
model_data_st<-dta_Shp_st[dta_Shp_st$reu_id %in% model_id,]

#test
summary(dta_Shp$Elevation)
summary(dta_Shp_st$Elevation)
summary(psm_Pairs$Elevation)
summary(psm_Pairs_st$Elevation)
summary(model_data$Elevation)
summary(model_data_st$Elevation)
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
          title="Regression Results", type="html", 
          star.cutoffs = c(0.05, 0.01, 0.001),
          omit.stat=c("f","ser"), align=TRUE)


#------------------------------
#Create summary stats table by treatment and pairs
#Include normalized difference in means
#------------------------------

stat_vars <- c("TrtBin","terrai_are","Pop_1990", "MeanT_1995", "MeanP_1995", 
               "Slope","Elevation","MaxL_1995","MeanL_1995","Riv_Dist","Road_dist")

#--
## *Never demarcated (all) vs. ever demarcated (all)*

#subset dataset with all communities for selected variables
dta_sub<-dta_Shp@data[stat_vars]
#create dataset with means, transpose, turn into dataframe
dta_Shp_stats <-aggregate(dta_sub, by=list(dta_sub$TrtBin), 
                    FUN=mean, na.rm=TRUE)
dta_Shp_stats1<-t(dta_Shp_stats)
dta_Shp_stats2<-data.frame(dta_Shp_stats1)
# Create a variable name column from rownames
dta_Shp_stats2$varname = rownames(dta_Shp_stats2)
# Reset the rownames of the original data
rownames(dta_Shp_stats2) = NULL
#rename columns to reflect demarcation status
names(dta_Shp_stats2)[names(dta_Shp_stats2)=="X1"] <- "PPTAL_nondem_all"
names(dta_Shp_stats2)[names(dta_Shp_stats2)=="X2"] <- "PPTAL_dem_all"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
dta_Shp_stats3<-dta_Shp_stats2[!(dta_Shp_stats2$varname %in% rows),]
#reorder rows alphabetically and assign variable id
dta_Shp_stats4<-dta_Shp_stats3[order(dta_Shp_stats3$varname),]
dta_Shp_stats4$var_id<-seq.int(nrow(dta_Shp_stats4))
dta_Shp_stats<-dta_Shp_stats4

##add normalized means and differences

#subset standardized dataset with all communities for selected variables
dta_sub_st<-dta_Shp_st[stat_vars]
#create dataset with means, transpose, turn into dataframe
dta_Shp_stats_st <-aggregate(dta_sub_st, by=list(dta_sub_st$TrtBin), 
                          FUN=mean, na.rm=TRUE)
dta_Shp_stats1_st<-t(dta_Shp_stats_st)
dta_Shp_stats2_st<-data.frame(dta_Shp_stats1_st)
# Create a variable name column from rownames
dta_Shp_stats2_st$varname = rownames(dta_Shp_stats2_st)
# Reset the rownames of the original data
rownames(dta_Shp_stats2_st) = NULL
#rename columns to reflect demarcation status
names(dta_Shp_stats2_st)[names(dta_Shp_stats2_st)=="X1"] <- "PPTAL_nondem_all_st"
names(dta_Shp_stats2_st)[names(dta_Shp_stats2_st)=="X2"] <- "PPTAL_dem_all_st"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
dta_Shp_stats3_st<-dta_Shp_stats2_st[!(dta_Shp_stats2_st$varname %in% rows),]
#create normalized difference var
dta_Shp_stats_st<-dta_Shp_stats3_st
dta_Shp_stats_st$ndiff_all<-abs(dta_Shp_stats_st$PPTAL_nondem_all_st-dta_Shp_stats_st$PPTAL_dem_all_st)
summ_stats<-merge(dta_Shp_stats,dta_Shp_stats_st)
#--

#--
## *All communities*

#create dataset with means, transpose, turn into dataframe
#create fake group variable to use aggregate function
dta_sub$group<-1
allcomms_stats <-aggregate(dta_sub,by= list(dta_sub$group),
                          FUN=mean, na.rm=TRUE)
allcomms_stats1<-t(allcomms_stats)
allcomms_stats2<-data.frame(allcomms_stats1)
# Create a variable name column from rownames
allcomms_stats2$varname = rownames(allcomms_stats2)
# Reset the rownames of the original data 
rownames(allcomms_stats2) = NULL
#drop unneeded rows
rows1<-c("TrtBin","Group.1","group")
allcomms_stats3<-allcomms_stats2[!(allcomms_stats2$varname %in% rows1),]

#merge all community stats with stats for ever vs. never demarcated
allcomms_stats<-allcomms_stats3
summ_stats1<-merge(summ_stats,allcomms_stats)
#--

#--
## *Matched without Replacement Stats*

#subset data for matched without replacement
pairs_sub<-psm_Pairs@data[stat_vars]
#create dataset with means, transpose, turn into dataframe
pairs_stats <-aggregate(pairs_sub, by=list(pairs_sub$TrtBin), 
                          FUN=mean, na.rm=TRUE)
pairs_stats1<-t(pairs_stats)
pairs_stats2<-data.frame(pairs_stats1)
# Create a variable name column from rownames
pairs_stats2$varname = rownames(pairs_stats2)
# Reset the rownames of the original data
rownames(pairs_stats2) = NULL
#rename columns to reflect demarcation status
names(pairs_stats2)[names(pairs_stats2)=="X1"] <- "nondem_worep"
names(pairs_stats2)[names(pairs_stats2)=="X2"] <- "dem_worep"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
pairs_stats3<-pairs_stats2[!(pairs_stats2$varname %in% rows),]
#merge
pairs_stats<-pairs_stats3
summ_stats2<-merge(summ_stats1,pairs_stats)

## Add normalized differences
#subset standardized dataset with pair communities for selected variables
pairs_sub_st<-psm_Pairs_st[stat_vars]
#create dataset with means, transpose, turn into dataframe
pairs_stats_st <-aggregate(pairs_sub_st, by=list(pairs_sub_st$TrtBin), 
                             FUN=mean, na.rm=TRUE)
pairs_stats1_st<-t(pairs_stats_st)
pairs_stats2_st<-data.frame(pairs_stats1_st)
# Create a variable name column from rownames
pairs_stats2_st$varname = rownames(pairs_stats2_st)
# Reset the rownames of the original data
rownames(pairs_stats2_st) = NULL
#rename columns to reflect demarcation status
names(pairs_stats2_st)[names(pairs_stats2_st)=="X1"] <- "nondem_worep_st"
names(pairs_stats2_st)[names(pairs_stats2_st)=="X2"] <- "dem_worep_st"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
pairs_stats3_st<-pairs_stats2_st[!(pairs_stats2_st$varname %in% rows),]
#create normalized difference var
pairs_stats_st<-pairs_stats3_st
pairs_stats_st$ndiff_worep<-abs(pairs_stats_st$nondem_worep_st-pairs_stats_st$dem_worep_st)
summ_stats3<-merge(summ_stats2,pairs_stats_st)
summ_stats<-summ_stats3

## *Matched with Replacement Stats*
#subset data for matched with replacement 
match_sub<-model_data[stat_vars]
#create dataset with means, transpose, turn into dataframe
match_stats <-aggregate(match_sub, by=list(match_sub$TrtBin), 
                        FUN=mean, na.rm=TRUE)
match_stats1<-t(match_stats)
match_stats2<-data.frame(match_stats1)
# Create a variable name column from rownames
match_stats2$varname = rownames(match_stats2)
# Reset the rownames of the original data
rownames(match_stats2) = NULL
#rename columns to reflect demarcation status
names(match_stats2)[names(match_stats2)=="X1"] <- "nondem_rep"
names(match_stats2)[names(match_stats2)=="X2"] <- "dem_rep"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
match_stats3<-match_stats2[!(match_stats2$varname %in% rows),]
#merge
match_stats<-match_stats3
summ_stats3<-merge(summ_stats,match_stats)
summ_stats<-summ_stats3

## Add normalized differences
#subset standardized dataset with pair communities for selected variables
match_sub_st<-model_data_st[stat_vars]
#create dataset with means, transpose, turn into dataframe
match_stats_st <-aggregate(match_sub_st, by=list(match_sub_st$TrtBin), 
                           FUN=mean, na.rm=TRUE)
match_stats1_st<-t(match_stats_st)
match_stats2_st<-data.frame(match_stats1_st)
# Create a variable name column from rownames
match_stats2_st$varname = rownames(match_stats2_st)
# Reset the rownames of the original data
rownames(match_stats2_st) = NULL
#rename columns to reflect demarcation status
names(match_stats2_st)[names(match_stats2_st)=="X1"] <- "nondem_rep_st"
names(match_stats2_st)[names(match_stats2_st)=="X2"] <- "dem_rep_st"
#drop unneeded rows
rows<-c("TrtBin","Group.1")
match_stats3_st<-match_stats2_st[!(match_stats2_st$varname %in% rows),]
#create normalized difference var
match_stats_st<-match_stats3_st
match_stats_st$ndiff_rep<-abs(match_stats_st$nondem_rep_st-match_stats_st$dem_rep_st)
summ_stats4<-merge(summ_stats,match_stats_st)
summ_stats<-summ_stats4

#Formatting for summary stats table
stats_table<-summ_stats
names(stats_table)[names(stats_table)=="allcomms_stats1"]<-"allPPTAL"
stats_table<-stats_table[,-grep("(_st)",names(stats_table))]
stats_table$varname_title<-c("Elevation (m)","NDVI Max, 1995","NDVI Mean, 1995","Mean Precipitation, 1995",
                             "Mean Temperature, 1995", "Population Density, 1990","River Distance (m)",
                             "Road Distance (m)","Slope (degrees)","Area (hectares)")
stats_table<-stats_table[c(1,13,6,2,3,5,7:12)]
stats_table<-stats_table[order(stats_table$varname_title),]

#output table using stargazer, then import into excel and format column names and add number of observations
stargazer(stats_table, type = "html", summary = FALSE, rownames = FALSE)

#to check values
#Unmatched
describeBy(dta_Shp_st$terrai_are, dta_Shp_st$TrtBin)
#Matched w/out replacement
describeBy(psm_Pairs_st$terrai_are, psm_Pairs_st$TrtBin)
#Matched w/replacement
describeBy(model_data$terrai_are, model_data$TrtBin)

## -------
# Identify 2010 Max NDVI values for Meta Analysis Paper, done in August 2018
## -------

#UNMATCHED
#this should match Table 5 from final JEEM manuscipt
describeBy(dta_Shp$MaxL_1995, dta_Shp_st$TrtBin)
#this is the new value for all communities, and for never vs. ever
mean(dta_Shp$MaxL_2010)
describeBy(dta_Shp$MaxL_2010, dta_Shp_st$TrtBin)

# MATCHED W/OUT REPLACEMENT
#this should match Table 6 from final JEEM manuscript
describeBy(psm_Pairs$MaxL_1995, psm_Pairs$TrtBin)
#this is the new value for never vs. ever
psm_Pairs_never<-psm_Pairs[psm_Pairs$TrtBin==0,]
psm_Pairs_ever<-psm_Pairs[psm_Pairs$TrtBin==1,]
mean(psm_Pairs_never$MaxL_2010)
mean(psm_Pairs_ever$MaxL_2010)

# MATCHED W/REPLACEMENT
#this should match Table 6 from final JEEM manuscript
describeBy(model_data$MaxL_1995, model_data$TrtBin)
#this is the new 2010 value for never vs. ever
model_data_never<-model_data[model_data$TrtBin==0,]
model_data_ever<-model_data[model_data$TrtBin==1,]
mean(model_data_never$MaxL_2010)
mean(model_data_ever$MaxL_2010)

