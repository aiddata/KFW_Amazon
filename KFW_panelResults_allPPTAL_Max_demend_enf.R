#-------------------------------------------------
#-------------------------------------------------
#Panel Models - KFW Community-Level
#Demarcated and Non-Demarcated Communities
#Max NDVI NDVI (LTDR)
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

#Drop non-PPTAL communities (but leave in those that were never demarcated)
dta_Shp$NA_check <- 0
dta_Shp$NA_check[is.na(dta_Shp$reu_id)] <- 1
dta_Shp <- dta_Shp[dta_Shp$NA_check != 1,]

#no sample trimming
#will create demarcation treatment variable after panel dataset is created

#-------------------------------------------------
#-------------------------------------------------
#Convert from a wide-form dataset for the Cross-sectional 
#to a long-form dataset for the panel model.
#-------------------------------------------------
#-------------------------------------------------

#Drop data for unused years (after 2010) in order to allow "reshape" to work
kfw_reshape <- dta_Shp
kfw_reshape1<-kfw_reshape[,-grep("(2011)",names(kfw_reshape))]
kfw_reshape2<-kfw_reshape1[,-grep("(2012)",names(kfw_reshape1))]
kfw_reshape3<-kfw_reshape2[,-grep("(2013)",names(kfw_reshape2))]
kfw_reshape4<-kfw_reshape3[,-grep("(2014)",names(kfw_reshape3))]
kfw_reshape5<-kfw_reshape4[,-grep("(2015)",names(kfw_reshape4))]
kfw_reshape6<-kfw_reshape5[,-grep("(2020)",names(kfw_reshape5))]

kfw_wide<-kfw_reshape6

#Order variables chronologically to allow reshape to work properly
kfw_wide<-kfw_wide[,order(names(kfw_wide))]

#Identify variables where values will change yearly in panel dataset
MeanT_<-grep("MeanT_",names(kfw_wide))
MeanP_<-grep("MeanP_",names(kfw_wide))
MinT_<-grep("MinT_",names(kfw_wide))
MaxT_<-grep("MaxT_",names(kfw_wide))
MinP_<-grep("MinP_",names(kfw_wide))
MaxP_<-grep("MaxP_",names(kfw_wide))
MaxL_<-grep("MaxL_",names(kfw_wide))
MeanL_<-grep("MeanL_",names(kfw_wide))
MedL_<-grep("MedL_",names(kfw_wide))
Pop_<-grep("Pop_",names(kfw_wide))

#Convert to long form panel dataset
all_reshape <- c(MeanT_,MeanP_,MinT_,MaxT_,MinP_,MaxP_,MaxL_,MeanL_,MedL_,Pop_)
psm_Long <- reshape(kfw_wide@data, varying=all_reshape, direction="long",idvar="id",sep="_",timevar="Year")

#--------------------------
# Construct Treatment Vars
#--------------------------
#Create years to demarcation
psm_Long$yrtodem <- NA
psm_Long$yrtodem=psm_Long$Year - psm_Long$demend_y

#Create demarcation treatment variable, using demend_y
#0 in years prior to demarcation, turns to 1 in year of demarcation
psmtest3 <- psm_Long
psmtest3$trtdem <- 0
psmtest3$trtdem[which(psmtest3$Year<psmtest3$demend_y)]<-0
psmtest3$trtdem[which(psmtest3$Year>=psmtest3$demend_y)]<-1

psm_Long <- psmtest3

#Correct communities with enforcement start date prior to demarcation
#Only 1 community, reu_id==84
psmtest5 <- psm_Long
psmtest5$enfdiff= psmtest5$enforce_st - psmtest5$demend_y
summary(psmtest5$enfdiff)
psmenf <- subset(psmtest5, psmtest5$enfdiff<0)
table(psmenf$reu_id)

psmtest5$enforce_st[which(psmtest5$reu_id==84)]<-2007

#Create enforcement treatment var
#0 before enforcement support begins, turns to 1 in the year that enforcement support starts
psmtest5$trtenf <- 0
psmtest5$trtenf[which(psmtest5$Year>=psmtest5$enforce_st)]<-1
table(psmtest5$trtenf,psmtest5$trtdem)
#Should only apply to communities that have been demarcated
# For non-demarcated communities, set all values back to 0
psmtest5$trtenf[is.na(psmtest5$demend_y)]<-0
table(psmtest5$trtenf,psmtest5$trtdem)

psm_Long <- psmtest5

#Create post-2004 dummy
psm_Long$Post2004 <- 0
psm_Long$Post2004[psm_Long$Year >= 2004] <- 1

#----------------------
# Regression Models
#----------------------


pModelMax_A <- "MaxL ~ trtdem + trtenf + factor(reu_id)"
pModelMax_B <- "MaxL ~ trtdem + trtenf + Pop + MeanT + MeanP + MaxT + MaxP + MinT + MinP  + factor(reu_id) "
pModelMax_C <- "MaxL ~ trtdem + trtenf + Pop + MeanT + MeanP + MaxT + MaxP + MinT + MinP  + Year + factor(reu_id)"
pModelMax_C2 <- "MaxL ~ trtdem + trtenf + Pop + MeanT + MeanP+ MaxT + MaxP + MinT + MinP  + factor(Year) + factor(reu_id)"

pModelMax_D <- "MaxL ~ trtdem + Pop+ MeanT + MeanP+ MaxT + MaxP + MinT + MinP + factor(Year) + factor(reu_id) + Post2004 + Post2004*trtdem + Post2004*trtdem*Road_dist + Post2004*Road_dist"
pModelMax_E <- "MaxL ~ trtdem + trtenf+ Pop + MeanT + MeanP + MaxT + MaxP + MinT + MinP + factor(Year) + factor(reu_id) + Post2004 + Post2004*trtdem + Post2004*trtdem*Road_dist + Post2004*Road_dist"


pModelMax_A_fit <- Stage2PSM(pModelMax_A ,psm_Long,type="cmreg", table_out=TRUE,opts=c("reu_id","Year"))
pModelMax_B_fit <- Stage2PSM(pModelMax_B ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C_fit <- Stage2PSM(pModelMax_C ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_C2_fit <- Stage2PSM(pModelMax_C2 ,psm_Long,type="cmreg", table_out=TRUE,opts=c("reu_id","Year"))

pModelMax_D_fit <- Stage2PSM(pModelMax_D ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))
pModelMax_E_fit <- Stage2PSM(pModelMax_E ,psm_Long,type="cmreg", table_out=TRUE, opts=c("reu_id","Year"))


## stargazer output with variable labels

#Used for JEEM 2nd resubmission
#imported html into excel and did some reformatting within
stargazer(pModelMax_A_fit$cmreg,pModelMax_B_fit$cmreg,pModelMax_C_fit$cmreg,
          pModelMax_C2_fit$cmreg,pModelMax_D_fit$cmreg,pModelMax_E_fit$cmreg,
          type="html", align=TRUE,
          omit=c("factor","Constant"),
          covariate.labels=c("Treatment (Demarcation)","Treatment (Demarcation + Enforcement Support)","Population","Mean Temp",
                             "Mean Precip","Max Temp","Max Precip","Min Temp","Min Precip","Year",
                             "Treatment(Dem)*Post 2004","Post 2004*Road Dist","Treatment(Dem)*Road Dist",
                             "Treatment(Dem)*Post 2004*Road Dist"),
          omit.stat=c("f","ser"),
          add.lines=list(c("Observations","4379","4379","4379","4379","4379","4379"),
                         c("Community Fixed Effects?","Yes","Yes","Yes","Yes","Yes","Yes"),
                         c("Year Fixed Effects?","No","No","No","Yes","Yes","Yes")),
          title="Regression Results",
          dep.var.labels=c("Max NDVI"))


#-----------------
#Scratch/Archived
#-----------------

#Create arc of deforestation dummy
psm_Long$arc<-NA
psm_Long$arc[which(psm_Long$UF=="PA")] <-1
psm_Long$arc[which(psm_Long$UF!="PA")]<-0

#Create dummy for each year post-demarcation
DemYears <- summaryBy(TrtMnt_demend_y~reu_id, data=psm_Long,FUN=c(mean,sum))
psm_Long2 <- psm_Long
psm_Long3 <- merge(psm_Long2, DemYears, by="reu_id")
psm_Long<-psm_Long3


pModelMax_A <- lm(MaxL_ ~ trtdem + trtenf , data=psm_Long)
clusterA<-cluster.vcov(pModelMax_A,cbind(psm_Long$Year,psm_Long$id), force_posdef=TRUE)
pModelMax_A_fit <- coeftest(pModelMax_A,clusterA)


Model4<- lm(lfreq ~ trt_dem + ifreq+
              Pop + 
              ntl +
              factor(year) + factor(id),
            data=panel_data)
cluster4 <- cluster.vcov(Model4, cbind(panel_data$year, panel_data$id), force_posdef=TRUE)
pModelMax_A_fit <- coeftest(Model4, cluster4)

