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
## Create panel dataset to create yearly NDVI for each community
#-------------------------------------------------
#-------------------------------------------------

varList = c("MaxL_")
psm_Long <- BuildTimeSeries(dta=dta_Shp,idField="reu_id",varList_pre=varList,1982,2010,colYears=c("demend_y"),interpYears=c("TrtBin"))
psm_Long$Year <- as.numeric(psm_Long$Year)





#-------------------------------------------------
#-------------------------------------------------
## View Time Series and GGPlot to create NDVI graphs
#-------------------------------------------------
#-------------------------------------------------

#Original Figure 5 (Dem vs. Non-Dem, with means as well as all comms showing)
ViewTimeSeries(dta=dta_Shp,IDfield="reu_id",TrtField="TrtBin",idPre="MaxL_[0-9][0-9][0-9][0-9]")

#Figure 5 for JEEM 2nd Resubmission (mean only showing)
ggplot(data = psm_Long, aes(x=Year, y=MaxL_, group=reu_id,linetype=factor(TrtBin))) + 
  geom_line(size=0, linetype=2) +
  stat_summary(fun.y=mean,aes(x=Year,y=MaxL_,group=TrtBin),data=psm_Long,geom='line',size=1.5)+
  stat_summary(fun.data="mean", conf.int=0.95, geom="errorbar", 
               width=0.1, linetype=2, colour="blue")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="Year",y="Max NDVI",linetype="Demarcation")

#scratch
my_ci <- function(x) data.frame(
  y=mean(x), 
  ymin=mean(x) - 1.96 * sd(x), 
  ymax=mean(x) + 1.96 * sd(x)
)
my_mean<-function(x) data.frame(
  y=mean(x)
)

ggplot(psm_Long, aes(Year, MaxL_)) + 
  #geom_point(color="red") +
  stat_summary(fun.data="my_ci", geom="smooth")

ggplot(psm_Long, aes(x=Year, y=MaxL_,linetype=factor(TrtBin))) +
  stat_summary(data=dem,fun.data="my_mean",geom="line")+
  stat_summary(data=nondem,fun.data="my_ci", geom="smooth")
nondem<-psm_Long[psm_Long$TrtBin==0,]
dem<-psm_Long[psm_Long$TrtBin==1,]

ggplot(data = psm_Long, aes(x=Year, y=MaxL_, group=reu_id,linetype=factor(TrtBin))) + 
  geom_line(size=0, linetype=2) +
  stat_summary(fun.y=mean,aes(x=Year,y=MaxL_,group=TrtBin),data=psm_Long,geom='line',size=1.5)+
  stat_summary(fun.data="my_ci", geom="smooth")+
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="Year",y="Max NDVI",linetype="Demarcation")



# time series graph with year effects stripped out

reg=lm(MaxL_ ~ factor(Year), data=psm_Long)
psm_Long$resid <- residuals(reg)
plot(psm_Long$resid)

psm_Long2 <- psm_Long
psm_Long2 <- psm_Long[c("Year","resid","reu_id","TrtBin")]

#Figure 6 for JEEM 2nd resubmission (mean only showing)
ggplot(data = psm_Long2, aes(x=Year, y=resid,linetype=factor(TrtBin))) + 
  geom_line(size=0, linetype=2) +
  stat_summary(fun.y=mean,aes(x=Year, y=resid,group=TrtBin),data=psm_Long2,geom='line',size=1.5) +
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  labs(x="Year",y="Residual",linetype="Demarcation")


## Scratch

#Original Figure 6 (Dem vs. non dem, with year effects stripped out, means + all comm values showing)
ggplot(data = psm_Long2, aes(x=Year, y=resid, group=reu_id,colour=factor(TrtBin))) + 
  geom_line(size=.5, linetype=2) +
  stat_summary(fun.y=mean,aes(x=Year, y=resid, group=TrtBin,colour=factor(TrtBin)),data=psm_Long2,geom='line',size=1.5)+
  theme(axis.text.x=element_text(angle=90,hjust=1))


