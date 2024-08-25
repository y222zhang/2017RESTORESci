rm(SS_output)
rm(SS_plots)
rm(list=ls())

load_libraries<-function() {
#  library(PBSmodelling)
  library(snowfall)
  library(parallel)
  library(snow)
  library(foreach)
  library(doSNOW)
  library(PBSadmb)
  library(plyr)
  library(dplyr)
  library(tidyr)
  library(r4ss)
  library(matrixStats)
}
load_libraries()

#library("devtools")
#devtools::install_github("r4ss/r4ss") #, ref="v1.23.1")
options(digits=16)

setwd("~/Dropbox/ (Whitespace Conflict)/FIU/RESEARCH/red snapper MSE/Dan 2018 red snapper assessment files/OFL")

#Part 1 record key parameters as soon as the zip file was read in
## directories
direct <- getwd()

##basic readin
base <- SS_output(dir = direct,printstats = T, covar=T, cormax=0.70, forecast=F,printhighcor=50, printlowcor=50)

##start year and end year
yr_start<-base$startyr
#endyear was used in Step 2
yr_end<-base$endyr

##age readin and setting initial population
#The following 5 rows were used in Step 2 and 3
#stock_1: East, stock_2: West
age_1<-base$agebins
age_2<-age_1
stock_1_mean<-c(t(base$natage[(base$natage$Area==1)&(base$natage$Yr==yr_end)&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))
stock_2_mean<-c(t(base$natage[(base$natage$Area==2)&(base$natage$Yr==yr_end)&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))
#sum(stock_1_mean) #matches the SA report 4.1 Spawning Output unit 1000s
#sum(stock_2_mean) #matches the SA report 4.1 Spawning Output unit 1000s

#parepare for MSE init population
iniPopu<-cbind(age_1,stock_1_mean,stock_2_mean)
colnames(iniPopu) <- c("age", "stock_1_mean", "stock_2_mean")

##biological parameters weight unit kg
# Number of eggs for each individual

fec_at_age<-c(t(base$ageselex[(base$ageselex$Factor=="Fecund")&(base$ageselex$Yr==yr_end),as.character(base$agebins)]))
fec_at_age_1<-fec_at_age
fec_at_age_2<-fec_at_age

spawning_output_1<-stock_1_mean*fec_at_age_1
spawning_output_2<-stock_2_mean*fec_at_age_2
#sum(spawning_output_1) #matches the SA report 4.1 Spawning Output 
#sum(spawning_output_2) #matches the SA report 4.1 Spawning Output 

# Biomass for each individual, unit mt # use weight at age at the beginning of the year
weight_at_age<-c(t(base$mean_body_wt[(base$mean_body_wt$Yr==yr_start),as.character(base$agebins)]))
weight_at_age_1<-weight_at_age
weight_at_age_2<-weight_at_age

#B_at_age_1<-stock_1_mean*weight_at_age_1 #biomass unit 1000*kg=mt
#B_at_age_2<-stock_2_mean*weight_at_age_2 #biomass unit 1000*kg=mt
#sum(B_at_age_1) #matches the SA report 4.1 Spawning Output unit mt
#sum(B_at_age_2) #matches the SA report 4.1 Spawning Output unit mt

bioPara<-cbind(age_1, weight_at_age_1,fec_at_age_1,weight_at_age_2,fec_at_age_2)
colnames(bioPara) <- c("age", "weight_at_age_Area1", "fec_at_age_Area1", "weight_at_age_Area2", "fec_at_age_Area2")

#Natural mortality
#Age-specific natural mortality rates (M) for Gulf of Mexico red snapper assuming a
#Lorenzen mortality curve rescaled to an average M = 0.0943. The column labeled M represents
#the average natural mortality experienced from July 1-June 30 (i.e., a birth year). The label Adj.
#M indicates the values used in the SS3 model to account for SS advancing age on January 1.

M<-c(t(base$M_at_age[(base$M_at_age$Year==yr_end),as.character(base$agebins)]))
for(i.M in 2:length(M)) {
  if(is.na(M[i.M])){
    M[i.M]=M[i.M-1]
  }
}

#if(unique(base$M_at_age$Bio_Pattern)==1){
  M_1<-M
  M_2<-M
#}

natM<-cbind(age_1, M_1, M_2)
colnames(natM) <- c("age","M_Area1", "M_Area2")

#Fraction before spawning, currently 0.5
season_factor<-base$seasfracs #unit year

#recruitment
#recruitment, read parameter
dat6<- base$parameters
steepness<-base$parameters[base$parameters$Label=="SR_BH_steep","Value"]
R0_late<-exp(base$parameters[base$parameters$Label=="SR_LN(R0)","Value"]) #unit 1000 R0 after 1984
R_offset_para<-base$parameters[base$parameters$Label=="SR_envlink","Value"]
R0_early<-R0_late*exp(R_offset_para)

SSB0_1<-base$Dynamic_Bzero[(base$Dynamic_Bzero$Era=="VIRG"),"SSB_area1"]
SSB0_2<-base$Dynamic_Bzero[(base$Dynamic_Bzero$Era=="VIRG"),"SSB_area2"]
SSB0<-SSB0_1+SSB0_2

sigma_R<-base$parameters[base$parameters$Label=="SR_sigmaR","Value"] #standard deviation of logged recruitment

#unit 1000s
#spawning_output_rec_1<-c(t(base$natage[(base$natage$Area==1)&(base$natage$Yr==(yr_end))&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))*fec_at_age_1
#spawning_output_rec_2<-c(t(base$natage[(base$natage$Area==2)&(base$natage$Yr==(yr_end))&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))*fec_at_age_2
#spawning_output_rec<-spawning_output_rec_1+spawning_output_rec_2

#unit 1000s
Rhist_1<-base$natage[(base$natage$Area==1)&(base$natage$Yr<=yr_end)&(base$natage$Yr>=yr_start)&(base$natage$`Beg/Mid`=="B"),"0"]
Rhist_2<-base$natage[(base$natage$Area==2)&(base$natage$Yr<=yr_end)&(base$natage$Yr>=yr_start)&(base$natage$`Beg/Mid`=="B"),"0"]

Rhist_early<-Rhist_1 + Rhist_2
Rhist_late<-Rhist_early[(length(Rhist_early)-(yr_end-1984)):length(Rhist_early)]

#Other information from stock assessment files

#Age length key
length_age_key<-base$ALK[,,1]
length_age_key_stock1<-length_age_key
length_age_key_stock2<-length_age_key

#read F
#Fishing fleets definitions (14) Directed fleet landings and discards (8) Bycatch fleets (discards only) (6)
#Commercial Vertical line (HL) E/W 1872-2016
#Commercial Longline (LL) E/W 1980-2016
#Recreational Private/Charter (MRFSS/MRIP) E/W 1950-2016
#Recreational Headboat (HBT) E/W 1950-2016
#Commercial Closed Season or zero IFQ allocation (C_Closed) E/W 1991-2016
#Recreational Closed Season (R_Closed) E/W 1997-2016
#Shrimp Bycatch (SHR) E/W 1950/1946 (respectively)-2015

#F in the last three years, add up together not the same as stock assessment Table 4.2 due to the different definations.
hl_e_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_1']))
hl_w_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_2']))
ll_e_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_3']))
ll_w_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_4']))
mrip_e_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_5']))
mrip_w_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_6']))
hbt_e_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_7']))
hbt_w_pred_F  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_8']))
comm_closed_e_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_9']))
comm_closed_w_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_10']))
rec_closed_e_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_11']))
rec_closed_w_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_12']))
shrimp_e_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'F:_13']))
shrimp_w_pred_F <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'F:_14']))

hl_e_pred_F_ave <- mean(hl_e_pred_F)
hl_w_pred_F_ave <- mean(hl_w_pred_F)
ll_e_pred_F_ave <- mean(ll_e_pred_F) 
ll_w_pred_F_ave <- mean(ll_w_pred_F)
mrip_e_pred_F_ave <- mean(mrip_e_pred_F)
mrip_w_pred_F_ave <- mean(mrip_w_pred_F)
hbt_e_pred_F_ave <- mean(hbt_e_pred_F)
hbt_w_pred_F_ave <- mean(hbt_w_pred_F)
comm_closed_e_pred_F_ave <- mean(comm_closed_e_pred_F)
comm_closed_w_pred_F_ave <- mean(comm_closed_w_pred_F)
rec_closed_e_pred_F_ave <- mean(rec_closed_e_pred_F)
rec_closed_w_pred_F_ave <- mean(rec_closed_w_pred_F)
shrimp_e_pred_F_ave <- mean(shrimp_e_pred_F)
shrimp_w_pred_F_ave <- mean(shrimp_w_pred_F)

#Last year selectivity
hl_e_selex <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,1), as.character(base$agebins)]))
hl_w_selex   <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,2), as.character(base$agebins)]))
ll_e_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,3), as.character(base$agebins)]))
ll_w_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,4), as.character(base$agebins)]))
mrip_e_selex   <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,5), as.character(base$agebins)]))
mrip_w_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,6), as.character(base$agebins)]))
hbt_e_selex   <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,7), as.character(base$agebins)]))
hbt_w_selex   <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,8), as.character(base$agebins)]))
comm_closed_e_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,9), as.character(base$agebins)]))
comm_closed_w_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,10), as.character(base$agebins)]))
rec_closed_e_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,11), as.character(base$agebins)]))
rec_closed_w_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,12), as.character(base$agebins)]))
shrimp_e_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,13), as.character(base$agebins)]))
shrimp_w_selex  <- c(t(base$ageselex[(base$ageselex$Yr==yr_end) & is.element(base$ageselex$Factor,"Asel") & is.element(base$ageselex$Fleet,14), as.character(base$agebins)]))

#Last year retention rate for fleet 9-14, retention rate 0.
hl_e_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,1), 6:dim(base$sizeselex)[2]]
#hl_e_retention_multiplier<-hl_e_retention_len[length(hl_e_retention_len)]
hl_w_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,2), 6:dim(base$sizeselex)[2]]
#hl_w_retention_multiplier<-hl_w_retention_len[length(hl_w_retention_len)]

ll_e_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,3), 6:dim(base$sizeselex)[2]]
#ll_e_retention_multiplier<-ll_e_retention[length(ll_e_retention)]
ll_w_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,4), 6:dim(base$sizeselex)[2]]
#ll_w_retention_multiplier<-ll_w_retention[length(ll_w_retention)]

mrip_e_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,5), 6:dim(base$sizeselex)[2]]
#mrip_e_retention_multiplier<-mrip_e_retention[length(mrip_e_retention)]
mrip_w_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,6), 6:dim(base$sizeselex)[2]]
#mrip_w_retention_multiplier<-mrip_w_retention[length(mrip_w_retention)]

hbt_e_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,7), 6:dim(base$sizeselex)[2]]
#hbt_e_retention_multiplier<-hbt_e_retention[length(hbt_e_retention)]
hbt_w_retention_len<-base$sizeselex[(base$sizeselex$Yr==yr_end) & is.element(base$sizeselex$Factor,"Ret") & is.element(base$sizeselex$Fleet,8), 6:dim(base$sizeselex)[2]]
#hbt_w_retention_multiplier<-hbt_w_retention[length(hbt_w_retention)]

#fisheries status
#Catch number in the last three years unit 1000s
#hl_e_pred_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'retain(N):_1']))
#hl_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'retain(N):_2']))
#ll_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'retain(N):_3']))
#ll_w_pred_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'retain(N):_4']))
#mrip_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'retain(N):_5']))
#mrip_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'retain(N):_6']))
#hbt_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'retain(N):_7']))
#hbt_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'retain(N):_8']))
hl_e_pred_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_1']))
hl_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_2']))
ll_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_3']))
ll_w_pred_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_4']))
mrip_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_5']))
mrip_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_6']))
hbt_e_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_7']))
hbt_w_pred_cat_N  <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_8']))
comm_closed_e_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_9']))
comm_closed_w_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_10']))
rec_closed_e_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_11']))
rec_closed_w_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_12']))
shrimp_e_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,1) & is.element(base$timeseries$Era,'TIME'),'dead(N):_13']))
shrimp_w_cat_N <- c(t(base$timeseries[(base$timeseries$Yr>=(yr_end-2)) & (base$timeseries$Yr<=yr_end) & is.element(base$timeseries$Area,2) & is.element(base$timeseries$Era,'TIME'),'dead(N):_14']))

#Ratio_hl_east<-mean(hl_e_pred_cat_N/(hl_e_pred_cat_N+ll_e_pred_cat_N+hl_w_pred_cat_N+ll_w_pred_cat_N))
#Ratio_ll_east<-mean(ll_e_pred_cat_N/(hl_e_pred_cat_N+ll_e_pred_cat_N+hl_w_pred_cat_N+ll_w_pred_cat_N))
#Ratio_hl_west<-mean(hl_w_pred_cat_N/(hl_e_pred_cat_N+ll_e_pred_cat_N+hl_w_pred_cat_N+ll_w_pred_cat_N))
#Ratio_ll_west<-1-Ratio_hl_east-Ratio_ll_east-Ratio_hl_west

stock_1_N<-base$natage[(base$natage$Area==1)&(base$natage$Yr>=(yr_end-2))&(base$natage$Yr<=yr_end)&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]
stock_2_N<-base$natage[(base$natage$Area==2)&(base$natage$Yr>=(yr_end-2))&(base$natage$Yr<=yr_end)&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]

total_catch_N <-  (hl_e_pred_cat_N + hl_w_pred_cat_N + ll_e_pred_cat_N + ll_w_pred_cat_N
                   + mrip_e_pred_cat_N + mrip_w_pred_cat_N + hbt_e_pred_cat_N + hbt_w_pred_cat_N
                   + comm_closed_e_cat_N + comm_closed_w_cat_N + rec_closed_e_cat_N + rec_closed_w_cat_N
                   + shrimp_e_cat_N + shrimp_w_cat_N)

sum_N <- c(t(rowSums(stock_1_N)+rowSums(stock_2_N)))

Current_F<-mean(total_catch_N/sum_N) #same as Table 5.2
Current_SSB<-sum(spawning_output_1+spawning_output_2)

#Other parameters might be useful
midyearweight_at_age<-c(t(base$ageselex[(base$ageselex$Yr==yr_start)&is.element(base$ageselex$Factor,"bodywt")&is.element(base$ageselex$Fleet,"1"),as.character(base$agebins)]))
midyearweight_at_age_1<-midyearweight_at_age
midyearweight_at_age_2<-midyearweight_at_age

stock_1_mean_preyear<-c(t(base$natage[(base$natage$Area==1)&(base$natage$Yr==(yr_end-1))&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))
stock_2_mean_preyear<-c(t(base$natage[(base$natage$Area==2)&(base$natage$Yr==(yr_end-1))&(base$natage$`Beg/Mid`=="B"),as.character(base$agebins)]))

SSB_preyear_1<-stock_1_mean*fec_at_age_1
SSB_preyear_2<-stock_2_mean*fec_at_age_2

#Also get biological reference point
#SSB_MSY_BRP # of eggs, F_MSY_BRP per year
SSB_MSY_BRP<-1.23e+12 #unit 1000 eggs 
F_MSY_BRP<-0.0588

MSST<-0.5*SSB_MSY_BRP
MFMT<-F_MSY_BRP

#Position of the year 2016 in the traffic light
Current_F_ratio<-Current_F/MFMT
Current_SSB_ratio<-Current_SSB/SSB_MSY_BRP

#Options, in the webpage, these value should be read from the database
#Stock Assessment Model Input:
model_switch<-1 #newly added # 1 means SS

#General Input:
time_step_switch<-2 # 1: half year, 2: 1 year.  12312018 currently only year.
project_start_year<-2016 #  12312018 suppose can start from any date, but currently only from the beginning of the last stock assessment year
StockAssInterval<-5 # unit year
Runtime_long<-20 # unit year
plusage<-20 #currently can not change
age_1<-0:plusage # This is the largest age bin (plus age). Should be generate after "agebins" is defined
age_2<-age_1 # This is the largest age bin (plus age). Should be generate after "agebins" is defined
Simrun_Num<-100
IniN_Ess_Num<-1000
seed_switch<-2 #1: system clock 2: default, 3: self-defind
if(seed_switch==1){#use system clock
}else if(seed_switch==2){#use default file
  seed_input<-read.csv(file.choose(),header = F)
}else if(seed_switch==3){#use self uploaded path
  seed_input<-read.csv(file.choose(),header = F)
}

#Mixing Pattern:
mixing_pattern_switch<-1 # 1: no mixing, 2: constant, 3: same over year, 4: time varies 12312018 currently only no mixing

if(mixing_pattern_switch==1){
  #no mixing
}else if(mixing_pattern_switch==2){
  #use constant
}

#Initial Population:
#stock_1_mean<-stock_1_mean #read in from database
#stock_2_mean<-stock_2_mean #read in from database
stock_1_dis_temp<-stock_1_mean/sum(stock_1_mean) #Initial distribution, accumulative percentage
stock_2_dis_temp<-stock_2_mean/sum(stock_2_mean) #Initial distribution, accumulative percentage

stock_1_dis<-stock_1_dis_temp
for(i.dist in 2:length(stock_1_mean)){#
  stock_1_dis[i.dist]<-stock_1_dis[i.dist-1]+stock_1_dis_temp[i.dist]
}
stock_2_dis<-stock_2_dis_temp
for(i.dist in 2:length(stock_2_mean)){#
  stock_2_dis[i.dist]<-stock_2_dis[i.dist-1]+stock_2_dis_temp[i.dist]
}

cv_N_1<-0.2
cv_N_2<-0.2

#Biological Parameters:
weight_at_age_1<-weight_at_age_1
fec_at_age_1<-fec_at_age_1
weight_at_age_2<-weight_at_age_2
fec_at_age_2<-fec_at_age_2

#Natural Mortality:
M_switch<-3 #1: high M, 2: low M, and 3: default current M
if(M_switch==1){
  M_1<-M_1*1.5
  M_2<-M_2*1.5
}else if(M_switch==2){
  M_1<-M_1*0.5
  M_2<-M_2*0.5
}else if(M_switch==3){
  M_1<-M_1
  M_2<-M_2
}

cv_M_1<-0.2
cv_M_2<-0.2

#Recruitment
season_factor<-0.5 #unit year
sigma_R<-0.3
stock_1_rec_ratio<-0.23
stock_2_rec_ratio<-1-stock_1_rec_ratio
Rratio_Ess_Num<-1000

#the following are interactions, can be done with javascript
#set switch for recruitment: switch_R_pattern 1: history  2: estimate from equation 
#set switch for recruitment: switch_R_formular_year & switch_R_history_year: 1: include years before 1984   2:exclude years before 1984
switch_R_pattern<-2 #corresponding to the "recruitTypeStock1" switch in the database
switch_formular_modified<-1 #corresponding to the "formulaStock1" switch in the database, currently no other choice
switch_R_formular_year<-2 #corresponding to the "fromFmlStock1" switch in the database
switch_R_history_year<-2 #corresponding to the "fromHisStock1" switch in the database
switch_percentile<-2 #corresponding to the "historySt1" switch in the database
switch_percentile_early<-2 #corresponding to the "historySt1_early" switch in the database

#assuming R unit 1000, the following four lines are get from the database.
R_mean_early<-146168.88
R_mean_late<-116746.32
R_percentile_early<-146168.88 
R_percentile_late<-116746.32

if(switch_R_pattern==1){
  if(switch_R_history_year==1){
    if(switch_percentile_early==1){
      Rhist<-R_mean_early
    }else{
      Rhist<-R_percentile_early
    }
  }else if(switch_R_history_year==2){
    if(switch_percentile==1){
      Rhist<-R_mean_late
    }else{
      Rhist<-R_percentile_late
    }
  }
}else if(switch_R_pattern==2){
  if(switch_formular_modified==1){
    if(switch_R_formular_year==1){
      R0<-R0_early #R0_early is saved in the database
    }else if (switch_R_formular_year==2){
      R0<-R0_late  #R0_late is saved in the database
    }
    steepness<-0.99 #from the database
    SSB0<-4.71534e+12 #from the database
  }
} 
  
#Management Options I:
SSB_MSY_BRP<-1.23e+12
F_MSY_BRP<-0.0588

HCR_pattern<-2 #HCR pattern 1: constant C, 2: constant F

Harvest_level<-0.0564 #FOY 0.0564 F26%SPR 0.0588

imple_error<-0.2

#Management option II:
Ratio_rec<-0.49
Ratio_comm<-1-Ratio_rec # ratio of recreational sector take from the annual catch target

#the following 4 lines are read in from the database 2014-2016, but the actually values come from stock assessment output, see previous lines 221-224. Only used for demostrate the results.
#Ratio_hl_east<-0.4036
#Ratio_ll_east<-0.0266
#Ratio_hl_west<-0.5657
#Ratio_ll_west<-0.0041

Forhire_ratio<-0.423
Private_ratio<-1-Forhire_ratio

Ratio_private_east<-0.802   #using data 2010-2014
Ratio_private_west<-1-Ratio_private_east

Ratio_forhire_east<-0.722   #using data 2010-2014
Ratio_forhire_west<-1-Ratio_forhire_east

#Ratio_headboat_E<-0.307     #using data 2010-2014
#Ratio_chartboat_E<-1-Ratio_headboat_E
#Ratio_headboat_W<-0.757     #using data 2010-2014
#Ratio_chartboat_W<-1-Ratio_headboat_W

ABC_pvalue<-0.4 #has changed to 0.4 in the interface

Comme_buffer<-0
Private_buffer<-0.2
Forhire_buffer<-0.09

#Management Options III:
#legal size
Comme_min_size_in<-13 # unit inch
Recre_min_size_in<-16  #unit inch

#bag limit
Private_bag_limit<- 2 # unit Number per bag
Forhire_bag_limit<- 2 # unit Number per bag

#Management Options IV:
#discard fractions unit % Shrimp E and Shrimp W are 1
Recre_discard_E_open<-0.118 #fleet 5 & 7
Recre_discard_E_closed<-0.118 #fleet 11
Recre_discard_W_open<-0.118 #fleet 6 & 8
Recre_discard_W_closed<-0.118 #fleet 12
Comme_hl_discard_E_open<-0.56 #fleet 1
Comme_hl_discard_W_open<-0.6 #fleet 2
Comme_ll_discard_E_open<-0.64 #fleet 3
Comme_ll_discard_W_open<-0.81 #fleet 4
Comme_discard_E_closed<-0.55 #fleet 9
Comme_discard_W_closed<-0.74 #fleet 10

#the following MRIP and HBT discard are determined by Recre_discard_E_open & Recre_discard_W_open
#Recre_MRIP_discard_E_open<-Recre_discard_E_open #fleet 5
#Recre_MRIP_discard_W_open<-Recre_discard_W_open #fleet 6
#Recre_HBT_discard_E_open<-Recre_discard_E_open #fleet 7
#Recre_HBT_discard_W_open<-Recre_discard_W_open #fleet 8

Recre_discard_E_open_CV<-0.05
Recre_discard_E_closed_CV<-0.05
Recre_discard_W_open_CV<-0.05
Recre_discard_W_closed_CV<-0.05
Comme_ll_discard_E_open_CV<-0.05
Comme_ll_discard_W_open_CV<-0.05
Comme_hl_discard_E_open_CV<-0.05
Comme_hl_discard_W_open_CV<-0.05
Comme_discard_E_closed_CV<-0.05
Comme_discard_W_closed_CV<-0.05

#the following MRIP and HBT discard are determined by Recre_discard_E_open & Recre_discard_W_open
#Recre_MRIP_discard_E_open_CV<-Recre_discard_E_open_CV #fleet 5
#Recre_MRIP_discard_W_open_CV<-Recre_discard_W_open_CV #fleet 6
#Recre_HBT_discard_E_open_CV<-Recre_discard_E_open_CV #fleet 8
#Recre_HBT_discard_W_open_CV<-Recre_discard_W_open_CV #fleet 8

#Mangement Option V: 
#Catch Rate
#get from the database
Forhire_federal_catch_rate_lb<-46077 #estimated considered the multiplier; unit lb per day
Private_AL_catch_rate_lb<-35225 #unit lb per day
Private_FL_catch_rate_lb<-50176 #unit lb per day
Private_LA_catch_rate_lb<-7783 #unit lb per day
Private_MS_catch_rate_lb<-1823 #unit lb per day
Forhire_MS_catch_rate_lb<-176 #unit lb per day
Private_TX_catch_rate_lb<-2880 #unit lb per day
Forhire_TX_catch_rate_lb<-225 #unit lb per day

Forhire_season_switch<-1
Private_season_switch<-1

if(Forhire_season_switch==1) {
  #determined by ACT
}  else if(Forhire_season_switch==2) {#input by user
  Forhire_season_length<-62 #unit days
}

Forhire_MS_season_length_mean<-17 #unit days
Forhire_TX_season_length<-365 #unit days for state water

if(Private_season_switch==1) {
  #season determined by ACT
}else if(Private_season_switch==2) {#input by user
  Private_AL_season_length<-36 #unit days
  Private_FL_season_length<-38 #unit days
  Private_LA_season_length<-109 #unit days
  Private_MS_season_length<-81 #unit days
  Private_TX_season_length<-62 #unit days
}

#Mangement Option VI:
Penalty_switch<-0  #0 no action, 1 pentalty to the next year
Carryover_switch<-1 #0 no action, 1 carryover to the next year but no over 95% of OFL, 2 to the next year no over 50% between ABC and OFL

#Options for output calculation from data base also:
Private_quota_AL<-0.26298
Private_quota_FL<-0.44822
Private_quota_LA<-0.1912
Private_quota_MS<-0.0355
Private_quota_TX<-0.0621

#Other parameters
#midyearweight_at_age_1<-midyearweight_at_age_1
#midyearweight_at_age_2<-midyearweight_at_age_2

SSB_preyear_1<-SSB_preyear_1
SSB_preyear_2<-SSB_preyear_2

#save.image("MSEspace1.RData")

#some calculation before the ACL simulation
#retention rate, it relates with the legal size, therefore, backcalculate the retention rate, ref min_size commerical 13 in (33.02 cm), recreational 16 in (40.64)
#although can use getparameters<-SS_parlines(dir = direct, ctlfile = "control.ss_new", version = "3.30", verbose = TRUE, active = FALSE) to get the parameters, it will be too complicated
ref_Comme_min_size_cm<-33.02 #unit cm
ref_Recre_min_size_cm<-40.64 #unit cm
Comme_min_size_cm<-Comme_min_size_in*2.54 # unit inch to cm
Recre_min_size_cm<-Recre_min_size_in*2.54  #unit inch to cm

#thre_bin_comme<-as.numeric(rownames(length_age_key))[1]+2-2*which.max(as.numeric(rownames(length_age_key))<=Comme_min_size_cm)
#thre_bin_recre<-as.numeric(rownames(length_age_key))[1]+2-2*which.max(as.numeric(rownames(length_age_key))<=Recre_min_size_cm)

act_hl_e_retention_len<-hl_e_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Comme_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Comme_min_size_cm)))
act_hl_e_retention<-c(as.matrix(act_hl_e_retention_len)%*%length_age_key_stock1[nrow(length_age_key_stock1):1,])
act_hl_w_retention_len<-hl_w_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Comme_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Comme_min_size_cm)))
act_hl_w_retention<-c(as.matrix(act_hl_w_retention_len)%*%length_age_key_stock2[nrow(length_age_key_stock2):1,])

act_ll_e_retention_len<-ll_e_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Comme_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Comme_min_size_cm)))
act_ll_e_retention<-c(as.matrix(act_ll_e_retention_len)%*%length_age_key_stock1[nrow(length_age_key_stock1):1,])
act_ll_w_retention_len<-ll_w_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Comme_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Comme_min_size_cm)))
act_ll_w_retention<-c(as.matrix(act_ll_w_retention_len)%*%length_age_key_stock2[nrow(length_age_key_stock2):1,])

act_mrip_e_retention_len<-mrip_e_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Recre_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Recre_min_size_cm)))
act_mrip_e_retention<-c(as.matrix(act_mrip_e_retention_len)%*%length_age_key_stock1[nrow(length_age_key_stock1):1,])
act_mrip_w_retention_len<-mrip_w_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Recre_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Recre_min_size_cm)))
act_mrip_w_retention<-c(as.matrix(act_mrip_w_retention_len)%*%length_age_key_stock2[nrow(length_age_key_stock2):1,])

act_hbt_e_retention_len<-hbt_e_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Recre_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Recre_min_size_cm)))
act_hbt_e_retention<-c(as.matrix(act_hbt_e_retention_len)%*%length_age_key_stock1[nrow(length_age_key_stock1):1,])
act_hbt_w_retention_len<-hbt_w_retention_len*(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-ref_Recre_min_size_cm)))/(1+exp(-(rev(as.numeric(rownames(length_age_key)))+1-Recre_min_size_cm)))
act_hbt_w_retention<-c(as.matrix(act_hbt_w_retention_len)%*%length_age_key_stock2[nrow(length_age_key_stock2):1,])

#estimate direct fishing F
comp_comm_sel_e <- hl_e_pred_F_ave * hl_e_selex * act_hl_e_retention + ll_e_pred_F_ave * ll_e_selex * act_ll_e_retention 
comp_comm_sel_w <- hl_w_pred_F_ave * hl_w_selex * act_hl_w_retention + ll_w_pred_F_ave * ll_w_selex * act_ll_w_retention

comp_recr_sel_e <- mrip_e_pred_F_ave * mrip_e_selex * act_mrip_e_retention + hbt_e_pred_F_ave* hbt_e_selex * act_hbt_e_retention
comp_recr_sel_w <- mrip_w_pred_F_ave * mrip_w_selex * act_mrip_w_retention + hbt_w_pred_F_ave* hbt_w_selex * act_hbt_w_retention

#estimate open season discards for a directed fleet
comp_comm_death_discard_e <- hl_e_pred_F_ave * hl_e_selex * (1-act_hl_e_retention) * Comme_hl_discard_E_open + ll_e_pred_F_ave * ll_e_selex * (1-act_ll_e_retention) * Comme_ll_discard_E_open
comp_comm_death_discard_w <- hl_w_pred_F_ave * hl_w_selex * (1-act_hl_w_retention) * Comme_hl_discard_W_open + ll_w_pred_F_ave * ll_w_selex * (1-act_ll_w_retention) * Comme_ll_discard_W_open

comp_recr_death_discard_e <- mrip_e_pred_F_ave * mrip_e_selex * (1-act_mrip_e_retention) * Recre_discard_E_open + hbt_e_pred_F_ave* hbt_e_selex * (1-act_hbt_e_retention) * Recre_discard_E_open
comp_recr_death_discard_w <- mrip_w_pred_F_ave * mrip_w_selex * (1-act_mrip_w_retention) * Recre_discard_W_open + hbt_w_pred_F_ave* hbt_w_selex * (1-act_hbt_w_retention) * Recre_discard_W_open

comp_bycatch_closed_e <- rec_closed_e_pred_F_ave * rec_closed_e_selex * Recre_discard_E_closed + comm_closed_e_pred_F_ave * comm_closed_e_selex * Comme_discard_E_closed + shrimp_e_pred_F_ave * shrimp_e_selex
comp_bycatch_closed_w <- rec_closed_w_pred_F_ave * rec_closed_w_selex * Recre_discard_W_closed + comm_closed_w_pred_F_ave * comm_closed_w_selex * Comme_discard_W_closed + shrimp_w_pred_F_ave * shrimp_w_selex 

#estimate F functions
Function_relF<-function(F0) (true_quote_short-sum(true_N_1_short*(1-exp(-F0*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_death_discard_e+comp_recr_death_discard_e) - comp_bycatch_closed_e - M_1_true))*
                                              (F0*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_death_discard_e+comp_recr_death_discard_e) + comp_bycatch_closed_e)/
                                              (F0*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_death_discard_e+comp_recr_death_discard_e) + comp_bycatch_closed_e + M_1_true)) -
                                              sum(true_N_2_short*(1-exp(-F0*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_death_discard_w+comp_recr_death_discard_w) - comp_bycatch_closed_w - M_2_true))*
                                              (F0*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_death_discard_w+comp_recr_death_discard_w) + comp_bycatch_closed_w)/
                                              (F0*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_death_discard_w+comp_recr_death_discard_w) + comp_bycatch_closed_w + M_2_true)))^2

Function_comm_imple<-function(F1)(Comme_planned_catch-sum(true_N_1*(1-exp(-F1*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_discard_e_death_open_true +comp_recr_discard_e_death_open_true) - comp_discard_e_death_closed_true - comp_bycatch_e_true - M_1_true))*
                                                            (F1*(comp_comm_sel_e))/
                                                            (F1*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_discard_e_death_open_true +comp_recr_discard_e_death_open_true) + comp_discard_e_death_closed_true + comp_bycatch_e_true + M_1_true) * weight_at_age_1 + 
                                                            true_N_2*(1-exp(-F1*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_discard_w_death_open_true +comp_recr_discard_w_death_open_true) - comp_discard_w_death_closed_true - comp_bycatch_w_true - M_2_true))*
                                                            (F1*(comp_comm_sel_w))/
                                                            (F1*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_discard_w_death_open_true +comp_recr_discard_w_death_open_true) + comp_discard_w_death_closed_true + comp_bycatch_w_true + M_2_true) * weight_at_age_2))^2

Function_recr_imple<-function(F2)(Recre_planned_catch-sum(true_N_1*(1-exp(-F2*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_discard_e_death_open_true +comp_recr_discard_e_death_open_true) - comp_discard_e_death_closed_true - comp_bycatch_e_true - M_1_true))*
                                                            (F2*(comp_recr_sel_e))/
                                                            (F2*(comp_comm_sel_e+comp_recr_sel_e + comp_comm_discard_e_death_open_true +comp_recr_discard_e_death_open_true) + comp_discard_e_death_closed_true + comp_bycatch_e_true + M_1_true) * weight_at_age_1 +
                                                            true_N_2*(1-exp(-F2*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_discard_w_death_open_true +comp_recr_discard_w_death_open_true) - comp_discard_w_death_closed_true - comp_bycatch_w_true - M_2_true))*
                                                            (F2*(comp_recr_sel_w))/
                                                            (F2*(comp_comm_sel_w+comp_recr_sel_w + comp_comm_discard_w_death_open_true +comp_recr_discard_w_death_open_true) + comp_discard_w_death_closed_true + comp_bycatch_w_true + M_2_true) * weight_at_age_2))^2

Runtime_short<-StockAssInterval
Simrun_Num_short<-100
####Below is to project for a short term to determin ABC/ACL: assume the name of the parameters are the same. Run 100 times
Catch_pred<-matrix(rep(0,Simrun_Num_short*Runtime_short), ncol=Runtime_short)
ass_err<-0.05
#save.image("test1.RData")

##true accountability measures
#R_1<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#R_2<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#SSB_1<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#SSB_2<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#total_SSB<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#totalcatch_1_imple<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#totalcatch_2_imple<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#total_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#comm_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#comm_catch_1<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#comm_catch_2<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#recr_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Forhire_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Forhire_catch_1<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Forhire_catch_2<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Private_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Private_catch_1<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#Private_catch_2<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#comm_discards<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#recr_discards<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#comm_dislandratio<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#recr_dislandratio<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_fed_forhire_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_AL_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_FL_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_LA_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_MS_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
##true_forhire_MS_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_TX_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
##true_forhire_TX_season_length<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_AL_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_FL_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_LA_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_MS_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
##true_forhire_MS_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#true_private_TX_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
##true_forhire_TX_catch<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#F_general<-matrix(rep(0,Simrun_Num*(Runtime_long+2)), ncol=(Runtime_long+2))
##AM_comm_overunderage<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#AM_comm_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_forhire_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_AL_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_FL_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_LA_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_MS_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
#AM_TX_overunderage<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1))
##AM_forhire_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
##AM_AL_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
##AM_FL_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
##AM_LA_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
##AM_MS_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
##AM_TX_flag<-matrix(rep(0,Simrun_Num*(Runtime_long+1)), ncol=(Runtime_long+1)) #default 0, 1 overage -1 underage
#Year_green<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#F_ratio<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)
#SSB_ratio<-matrix(rep(0,Simrun_Num*Runtime_long), ncol=Runtime_long)

#parallel temp
R_1_paratemp<-rep(0,(Runtime_long+1))
R_2_paratemp<-rep(0,(Runtime_long+1))
SSB_1_paratemp<-rep(0,(Runtime_long+1))
SSB_2_paratemp<-rep(0,(Runtime_long+1))
total_SSB_paratemp<-rep(0,Runtime_long)
totalcatch_1_imple_paratemp<-rep(0,Runtime_long)
totalcatch_2_imple_paratemp<-rep(0,Runtime_long)
total_catch_paratemp<-rep(0,Runtime_long)
comm_catch_paratemp<-rep(0,Runtime_long)
comm_catch_1_paratemp<-rep(0,Runtime_long)
comm_catch_2_paratemp<-rep(0,Runtime_long)
recr_catch_paratemp<-rep(0,Runtime_long)
Forhire_catch_paratemp<-rep(0,Runtime_long)
Forhire_catch_1_paratemp<-rep(0,Runtime_long)
Forhire_catch_2_paratemp<-rep(0,Runtime_long)
Private_catch_paratemp<-rep(0,Runtime_long)
Private_catch_1_paratemp<-rep(0,Runtime_long)
Private_catch_2_paratemp<-rep(0,Runtime_long)
comm_discards_paratemp<-rep(0,Runtime_long)
recr_discards_paratemp<-rep(0,Runtime_long)
comm_dislandratio_paratemp<-rep(0,Runtime_long)
recr_dislandratio_paratemp<-rep(0,Runtime_long)
true_fed_forhire_season_length_paratemp<-rep(0,Runtime_long)
true_private_AL_season_length_paratemp<-rep(0,Runtime_long)
true_private_FL_season_length_paratemp<-rep(0,Runtime_long)
true_private_LA_season_length_paratemp<-rep(0,Runtime_long)
true_private_MS_season_length_paratemp<-rep(0,Runtime_long)
true_private_TX_season_length_paratemp<-rep(0,Runtime_long)
true_private_AL_catch_paratemp<-rep(0,Runtime_long)
true_private_FL_catch_paratemp<-rep(0,Runtime_long)
true_private_LA_catch_paratemp<-rep(0,Runtime_long)
true_private_MS_catch_paratemp<-rep(0,Runtime_long)
true_private_TX_catch_paratemp<-rep(0,Runtime_long)
F_general_paratemp<-rep(0,(Runtime_long+2))
AM_comm_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_forhire_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_AL_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_FL_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_LA_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_MS_overunderage_paratemp<-rep(0,(Runtime_long+1))
AM_TX_overunderage_paratemp<-rep(0,(Runtime_long+1))
Year_green_paratemp<-rep(0,Runtime_long)
F_ratio_paratemp<-rep(0,Runtime_long)
SSB_ratio_paratemp<-rep(0,Runtime_long)

imple_error2<-0.05 #CV for catch rate

switch_standard<-1

#set up parallel 
no_cores <- detectCores()  #determine number of cores on computer
cl<-makeCluster(no_cores)  #setup the size of your parallel cluster based on number of cores
registerDoSNOW(cl)         #register cluster to get it ready for parallel processing

Sys.time()
ls=foreach(i.run=1:Simrun_Num) %dopar% {
  
  if(seed_switch==2|seed_switch==3){#use seed
    set.seed(as.numeric(seed_input[i.run]))
  }
  
  #initial N Distribution, add random noise
  temp1<-runif(IniN_Ess_Num)
  temp2<-runif(IniN_Ess_Num)
  count_1_imple<-rep(0,length(stock_1_mean))
  for (i.ess in 1:IniN_Ess_Num){
    j<-which.min(temp1[i.ess]>stock_1_dis)
    count_1_imple[j]<-count_1_imple[j]+1
  }
  true_N_1_dis<-count_1_imple/sum(count_1_imple)
  
  count_2_imple<-rep(0,length(stock_2_mean))
  for (i.ess in 1:IniN_Ess_Num){
    j<-which.min(temp2[i.ess]>stock_2_dis)
    count_2_imple[j]<-count_2_imple[j]+1
  }
  true_N_2_dis<-count_2_imple/sum(count_2_imple)
  
  #initial N unit 1000s, add random noise
  mean_true_N_1_total<-log(sum(stock_1_mean))-log(1+cv_N_1^2)/2
  sd_true_N_1_total<-sqrt(log(1+cv_N_1^2))
  true_N_1_total<-rlnorm(1,mean_true_N_1_total,sd_true_N_1_total)
  while ((true_N_1_total>sum(stock_1_mean)*2)|(true_N_1_total<sum(stock_1_mean)/2)){
    true_N_1_total<-rlnorm(1,mean_true_N_1_total,sd_true_N_1_total)
  }
  
  mean_true_N_2_total<-log(sum(stock_2_mean))-log(1+cv_N_2^2)/2
  sd_true_N_2_total<-sqrt(log(1+cv_N_2^2))
  true_N_2_total<-rlnorm(1,mean_true_N_2_total,sd_true_N_2_total)
  while ((true_N_2_total>sum(stock_2_mean)*2)|(true_N_2_total<sum(stock_2_mean)/2)){
    true_N_2_total<-rlnorm(1,mean_true_N_2_total,sd_true_N_2_total)
  }
  
  true_N_1<-true_N_1_total*true_N_1_dis
  true_N_2<-true_N_2_total*true_N_2_dis
  
  R_1_paratemp[1]<-true_N_1[1] # initial recruitment unit 1000s
  R_2_paratemp[1]<-true_N_2[1] # initial recruitment unit 1000s
  
  SSB_1_paratemp[1]<-sum(SSB_preyear_1) # initial SSB unit 1000 eggs
  SSB_2_paratemp[1]<-sum(SSB_preyear_2) # initial SSB unit 1000 eggs
  
  F_general_paratemp[1]<-0.041
  F_general_paratemp[2]<-0.052
  
  overunderage<-0
  
  for (i.runtime in 1:Runtime_long){# Both long term and short term projections are all 100 years, "short" only means OFL and ABC projection.
   
    ##estimate recruitment
    #estimate SSB for the next year
    true_SSB_1<-true_N_1 * fec_at_age_1
    true_SSB_2<-true_N_2 * fec_at_age_2
    
    SSB_1_paratemp[i.runtime+1]<-sum(true_SSB_1)
    SSB_2_paratemp[i.runtime+1]<-sum(true_SSB_2)
    
    if(switch_R_pattern==1){
      tempR_mean<-as.numeric(Rhist)
    }else if(switch_R_pattern==2){
      tempR_mean<-4*steepness*R0*(SSB_1_paratemp[i.runtime+1]+SSB_2_paratemp[i.runtime+1])/((1-steepness)*SSB0+(5*steepness-1)*(SSB_1_paratemp[i.runtime+1]+SSB_2_paratemp[i.runtime+1]))
    } #R unit 1000s, SSB unit 1000 eggs
    
    #sigma_R here is the CV, sorry for the confusing name
    mean_tempR<-log(tempR_mean)-log(1+sigma_R^2)/2
    sd_tempR<-sqrt(log(1+sigma_R^2))
    tempR<-rlnorm(1,mean_tempR,sd_tempR)
    while ((tempR>tempR_mean*2)|(tempR<tempR_mean/2)){
      tempR<-rlnorm(1,mean_tempR,sd_tempR)
    }
    
    count_3_imple<-0
    temp3<-runif(Rratio_Ess_Num)
    for (i.ess2 in 1:Rratio_Ess_Num){
      if(temp3[i.ess2]<stock_1_rec_ratio)
        count_3_imple<-count_3_imple+1
    }
    true_stock_1_rec_ratio<-count_3_imple/Rratio_Ess_Num
    
    true_R_1<-tempR*stock_1_rec_ratio
    true_R_2<-tempR*stock_2_rec_ratio
    
    R_1_paratemp[i.runtime+1]<-true_R_1
    R_2_paratemp[i.runtime+1]<-true_R_2
    
    B_at_age_1<-true_N_1*weight_at_age_1 #biomass unit 1000*kg=mt
    B_at_age_2<-true_N_2*weight_at_age_2 #biomass unit 1000*kg=mt
    
    #Check whether it is a stock assessment year to update ABC or not
    if(i.runtime%%StockAssInterval==1){
      if(switch_standard==1){#standard version use short-cut MSE
        stock_1_mean_short<-true_N_1
        stock_2_mean_short<-true_N_2
          
        #N distribution
        stock_1_dis_temp<-stock_1_mean_short/sum(stock_1_mean_short)
        stock_1_dis_short<-stock_1_dis_temp
        for(i.dist in 2:length(stock_1_mean_short)){#accumulative percentage
          stock_1_dis_short[i.dist]<-stock_1_dis_short[i.dist-1]+stock_1_dis_temp[i.dist]
        }
        
        stock_2_dis_temp<-stock_2_mean_short/sum(stock_2_mean_short)
        stock_2_dis_short<-stock_2_dis_temp
        for(i.dist in 2:length(stock_2_mean_short)){#accumulative percentage
          stock_2_dis_short[i.dist]<-stock_2_dis_short[i.dist-1]+stock_2_dis_temp[i.dist]
        }
        
        for(i.run_short in 1:Simrun_Num_short){
          #initial N Distribution
          #if(seed_switch==2|seed_switch==3){#use seed
          #  set.seed(as.numeric(seed_input[100+i.run_short]))
          #}
          
          temp1<-runif(IniN_Ess_Num)
          temp2<-runif(IniN_Ess_Num)
          
          count_1<-rep(0,length(stock_1_mean_short))
          for (i.ess in 1:IniN_Ess_Num){
            j<-which.min(temp1[i.ess]>stock_1_dis_short)
            count_1[j]<-count_1[j]+1
          }
          true_N_1_dis_short<-count_1/sum(count_1)
          
          count_2<-rep(0,length(stock_2_mean_short))
          for (i.ess in 1:IniN_Ess_Num){
            j<-which.min(temp2[i.ess]>stock_2_dis_short)
            count_2[j]<-count_2[j]+1
          }
          true_N_2_dis_short<-count_2/sum(count_2)
          
          #initial N unit 1000s
          mean_true_N_1_total_short<-log(sum(stock_1_mean_short))-log(1+ass_err^2)/2
          sd_true_N_1_total_short<-sqrt(log(1+ass_err^2))
          true_N_1_total_short<-rlnorm(1,mean_true_N_1_total_short,sd_true_N_1_total_short)
          while ((true_N_1_total_short>sum(stock_1_mean_short)*2)|(true_N_1_total_short<sum(stock_1_mean_short)/2)){
            true_N_1_total_short<-rlnorm(1,mean_true_N_1_total_short,sd_true_N_1_total_short)
          }
          
          mean_true_N_2_total_short<-log(sum(stock_2_mean_short))-log(1+ass_err^2)/2
          sd_true_N_2_total_short<-sqrt(log(1+ass_err^2))
          true_N_2_total_short<-rlnorm(1,mean_true_N_2_total_short,sd_true_N_2_total_short)
          while ((true_N_2_total_short>sum(stock_2_mean_short)*2)|(true_N_2_total_short<sum(stock_2_mean_short)/2)){
            true_N_2_total_short<-rlnorm(1,mean_true_N_2_total_short,sd_true_N_2_total_short)
          }
          
          true_N_1_short<-true_N_1_total_short*true_N_1_dis_short
          true_N_2_short<-true_N_2_total_short*true_N_2_dis_short
          
          for (i.runtime_short in 1:Runtime_short){#Runtime_short just mean shot-term pre-projection years.
            ##estimate recruitment
            #estimate SSB in the next year
            true_SSB_1_short<-true_N_1_short * fec_at_age_1
            true_SSB_2_short<-true_N_2_short * fec_at_age_2
            
            if(switch_R_pattern==1){
              tempR_mean_short<-as.numeric(Rhist) # input
            }else if (switch_R_pattern==2){
              tempR_mean_short<-4*steepness*R0*sum(true_SSB_1_short+true_SSB_2_short)/((1-steepness)*SSB0+(5*steepness-1)*sum(true_SSB_1_short+true_SSB_2_short))
            }#R unit 1000s SSB unit 1000 eggs
            
            #here sigma_R is the CV, sorry for the confusing name
            mean_tempR_short<-log(tempR_mean_short)-log(1+sigma_R^2)/2
            sd_tempR_short<-sqrt(log(1+sigma_R^2))
            tempR_short<-rlnorm(1,mean_tempR_short,sd_tempR_short)
            while ((tempR_short>tempR_mean_short*2)|(tempR_short<tempR_mean_short/2)){
              tempR_short<-rlnorm(1,mean_tempR_short,sd_tempR_short)
            }
            
            count_4<-0
            temp4<-runif(Rratio_Ess_Num)
            for (i.ess2 in 1:Rratio_Ess_Num){
              if(temp4[i.ess2]<stock_1_rec_ratio)
                count_4<-count_4+1
            }
            true_stock_1_rec_ratio_short<-count_4/Rratio_Ess_Num
            
            pred_R_1_short<-tempR_short*true_stock_1_rec_ratio_short
            pred_R_2_short<-tempR_short*(1-true_stock_1_rec_ratio_short)
            
            B_at_age_1_short<-true_N_1_short*weight_at_age_1 #biomass unit 1000*kg=mt
            B_at_age_2_short<-true_N_2_short*weight_at_age_2 #biomass unit 1000*kg=mt
            
            #Catch quote unit number which includes bycatch and discards
            Quote_total_short<-(sum(true_N_1_short)+sum(true_N_2_short))*Harvest_level
            mean_true_quote_short<-log(Quote_total_short)-log(1+imple_error^2)/2
            sd_true_quote_short<-sqrt(log(1+imple_error^2))
            true_quote_short<-rlnorm(1,mean_true_quote_short,sd_true_quote_short)
            while ((true_quote_short>Quote_total_short*2)|(true_quote_short<Quote_total_short/2)){
              true_quote_short<-rlnorm(1,mean_true_quote_short,sd_true_quote_short)
            }
            
            #estimate M
            mean_M_1_relative<-log(1)-log(1+cv_M_1^2)/2
            sd_M_1_relative<-sqrt(log(1+cv_M_1^2))
            M_1_relative<-rlnorm(1,mean_M_1_relative,sd_M_1_relative)
            while ((M_1_relative>2)|(M_1_relative<1/2)){
              M_1_relative<-rlnorm(1,mean_M_1_relative,sd_M_1_relative)
            }
            
            mean_M_2_relative<-log(1)-log(1+cv_M_2^2)/2
            sd_M_2_relative<-sqrt(log(1+cv_M_2^2))
            M_2_relative<-rlnorm(1,mean_M_2_relative,sd_M_2_relative)
            while ((M_2_relative>2)|(M_2_relative<1/2)){
              M_2_relative<-rlnorm(1,mean_M_2_relative,sd_M_2_relative)
            }
            M_1_true<-M_1_relative*M_1
            M_2_true<-M_2_relative*M_2
            
            relF<-max(optimize(Function_relF,lower=0,upper=1000,tol = 0.01)$minimum,0)
            
            #optim(0,Function_relF,method="Nelder-Mead",hessian=TRUE)$par
            
            total_F_1<-relF * comp_recr_sel_e + relF * comp_comm_sel_e + 
              relF * comp_recr_death_discard_e + relF * comp_comm_death_discard_e + 
              comp_bycatch_closed_e #retention 0 discard 1 for shirmp
            
            total_F_2<-relF * comp_recr_sel_w + relF * comp_comm_sel_w + 
              relF * comp_recr_death_discard_w + relF * comp_comm_death_discard_w + 
              comp_bycatch_closed_w
            
            mod_F_1<-relF * comp_recr_sel_e + relF * comp_comm_sel_e
            
            mod_F_2<-relF * comp_recr_sel_w + relF * comp_comm_sel_w
            
            true_N_1_nextyeartemp_short<-true_N_1_short*exp(-total_F_1 -M_1_true)
            true_N_2_nextyeartemp_short<-true_N_2_short*exp(-total_F_2 -M_2_true)
            
            catch_N_1_pred<-(true_N_1_short-true_N_1_nextyeartemp_short)*mod_F_1/(total_F_1 + M_1_true)
            catch_N_2_pred<-(true_N_2_short-true_N_2_nextyeartemp_short)*mod_F_2/(total_F_2 + M_2_true)
            
            total_catch_B_pred<-sum(catch_N_1_pred * weight_at_age_1 +catch_N_2_pred * weight_at_age_2) #unit mt
            Catch_pred[i.run_short,i.runtime_short]<-total_catch_B_pred
            
            true_N_1_short[2:(length(true_N_1_short)-1)]<-true_N_1_nextyeartemp_short[1:(length(true_N_1_short)-2)]
            true_N_2_short[2:(length(true_N_2_short)-1)]<-true_N_2_nextyeartemp_short[1:(length(true_N_2_short)-2)]
            true_N_1_short[length(true_N_1_short)]<-true_N_1_nextyeartemp_short[length(true_N_1_short)-1]+true_N_1_nextyeartemp_short[length(true_N_1_short)]
            true_N_2_short[length(true_N_2_short)]<-true_N_2_nextyeartemp_short[length(true_N_2_short)-1]+true_N_2_nextyeartemp_short[length(true_N_2_short)]
            
            true_N_1_short[1]<-pred_R_1_short
            true_N_2_short[1]<-pred_R_2_short
          }
        }
      }else{#professional version use full MSE
        #will call SS3
      }
      #Estimate commercial and recreational related ABC and OFL, or commercial and recreational related ACL
      ABC_planned<-as.numeric(quantile(rowMeans(Catch_pred),ABC_pvalue)) #unit 1000*kg=mt
      OFL_planned<-as.numeric(quantile(rowMeans(Catch_pred),0.5)) 
      #save.image("MSEspace2.RData")
    }
   
    if(overunderage>0){ #has carryover. Here, ACL_planned is considered as a new ABC_planned. "overunderage" is the carryover_switch
      if(Carryover_switch==1){
        ACL_planned<-max(min((ABC_planned+overunderage),OFL_planned*0.95),0)
      }else if(Carryover_switch==2){
        ACL_planned<-max(min((ABC_planned+overunderage),(OFL_planned+ABC_planned)/2),0)
      }else if(Carryover_switch==0){
        ACL_planned<-max(ABC_planned,0)
      }
      
      sum_overunderage_lastyear<-max(AM_comm_overunderage_paratemp[i.runtime],0)+max(AM_forhire_overunderage_paratemp[i.runtime],0)+max(AM_AL_overunderage_paratemp[i.runtime],0)+max(AM_FL_overunderage_paratemp[i.runtime],0)+max(AM_LA_overunderage_paratemp[i.runtime],0)+max(AM_MS_overunderage_paratemp[i.runtime],0)+max(AM_TX_overunderage_paratemp[i.runtime],0)
      
      comm_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_comm_overunderage_paratemp[i.runtime],0)
      forhire_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_forhire_overunderage_paratemp[i.runtime],0)
      AL_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_AL_overunderage_paratemp[i.runtime],0)
      FL_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_FL_overunderage_paratemp[i.runtime],0)
      LA_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_LA_overunderage_paratemp[i.runtime],0)
      MS_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_MS_overunderage_paratemp[i.runtime],0)
      TX_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear+0.000001)*max(AM_TX_overunderage_paratemp[i.runtime],0)
      
    }else if(overunderage<=0){#has penalty. Currently no penalty yet. This part is reserved for the future expansion.
      if(Penalty_switch==1){
        if((yr_start+i.runtime-1)<=2032){
          ACL_planned<-max(min((ABC_planned+overunderage),OFL_planned*0.95),0)
        }else if(SSB_ratio_paratemp[(i.runtime-1)]<0.5){
          ACL_planned<-max(min((ABC_planned+overunderage),OFL_planned*0.95),0)
        }
      }else{
      ACL_planned<-max(ABC_planned,0) 
      }
      #ACL_planned<-max(ABC_planned,0) 
      
      comm_overunderage_lastyear<-0
      sum_overunderage_lastyear<-min(AM_forhire_overunderage_paratemp[i.runtime],0)+min(AM_AL_overunderage_paratemp[i.runtime],0)+min(AM_FL_overunderage_paratemp[i.runtime],0)+min(AM_LA_overunderage_paratemp[i.runtime],0)+min(AM_MS_overunderage_paratemp[i.runtime],0)+min(AM_TX_overunderage_paratemp[i.runtime],0)
      forhire_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_forhire_overunderage_paratemp[i.runtime],0)
      AL_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_AL_overunderage_paratemp[i.runtime],0)
      FL_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_FL_overunderage_paratemp[i.runtime],0)
      LA_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_LA_overunderage_paratemp[i.runtime],0)
      MS_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_MS_overunderage_paratemp[i.runtime],0)
      TX_overunderage_lastyear<-(ACL_planned-ABC_planned)/(sum_overunderage_lastyear-0.000001)*min(AM_TX_overunderage_paratemp[i.runtime],0)
    }
    
    Quote_comm_befbuf<-ACL_planned * Ratio_comm #unit 1000*kg=mt
    Quote_rec_befbuf<-ACL_planned * Ratio_rec #unit 1000*kg=mt
    
    Quote_forhire_befbuf<-Quote_rec_befbuf * Forhire_ratio #unit 1000*kg=mt
    Quote_private_befbuf<-Quote_rec_befbuf * Private_ratio #unit 1000*kg=mt
    
    Quote_comm_aftbuf<-Quote_comm_befbuf*(1-Comme_buffer)
    Quote_forhire_aftbuf<-Quote_forhire_befbuf*(1-Forhire_buffer)
    Quote_private_aftbuf<-Quote_private_befbuf*(1-Private_buffer)
    
    Quote_private_aftbuf_AL<-Private_quota_AL*Quote_private_aftbuf
    Quote_private_aftbuf_FL<-Private_quota_FL*Quote_private_aftbuf
    Quote_private_aftbuf_LA<-Private_quota_LA*Quote_private_aftbuf
    Quote_private_aftbuf_MS<-Private_quota_MS*Quote_private_aftbuf
    Quote_private_aftbuf_TX<-Private_quota_TX*Quote_private_aftbuf
    
    #commercial catch quota determined by ACT # set a lower implementation error for commercial catch?
    mean_Plan_quote_comm<-log(Quote_comm_aftbuf)-log(1+imple_error^2)/2
    sd_Plan_quote_comm<-sqrt(log(1+imple_error^2))
    Plan_quote_comm<-rlnorm(1,mean_Plan_quote_comm,sd_Plan_quote_comm)
    while((Plan_quote_comm>Quote_comm_aftbuf*2)|(Plan_quote_comm<Quote_comm_aftbuf/2)){
      Plan_quote_comm<-rlnorm(1,mean_Plan_quote_comm,sd_Plan_quote_comm)
    }
    Plan_quote_comm<-min(Plan_quote_comm,Quote_comm_aftbuf)
    
    Forhire_federal_catch_rate_kg<-Forhire_federal_catch_rate_lb*0.453592 #unit kg per day
    mean_Plan_Forhire_federal_catch_rate_kg<-log(Forhire_federal_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Forhire_federal_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Forhire_federal_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_federal_catch_rate_kg,sd_Plan_Forhire_federal_catch_rate_kg)
    while((Plan_Forhire_federal_catch_rate_kg>Forhire_federal_catch_rate_kg*2)|(Plan_Forhire_federal_catch_rate_kg<Forhire_federal_catch_rate_kg/2)){
      Plan_Forhire_federal_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_federal_catch_rate_kg,sd_Plan_Forhire_federal_catch_rate_kg)
    }
    
    if(Forhire_season_switch==1) {
      #determined by ACT
      mean_Plan_quote_forhire<-log(Quote_forhire_aftbuf)-log(1+imple_error^2)/2
      sd_Plan_quote_forhire<-sqrt(log(1+imple_error^2))
      Plan_quote_forhire<-rlnorm(1,mean_Plan_quote_forhire,sd_Plan_quote_forhire)
      while((Plan_quote_forhire>Quote_forhire_aftbuf*2)|(Plan_quote_forhire<Quote_forhire_aftbuf/2)){
        Plan_quote_forhire<-rlnorm(1,mean_Plan_quote_forhire,sd_Plan_quote_forhire)
      }
    }else if(Forhire_season_switch==2) {#input by user
      Plan_quote_forhire<-Forhire_season_length*Plan_Forhire_federal_catch_rate_kg/1000 #unit 1000 kg
    }
    
    Private_AL_catch_rate_kg<-Private_AL_catch_rate_lb*0.453592
    Private_FL_catch_rate_kg<-Private_FL_catch_rate_lb*0.453592
    Private_LA_catch_rate_kg<-Private_LA_catch_rate_lb*0.453592
    Private_MS_catch_rate_kg<-Private_MS_catch_rate_lb*0.453592
    Forhire_MS_catch_rate_kg<-Forhire_MS_catch_rate_lb*0.453592
    Private_TX_catch_rate_kg<-Private_TX_catch_rate_lb*0.453592
    Forhire_TX_catch_rate_kg<-Forhire_TX_catch_rate_lb*0.453592
    
    #For the states with both state-forhire and private angling (MS and TX), we set the state for-hire season as MS about 20 days (with 5% CV), TX 365 days.
    #The catch rate are all assumed as 5% CV. Fishing season set as 20% CV
    #reference: MS: state forhire VS private angling 17* 176 (2992) VS 81* 1823 (147663). TX: state forhire VS private angling 365*225 (82125) VS 62 * 2880 (178560)
    mean_Forhire_MS_season_length_temp<-log(Forhire_MS_season_length_mean)-log(1+imple_error^2)/2
    sd_Forhire_MS_season_length_temp<-sqrt(log(1+imple_error^2))
    Forhire_MS_season_length_temp<-rlnorm(1,mean_Forhire_MS_season_length_temp,sd_Forhire_MS_season_length_temp)
    while((Forhire_MS_season_length_temp>Forhire_MS_season_length_mean*2)|(Forhire_MS_season_length_temp<Forhire_MS_season_length_mean/2)){
      Forhire_MS_season_length_temp<-rlnorm(1,mean_Forhire_MS_season_length_temp,sd_Forhire_MS_season_length_temp)
    }
    Forhire_MS_season_length<-round(Forhire_MS_season_length_temp)
    
    mean_Plan_Private_AL_catch_rate_kg<-log(Private_AL_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Private_AL_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Private_AL_catch_rate_kg<-rlnorm(1,mean_Plan_Private_AL_catch_rate_kg,sd_Plan_Private_AL_catch_rate_kg)
    while((Plan_Private_AL_catch_rate_kg>Private_AL_catch_rate_kg*2)|(Plan_Private_AL_catch_rate_kg<Private_AL_catch_rate_kg/2)){
      Plan_Private_AL_catch_rate_kg<-rlnorm(1,mean_Plan_Private_AL_catch_rate_kg,sd_Plan_Private_AL_catch_rate_kg)
    }
    
    mean_Plan_Private_FL_catch_rate_kg<-log(Private_FL_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Private_FL_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Private_FL_catch_rate_kg<-rlnorm(1,mean_Plan_Private_FL_catch_rate_kg,sd_Plan_Private_FL_catch_rate_kg)
    while((Plan_Private_FL_catch_rate_kg>Private_FL_catch_rate_kg*2)|(Plan_Private_FL_catch_rate_kg<Private_FL_catch_rate_kg/2)){
      Plan_Private_FL_catch_rate_kg<-rlnorm(1,mean_Plan_Private_FL_catch_rate_kg,sd_Plan_Private_FL_catch_rate_kg)
    }
    
    mean_Plan_Private_LA_catch_rate_kg<-log(Private_LA_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Private_LA_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Private_LA_catch_rate_kg<-rlnorm(1,mean_Plan_Private_LA_catch_rate_kg,sd_Plan_Private_LA_catch_rate_kg)
    while((Plan_Private_LA_catch_rate_kg>Private_LA_catch_rate_kg*2)|(Plan_Private_LA_catch_rate_kg<Private_LA_catch_rate_kg/2)){
      Plan_Private_LA_catch_rate_kg<-rlnorm(1,mean_Plan_Private_LA_catch_rate_kg,sd_Plan_Private_LA_catch_rate_kg)
    }
    
    mean_Plan_Private_MS_catch_rate_kg<-log(Private_MS_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Private_MS_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Private_MS_catch_rate_kg<-rlnorm(1,mean_Plan_Private_MS_catch_rate_kg,sd_Plan_Private_MS_catch_rate_kg)
    while((Plan_Private_MS_catch_rate_kg>Private_MS_catch_rate_kg*2)|(Plan_Private_MS_catch_rate_kg<Private_MS_catch_rate_kg/2)){
      Plan_Private_MS_catch_rate_kg<-rlnorm(1,mean_Plan_Private_MS_catch_rate_kg,sd_Plan_Private_MS_catch_rate_kg)
    }
    
    mean_Plan_Private_TX_catch_rate_kg<-log(Private_TX_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Private_TX_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Private_TX_catch_rate_kg<-rlnorm(1,mean_Plan_Private_TX_catch_rate_kg,sd_Plan_Private_TX_catch_rate_kg)
    while((Plan_Private_TX_catch_rate_kg>Private_TX_catch_rate_kg*2)|(Plan_Private_TX_catch_rate_kg<Private_TX_catch_rate_kg/2)){
      Plan_Private_TX_catch_rate_kg<-rlnorm(1,mean_Plan_Private_TX_catch_rate_kg,sd_Plan_Private_TX_catch_rate_kg)
    }
    
    mean_Plan_Forhire_MS_catch_rate_kg<-log(Forhire_MS_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Forhire_MS_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Forhire_MS_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_MS_catch_rate_kg,sd_Plan_Forhire_MS_catch_rate_kg)
    while((Plan_Forhire_MS_catch_rate_kg>Forhire_MS_catch_rate_kg*2)|(Plan_Forhire_MS_catch_rate_kg<Forhire_MS_catch_rate_kg/2)){
      Plan_Forhire_MS_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_MS_catch_rate_kg,sd_Plan_Forhire_MS_catch_rate_kg)
    }
    
    mean_Plan_Forhire_TX_catch_rate_kg<-log(Forhire_TX_catch_rate_kg)-log(1+imple_error2^2)/2
    sd_Plan_Forhire_TX_catch_rate_kg<-sqrt(log(1+imple_error2^2))
    Plan_Forhire_TX_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_TX_catch_rate_kg,sd_Plan_Forhire_TX_catch_rate_kg)
    while((Plan_Forhire_TX_catch_rate_kg>Forhire_TX_catch_rate_kg*2)|(Plan_Forhire_TX_catch_rate_kg<Forhire_TX_catch_rate_kg/2)){
      Plan_Forhire_TX_catch_rate_kg<-rlnorm(1,mean_Plan_Forhire_TX_catch_rate_kg,sd_Plan_Forhire_TX_catch_rate_kg)
    }

    if(Private_season_switch==1) {
      #season determined by ACT
      mean_Plan_quote_private_AL<-log(Quote_private_aftbuf_AL)-log(1+imple_error^2)/2
      sd_Plan_quote_private_AL<-sqrt(log(1+imple_error^2))
      Plan_quote_private_AL<-rlnorm(1,mean_Plan_quote_private_AL,sd_Plan_quote_private_AL)
      while((Plan_quote_private_AL>Quote_private_aftbuf_AL*2)|(Plan_quote_private_AL<Quote_private_aftbuf_AL/2)){
        Plan_quote_private_AL<-rlnorm(1,mean_Plan_quote_private_AL,sd_Plan_quote_private_AL)
      }
      
      mean_Plan_quote_private_FL<-log(Quote_private_aftbuf_FL)-log(1+imple_error^2)/2
      sd_Plan_quote_private_FL<-sqrt(log(1+imple_error^2))
      Plan_quote_private_FL<-rlnorm(1,mean_Plan_quote_private_FL,sd_Plan_quote_private_FL)
      while((Plan_quote_private_FL>Quote_private_aftbuf_FL*2)|(Plan_quote_private_FL<Quote_private_aftbuf_FL/2)){
        Plan_quote_private_FL<-rlnorm(1,mean_Plan_quote_private_FL,sd_Plan_quote_private_FL)
      }
      
      mean_Plan_quote_private_LA<-log(Quote_private_aftbuf_LA)-log(1+imple_error^2)/2
      sd_Plan_quote_private_LA<-sqrt(log(1+imple_error^2))
      Plan_quote_private_LA<-rlnorm(1,mean_Plan_quote_private_LA,sd_Plan_quote_private_LA)
      while((Plan_quote_private_LA>Quote_private_aftbuf_LA*2)|(Plan_quote_private_LA<Quote_private_aftbuf_LA/2)){
        Plan_quote_private_LA<-rlnorm(1,mean_Plan_quote_private_LA,sd_Plan_quote_private_LA)
      }
      
      mean_Plan_quote_private_MS<-log(Quote_private_aftbuf_MS)-log(1+imple_error^2)/2
      sd_Plan_quote_private_MS<-sqrt(log(1+imple_error^2))
      Plan_quote_private_MS<-rlnorm(1,mean_Plan_quote_private_MS,sd_Plan_quote_private_MS)
      while((Plan_quote_private_MS>Quote_private_aftbuf_MS*2)|(Plan_quote_private_MS<Quote_private_aftbuf_MS/2)){
        Plan_quote_private_MS<-rlnorm(1,mean_Plan_quote_private_MS,sd_Plan_quote_private_MS)
      }
      
      mean_Plan_quote_private_TX<-log(Quote_private_aftbuf_TX)-log(1+imple_error^2)/2
      sd_Plan_quote_private_TX<-sqrt(log(1+imple_error^2))
      Plan_quote_private_TX<-rlnorm(1,mean_Plan_quote_private_TX,sd_Plan_quote_private_TX)
      while((Plan_quote_private_TX>Quote_private_aftbuf_TX*2)|(Plan_quote_private_TX<Quote_private_aftbuf_TX/2)){
        Plan_quote_private_TX<-rlnorm(1,mean_Plan_quote_private_TX,sd_Plan_quote_private_TX)
      }
      
    }else if(Private_season_switch==2) {#input by user
      Plan_quote_private_AL<-(Plan_Private_AL_catch_rate_kg/1000)*Private_AL_season_length
      Plan_quote_private_FL<-(Plan_Private_FL_catch_rate_kg/1000)*Private_FL_season_length
      Plan_quote_private_LA<-(Plan_Private_LA_catch_rate_kg/1000)*Private_LA_season_length
      Plan_quote_private_MS<-(Plan_Private_MS_catch_rate_kg/1000)*Private_MS_season_length + (Plan_Forhire_MS_catch_rate_kg/1000)*Forhire_MS_season_length
      Plan_quote_private_TX<-(Plan_Private_TX_catch_rate_kg/1000)*Private_TX_season_length + (Plan_Forhire_TX_catch_rate_kg/1000)*Forhire_TX_season_length
    }
    
    Plan_quote_comm<-Plan_quote_comm+comm_overunderage_lastyear
    Plan_quote_forhire<-Plan_quote_forhire+forhire_overunderage_lastyear
    Plan_quote_private_AL<-Plan_quote_private_AL+AL_overunderage_lastyear
    Plan_quote_private_FL<-Plan_quote_private_FL+FL_overunderage_lastyear
    Plan_quote_private_LA<-Plan_quote_private_LA+LA_overunderage_lastyear
    Plan_quote_private_MS<-Plan_quote_private_MS+MS_overunderage_lastyear
    Plan_quote_private_TX<-Plan_quote_private_TX+TX_overunderage_lastyear
  
    Comme_planned_catch<-Plan_quote_comm
    Recre_planned_catch<-Plan_quote_forhire + Plan_quote_private_AL + Plan_quote_private_FL + Plan_quote_private_LA + Plan_quote_private_MS + Plan_quote_private_TX
    
    mean_M_1_relative<-log(1)-log(1+cv_M_1^2)/2
    sd_M_1_relative<-sqrt(log(1+cv_M_1^2))
    M_1_relative<-rlnorm(1,mean_M_1_relative,sd_M_1_relative)
    while ((M_1_relative>2)|(M_1_relative<1/2)){
      M_1_relative<-rlnorm(1,mean_M_1_relative,sd_M_1_relative)
    }
    
    mean_M_2_relative<-log(1)-log(1+cv_M_2^2)/2
    sd_M_2_relative<-sqrt(log(1+cv_M_2^2))
    M_2_relative<-rlnorm(1,mean_M_2_relative,sd_M_2_relative)
    while ((M_2_relative>2)|(M_2_relative<1/2)){
      M_2_relative<-rlnorm(1,mean_M_2_relative,sd_M_2_relative)
    }
    M_1_true<-M_1_relative*M_1
    M_2_true<-M_2_relative*M_2
    
    mean_Comme_hl_discard_E_open_true<-log(Comme_hl_discard_E_open)-log(1+Comme_hl_discard_E_open_CV^2)/2
    sd_Comme_hl_discard_E_open_true<-sqrt(log(1+Comme_hl_discard_E_open_CV^2))
    Comme_hl_discard_E_open_true<- rlnorm(1,mean_Comme_hl_discard_E_open_true, sd_Comme_hl_discard_E_open_true)
    while ((Comme_hl_discard_E_open_true>2*Comme_hl_discard_E_open)|(Comme_hl_discard_E_open_true<1/2*Comme_hl_discard_E_open)){
      Comme_hl_discard_E_open_true<- rlnorm(1,mean_Comme_hl_discard_E_open_true, sd_Comme_hl_discard_E_open_true)
    }
    
    mean_Comme_ll_discard_E_open_true<- log(Comme_ll_discard_E_open)-log(1+Comme_ll_discard_E_open_CV^2)/2
    sd_Comme_ll_discard_E_open_true<- sqrt(log(1+Comme_ll_discard_E_open_CV^2))
    Comme_ll_discard_E_open_true<- rlnorm(1,mean_Comme_ll_discard_E_open_true, sd_Comme_ll_discard_E_open_true)
    while ((Comme_ll_discard_E_open_true>2*Comme_ll_discard_E_open)|(Comme_ll_discard_E_open_true<1/2*Comme_ll_discard_E_open)){
      Comme_ll_discard_E_open_true<- rlnorm(1,mean_Comme_ll_discard_E_open_true, sd_Comme_ll_discard_E_open_true)
    }
    
    mean_Comme_hl_discard_W_open_true<-log(Comme_hl_discard_W_open)-log(1+Comme_hl_discard_W_open_CV^2)/2
    sd_Comme_hl_discard_W_open_true<-sqrt(log(1+Comme_hl_discard_W_open_CV^2))
    Comme_hl_discard_W_open_true<- rlnorm(1,mean_Comme_hl_discard_W_open_true, sd_Comme_hl_discard_W_open_true)
    while((Comme_hl_discard_W_open_true>2*Comme_hl_discard_W_open)|(Comme_hl_discard_W_open_true<1/2*Comme_hl_discard_W_open)){
      Comme_hl_discard_W_open_true<- rlnorm(1,mean_Comme_hl_discard_W_open_true, sd_Comme_hl_discard_W_open_true)
    }
    
    mean_Comme_ll_discard_W_open_true<-log(Comme_ll_discard_W_open)-log(1+Comme_ll_discard_W_open_CV^2)/2
    sd_Comme_ll_discard_W_open_true<-sqrt(log(1+Comme_ll_discard_W_open_CV^2))
    Comme_ll_discard_W_open_true<- rlnorm(1,mean_Comme_ll_discard_W_open_true, sd_Comme_ll_discard_W_open_true)
    while((Comme_ll_discard_W_open_true>2*Comme_ll_discard_W_open)|(Comme_ll_discard_W_open_true<1/2*Comme_ll_discard_W_open)){
      Comme_ll_discard_W_open_true<- rlnorm(1,mean_Comme_ll_discard_W_open_true, sd_Comme_ll_discard_W_open_true)
    }
    
    mean_Recre_discard_E_open_true<-log(Recre_discard_E_open)-log(1+Recre_discard_E_open_CV^2)/2
    sd_Recre_discard_E_open_true<-sqrt(log(1+Recre_discard_E_open_CV^2))
    Recre_discard_E_open_true<-rlnorm(1,mean_Recre_discard_E_open_true, sd_Recre_discard_E_open_true)
    while((Recre_discard_E_open_true>2*Recre_discard_E_open)|(Recre_discard_E_open_true<1/2*Recre_discard_E_open)){
      Recre_discard_E_open_true<-rlnorm(1,mean_Recre_discard_E_open_true, sd_Recre_discard_E_open_true)
    }
    
    mean_Recre_discard_W_open_true<-log(Recre_discard_W_open)-log(1+Recre_discard_W_open_CV^2)/2
    sd_Recre_discard_W_open_true<-sqrt(log(1+Recre_discard_W_open_CV^2))
    Recre_discard_W_open_true<-rlnorm(1,mean_Recre_discard_W_open_true, sd_Recre_discard_W_open_true)
    while((Recre_discard_W_open_true>2*Recre_discard_W_open)|(Recre_discard_W_open_true<1/2*Recre_discard_W_open)){
      Recre_discard_W_open_true<-rlnorm(1,mean_Recre_discard_W_open_true, sd_Recre_discard_W_open_true)
    }
    
    mean_Recre_discard_E_closed_true<-log(Recre_discard_E_closed)-log(1+Recre_discard_E_closed_CV^2)/2
    sd_Recre_discard_E_closed_true<-sqrt(log(1+Recre_discard_E_closed_CV^2))
    Recre_discard_E_closed_true<-rlnorm(1,mean_Recre_discard_E_closed_true, sd_Recre_discard_E_closed_true)
    while((Recre_discard_E_closed_true>2*Recre_discard_E_closed)|(Recre_discard_E_closed_true<1/2*Recre_discard_E_closed)){
      Recre_discard_E_closed_true<-rlnorm(1,mean_Recre_discard_E_closed_true, sd_Recre_discard_E_closed_true)
    }
    
    mean_Recre_discard_W_closed_true<-log(Recre_discard_W_closed)-log(1+Recre_discard_W_closed_CV^2)/2
    sd_Recre_discard_W_closed_true<-sqrt(log(1+Recre_discard_W_closed_CV^2))
    Recre_discard_W_closed_true<-rlnorm(1,mean_Recre_discard_W_closed_true, sd_Recre_discard_W_closed_true)
    while((Recre_discard_W_closed_true>2*Recre_discard_W_closed)|(Recre_discard_W_closed_true<1/2*Recre_discard_W_closed)){
      Recre_discard_W_closed_true<-rlnorm(1,mean_Recre_discard_W_closed_true, sd_Recre_discard_W_closed_true)
    }
    
    mean_Comme_discard_E_closed_true<-log(Comme_discard_E_closed)-log(1+Comme_discard_E_closed_CV^2)/2
    sd_Comme_discard_E_closed_true<-sqrt(log(1+Comme_discard_E_closed_CV^2))
    Comme_discard_E_closed_true<-rlnorm(1,mean_Comme_discard_E_closed_true, sd_Comme_discard_E_closed_true)
    while((Comme_discard_E_closed_true>2*Comme_discard_E_closed)|(Comme_discard_E_closed_true<1/2*Comme_discard_E_closed)){
      Comme_discard_E_closed_true<-rlnorm(1,mean_Comme_discard_E_closed_true, sd_Comme_discard_E_closed_true)
    }
    
    mean_Comme_discard_W_closed_true<-log(Comme_discard_W_closed)-log(1+Comme_discard_W_closed_CV^2)/2
    sd_Comme_discard_W_closed_true<-sqrt(log(1+Comme_discard_W_closed_CV^2))
    Comme_discard_W_closed_true<-rlnorm(1,mean_Comme_discard_W_closed_true, sd_Comme_discard_W_closed_true)
    while((Comme_discard_W_closed_true>2*Comme_discard_W_closed)|(Comme_discard_W_closed_true<1/2*Comme_discard_W_closed)){
      Comme_discard_W_closed_true<-rlnorm(1,mean_Comme_discard_W_closed_true, sd_Comme_discard_W_closed_true)
    }
    
    comp_comm_discard_e_open_true <- hl_e_pred_F_ave * hl_e_selex * (1-act_hl_e_retention) + ll_e_pred_F_ave * ll_e_selex * (1-act_ll_e_retention)
    comp_comm_discard_w_open_true <- hl_w_pred_F_ave * hl_w_selex * (1-act_hl_w_retention) + ll_w_pred_F_ave * ll_w_selex * (1-act_ll_w_retention)
    
    comp_comm_discard_e_death_open_true <- hl_e_pred_F_ave * hl_e_selex * (1-act_hl_e_retention) * Comme_hl_discard_E_open_true + ll_e_pred_F_ave * ll_e_selex * (1-act_ll_e_retention) * Comme_ll_discard_E_open_true
    comp_comm_discard_w_death_open_true <- hl_w_pred_F_ave * hl_w_selex * (1-act_hl_w_retention) * Comme_hl_discard_W_open_true + ll_w_pred_F_ave * ll_w_selex * (1-act_ll_w_retention) * Comme_ll_discard_W_open_true
    
    comp_recr_discard_e_open_true <- mrip_e_pred_F_ave * mrip_e_selex * (1-act_mrip_e_retention) + hbt_e_pred_F_ave* hbt_e_selex * (1-act_hbt_e_retention)
    comp_recr_discard_w_open_true <- mrip_w_pred_F_ave * mrip_w_selex * (1-act_mrip_w_retention) + hbt_w_pred_F_ave* hbt_w_selex * (1-act_hbt_w_retention)
    
    comp_recr_discard_e_death_open_true <- mrip_e_pred_F_ave * mrip_e_selex * (1-act_mrip_e_retention) * Recre_discard_E_open_true + hbt_e_pred_F_ave* hbt_e_selex * (1-act_hbt_e_retention) * Recre_discard_E_open_true
    comp_recr_discard_w_death_open_true <- mrip_w_pred_F_ave * mrip_w_selex * (1-act_mrip_w_retention) * Recre_discard_W_open_true + hbt_w_pred_F_ave* hbt_w_selex * (1-act_hbt_w_retention) * Recre_discard_W_open_true
    
    comp_discard_e_close_true <- rec_closed_e_pred_F_ave * rec_closed_e_selex + comm_closed_e_pred_F_ave * comm_closed_e_selex
    comp_discard_w_close_true <- rec_closed_w_pred_F_ave * rec_closed_w_selex + comm_closed_w_pred_F_ave * comm_closed_w_selex
    
    comp_discard_e_death_closed_true <- rec_closed_e_pred_F_ave * rec_closed_e_selex * Recre_discard_E_closed_true + comm_closed_e_pred_F_ave * comm_closed_e_selex * Comme_discard_E_closed_true
    comp_discard_w_death_closed_true <- rec_closed_w_pred_F_ave * rec_closed_w_selex * Recre_discard_W_closed_true + comm_closed_w_pred_F_ave * comm_closed_w_selex * Comme_discard_W_closed_true
    
    comp_bycatch_e_true <- shrimp_e_pred_F_ave * shrimp_e_selex
    comp_bycatch_w_true <- shrimp_w_pred_F_ave * shrimp_w_selex
    
    relF_comm_imple<-optimize(Function_comm_imple,lower=0,upper=1000,tol=0.01)$minimum
    relF_recr_imple<-optimize(Function_recr_imple,lower=0,upper=1000,tol=0.01)$minimum
    #relF_comm_imple<-optim(c(1, 1),Function_comm_imple,method="Nelder-Mead",hessian=TRUE)$par
    #relF_recr_imple<-optim(c(1, 1),Function_recr_imple,method="Nelder-Mead",hessian=TRUE)$par
    #relF_imple<-optim(c(1, 1),Function_relF_imple,method="L-BFGS-B",lower = c(0,0), upper = c(10,10),hessian=TRUE)$par
    
    comm_relF_imple<-max(relF_comm_imple,0)
    recr_relF_imple<-max(relF_recr_imple,0)
    
    total_F_1_imple<-recr_relF_imple * comp_recr_sel_e + comm_relF_imple * comp_comm_sel_e + 
      recr_relF_imple * comp_recr_discard_e_death_open_true + comm_relF_imple * comp_comm_discard_e_death_open_true + comp_discard_e_death_closed_true + comp_bycatch_e_true
    
    total_F_2_imple<-recr_relF_imple * comp_recr_sel_w + comm_relF_imple * comp_comm_sel_w + 
      recr_relF_imple * comp_recr_discard_w_death_open_true + comm_relF_imple * comp_comm_discard_w_death_open_true + comp_discard_w_death_closed_true + comp_bycatch_w_true 
    
    mod_F_1_imple<-recr_relF_imple * comp_recr_sel_e + comm_relF_imple * comp_comm_sel_e
    mod_F_2_imple<-recr_relF_imple * comp_recr_sel_w + comm_relF_imple * comp_comm_sel_w
    
    true_N_1_nextyeartemp<-true_N_1*exp(-total_F_1_imple -M_1_true)
    true_N_2_nextyeartemp<-true_N_2*exp(-total_F_2_imple -M_2_true)
    
    catch_N_1_imple<-(true_N_1-true_N_1_nextyeartemp)*mod_F_1_imple/(total_F_1_imple + M_1_true)
    catch_N_2_imple<-(true_N_2-true_N_2_nextyeartemp)*mod_F_2_imple/(total_F_2_imple + M_2_true)
    
    total_SSB_paratemp[i.runtime]<-sum(SSB_1_paratemp[i.runtime]+SSB_2_paratemp[i.runtime])
    totalcatch_1_imple_paratemp[i.runtime]<-sum(catch_N_1_imple * weight_at_age_1) #unit mt
    totalcatch_2_imple_paratemp[i.runtime]<-sum(catch_N_2_imple * weight_at_age_2) #unit mt
    total_catch_paratemp[i.runtime]<-totalcatch_1_imple_paratemp[i.runtime]+totalcatch_2_imple_paratemp[i.runtime]
    comm_catch_1_paratemp[i.runtime]<-sum((true_N_1-true_N_1_nextyeartemp)*comm_relF_imple * comp_comm_sel_e/(total_F_1_imple + M_1_true)*weight_at_age_1)
    comm_catch_2_paratemp[i.runtime]<-sum((true_N_2-true_N_2_nextyeartemp)*comm_relF_imple * comp_comm_sel_w/(total_F_2_imple + M_2_true)*weight_at_age_2)
    comm_catch_paratemp[i.runtime]<-comm_catch_1_paratemp[i.runtime]+comm_catch_2_paratemp[i.runtime]
    recr_catch_paratemp[i.runtime]<-sum((true_N_1-true_N_1_nextyeartemp)*recr_relF_imple * comp_recr_sel_e/(total_F_1_imple + M_1_true)*weight_at_age_1+
                                    (true_N_2-true_N_2_nextyeartemp)*recr_relF_imple * comp_recr_sel_w/(total_F_2_imple + M_2_true)*weight_at_age_2)
    Forhire_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_forhire/(Recre_planned_catch+0.000001)
    Forhire_catch_1_paratemp[i.runtime]<-Forhire_catch_paratemp[i.runtime]*Ratio_forhire_east
    Forhire_catch_2_paratemp[i.runtime]<-Forhire_catch_paratemp[i.runtime]*Ratio_forhire_west
      
    Private_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*(Recre_planned_catch-Plan_quote_forhire)/(Recre_planned_catch+0.000001)  
    Private_catch_1_paratemp[i.runtime]<-Private_catch_paratemp[i.runtime]*Ratio_private_east
    Private_catch_2_paratemp[i.runtime]<-Private_catch_paratemp[i.runtime]*Ratio_private_west
    
    comm_discards_paratemp[i.runtime]<-sum((true_N_1-true_N_1_nextyeartemp)*(comm_relF_imple * comp_comm_discard_e_death_open_true)/(total_F_1_imple + M_1_true)*weight_at_age_1+
                                        (true_N_2-true_N_2_nextyeartemp)*(comm_relF_imple * comp_comm_discard_w_death_open_true)/(total_F_2_imple + M_2_true)*weight_at_age_2)
    recr_discards_paratemp[i.runtime]<-sum((true_N_1-true_N_1_nextyeartemp)*(recr_relF_imple * comp_recr_discard_e_death_open_true)/(total_F_1_imple + M_1_true)*weight_at_age_1+
                                        (true_N_2-true_N_2_nextyeartemp)*(recr_relF_imple * comp_recr_discard_e_death_open_true)/(total_F_2_imple + M_2_true)*weight_at_age_2)
    comm_dislandratio_paratemp[i.runtime]<-comm_discards_paratemp[i.runtime]/(comm_catch_paratemp[i.runtime]+0.000001)
    recr_dislandratio_paratemp[i.runtime]<-recr_discards_paratemp[i.runtime]/(recr_catch_paratemp[i.runtime]+0.000001)
 
    true_private_AL_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_private_AL/(Recre_planned_catch+0.000001)
    true_private_FL_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_private_FL/(Recre_planned_catch+0.000001)
    true_private_LA_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_private_LA/(Recre_planned_catch+0.000001)
    true_private_MS_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_private_MS/(Recre_planned_catch+0.000001)
    true_private_TX_catch_paratemp[i.runtime]<-recr_catch_paratemp[i.runtime]*Plan_quote_private_TX/(Recre_planned_catch+0.000001)
    true_fed_forhire_season_length_paratemp[i.runtime]<-round(Forhire_catch_paratemp[i.runtime]/(Plan_Forhire_federal_catch_rate_kg/1000))
    true_private_AL_season_length_paratemp[i.runtime]<-round(true_private_AL_catch_paratemp[i.runtime]/(Plan_Private_AL_catch_rate_kg/1000))
    true_private_FL_season_length_paratemp[i.runtime]<-round(true_private_FL_catch_paratemp[i.runtime]/(Plan_Private_FL_catch_rate_kg/1000))
    true_private_LA_season_length_paratemp[i.runtime]<-round(true_private_LA_catch_paratemp[i.runtime]/(Plan_Private_LA_catch_rate_kg/1000))
    true_private_MS_season_length_paratemp[i.runtime]<-max(0,round((true_private_MS_catch_paratemp[i.runtime]-Plan_Forhire_MS_catch_rate_kg/1000*Forhire_MS_season_length)/(Plan_Private_MS_catch_rate_kg/1000)))
    true_private_TX_season_length_paratemp[i.runtime]<-max(0,round((true_private_TX_catch_paratemp[i.runtime]-Plan_Forhire_TX_catch_rate_kg/1000*Forhire_TX_season_length)/(Plan_Private_TX_catch_rate_kg/1000)))
    
    #The next 9 lines seems redundant, but may not hurt
    Quote_comm_aftbuf<-Quote_comm_befbuf*(1-Comme_buffer)
    Quote_forhire_aftbuf<-Quote_forhire_befbuf*(1-Forhire_buffer)
    Quote_private_aftbuf<-Quote_private_befbuf*(1-Private_buffer)
    
    Quote_private_aftbuf_AL<-Private_quota_AL*Quote_private_aftbuf
    Quote_private_aftbuf_FL<-Private_quota_FL*Quote_private_aftbuf
    Quote_private_aftbuf_LA<-Private_quota_LA*Quote_private_aftbuf
    Quote_private_aftbuf_MS<-Private_quota_MS*Quote_private_aftbuf
    Quote_private_aftbuf_TX<-Private_quota_TX*Quote_private_aftbuf
    #estimate AM overunderage for the next year
    #Quote_forhire_aftbuf+Quote_private_aftbuf-recr_catch_paratemp[i.runtime]
    AM_comm_overunderage_paratemp[i.runtime+1]<-Quote_comm_aftbuf-comm_catch_paratemp[i.runtime]
    AM_forhire_overunderage_paratemp[i.runtime+1]<-Quote_forhire_aftbuf-Forhire_catch_paratemp[i.runtime]
    AM_AL_overunderage_paratemp[i.runtime+1]<-Quote_private_aftbuf_AL-true_private_AL_catch_paratemp[i.runtime]
    AM_FL_overunderage_paratemp[i.runtime+1]<-Quote_private_aftbuf_FL-true_private_FL_catch_paratemp[i.runtime]
    AM_LA_overunderage_paratemp[i.runtime+1]<-Quote_private_aftbuf_LA-true_private_LA_catch_paratemp[i.runtime]
    AM_MS_overunderage_paratemp[i.runtime+1]<-Quote_private_aftbuf_MS-true_private_MS_catch_paratemp[i.runtime]
    AM_TX_overunderage_paratemp[i.runtime+1]<-Quote_private_aftbuf_TX-true_private_TX_catch_paratemp[i.runtime]
    
    overunderage<-ACL_planned-total_catch_paratemp[i.runtime] #overage for the next year
    
    F_general_paratemp[i.runtime+2]<-sum((true_N_1-true_N_1_nextyeartemp)*total_F_1_imple/(total_F_1_imple + M_1_true)+(true_N_2-true_N_2_nextyeartemp)*total_F_2_imple/(total_F_2_imple + M_2_true))/sum(true_N_1+true_N_2)
    
    true_N_1[2:(length(true_N_1)-1)]<-true_N_1_nextyeartemp[1:(length(true_N_1)-2)]
    true_N_2[2:(length(true_N_2)-1)]<-true_N_2_nextyeartemp[1:(length(true_N_2)-2)]
    true_N_1[length(true_N_1)]<-true_N_1_nextyeartemp[length(true_N_1)-1]+true_N_1_nextyeartemp[length(true_N_1)]
    true_N_2[length(true_N_2)]<-true_N_2_nextyeartemp[length(true_N_2)-1]+true_N_2_nextyeartemp[length(true_N_2)]
    
    true_N_1[1]<-true_R_1
    true_N_2[1]<-true_R_2
    
    F_ratio_paratemp[i.runtime]<-(F_general_paratemp[i.runtime]+F_general_paratemp[i.runtime+1]+F_general_paratemp[i.runtime+2])/3/MFMT
    SSB_ratio_paratemp[i.runtime]<-(SSB_1_paratemp[i.runtime]+SSB_2_paratemp[i.runtime])/SSB_MSY_BRP

    if(SSB_ratio_paratemp[i.runtime]>=1){
      if(F_ratio_paratemp[i.runtime]<=1){
        Year_green_paratemp[i.runtime]<-1 #1: green 
      }else{
        Year_green_paratemp[i.runtime]<--2 #2: grey green -2
      }
    }else{
      if(F_ratio_paratemp[i.runtime]<=1){
        Year_green_paratemp[i.runtime]<--3 #1: orange -3
      }else{
        Year_green_paratemp[i.runtime]<--4 #2: red -4
      }
    }
      
  }
  
  #save parallel temp
  write.csv(R_1_paratemp,paste0(i.run,"_R_1_paratemp.csv"))
  write.csv(R_2_paratemp,paste0(i.run,"_R_2_paratemp.csv"))
  write.csv(SSB_1_paratemp,paste0(i.run,"_SSB_1_paratemp.csv"))
  write.csv(SSB_2_paratemp,paste0(i.run,"_SSB_2_paratemp.csv"))
  write.csv(total_SSB_paratemp,paste0(i.run,"_total_SSB_paratemp.csv"))
  write.csv(totalcatch_1_imple_paratemp,paste0(i.run,"_totalcatch_1_imple_paratemp.csv"))
  write.csv(totalcatch_2_imple_paratemp,paste0(i.run,"_totalcatch_2_imple_paratemp.csv"))
  write.csv(total_catch_paratemp,paste0(i.run,"_total_catch_paratemp.csv"))
  write.csv(comm_catch_paratemp,paste0(i.run,"_comm_catch_paratemp.csv"))
  write.csv(comm_catch_1_paratemp,paste0(i.run,"_comm_catch_1_paratemp.csv"))
  write.csv(comm_catch_2_paratemp,paste0(i.run,"_comm_catch_2_paratemp.csv"))
  write.csv(recr_catch_paratemp,paste0(i.run,"_recr_catch_paratemp.csv"))
  write.csv(Forhire_catch_paratemp,paste0(i.run,"_Forhire_catch_paratemp.csv"))
  write.csv(Forhire_catch_1_paratemp,paste0(i.run,"_Forhire_catch_1_paratemp.csv"))
  write.csv(Forhire_catch_2_paratemp,paste0(i.run,"_Forhire_catch_2_paratemp.csv"))
  write.csv(Private_catch_paratemp,paste0(i.run,"_Private_catch_paratemp.csv"))
  write.csv(Private_catch_1_paratemp,paste0(i.run,"_Private_catch_1_paratemp.csv"))
  write.csv(Private_catch_2_paratemp,paste0(i.run,"_Private_catch_2_paratemp.csv"))
  write.csv(comm_discards_paratemp,paste0(i.run,"_comm_discards_paratemp.csv"))
  write.csv(recr_discards_paratemp,paste0(i.run,"_recr_discards_paratemp.csv"))
  write.csv(comm_dislandratio_paratemp,paste0(i.run,"_comm_dislandratio_paratemp.csv"))
  write.csv(recr_dislandratio_paratemp,paste0(i.run,"_recr_dislandratio_paratemp.csv"))
  write.csv(true_fed_forhire_season_length_paratemp,paste0(i.run,"_true_fed_forhire_season_length_paratemp.csv"))
  write.csv(true_private_AL_season_length_paratemp,paste0(i.run,"_true_private_AL_season_length_paratemp.csv"))
  write.csv(true_private_FL_season_length_paratemp,paste0(i.run,"_true_private_FL_season_length_paratemp.csv"))
  write.csv(true_private_LA_season_length_paratemp,paste0(i.run,"_true_private_LA_season_length_paratemp.csv"))
  write.csv(true_private_MS_season_length_paratemp,paste0(i.run,"_true_private_MS_season_length_paratemp.csv"))
  write.csv(true_private_TX_season_length_paratemp,paste0(i.run,"_true_private_TX_season_length_paratemp.csv"))
  write.csv(true_private_AL_catch_paratemp,paste0(i.run,"_true_private_AL_catch_paratemp.csv"))
  write.csv(true_private_FL_catch_paratemp,paste0(i.run,"_true_private_FL_catch_paratemp.csv"))
  write.csv(true_private_LA_catch_paratemp,paste0(i.run,"_true_private_LA_catch_paratemp.csv"))
  write.csv(true_private_MS_catch_paratemp,paste0(i.run,"_true_private_MS_catch_paratemp.csv"))
  write.csv(true_private_TX_catch_paratemp,paste0(i.run,"_true_private_TX_catch_paratemp.csv"))
  write.csv(F_general_paratemp,paste0(i.run,"_F_general_paratemp.csv"))
  write.csv(AM_comm_overunderage_paratemp,paste0(i.run,"_AM_comm_overunderage_paratemp.csv"))
  write.csv(AM_forhire_overunderage_paratemp,paste0(i.run,"_AM_forhire_overunderage_paratemp.csv"))
  write.csv(AM_AL_overunderage_paratemp,paste0(i.run,"_AM_AL_overunderage_paratemp.csv"))
  write.csv(AM_FL_overunderage_paratemp,paste0(i.run,"_AM_FL_overunderage_paratemp.csv"))
  write.csv(AM_LA_overunderage_paratemp,paste0(i.run,"_AM_LA_overunderage_paratemp.csv"))
  write.csv(AM_MS_overunderage_paratemp,paste0(i.run,"_AM_MS_overunderage_paratemp.csv"))
  write.csv(AM_TX_overunderage_paratemp,paste0(i.run,"_AM_TX_overunderage_paratemp.csv"))
  write.csv(Year_green_paratemp,paste0(i.run,"_Year_green_paratemp.csv"))
  write.csv(F_ratio_paratemp,paste0(i.run,"_F_ratio_paratemp.csv"))
  write.csv(SSB_ratio_paratemp,paste0(i.run,"_SSB_ratio_paratemp.csv"))
}

# end parallel job
stopCluster(cl) #end the cluster for parallel processing
Sys.time()

Function_readinpara<-function(filename){
  filedata<-NULL
  for(i.para in 1:Simrun_Num)
    filedata<-cbind(filedata,as.matrix(read.csv(paste0(i.para,"_",filename,"_paratemp.csv"),header=T))[,2])
  t(filedata)
}

Function_deletepara<-function(filename){
  for(i.para in 1:Simrun_Num)
    file.remove(paste0(i.para,"_",filename,"_paratemp.csv"))
}

R_1<-as.matrix(Function_readinpara("R_1"))
R_2<-as.matrix(Function_readinpara("R_2"))
SSB_1<-as.matrix(Function_readinpara("SSB_1"))
SSB_2<-as.matrix(Function_readinpara("SSB_2"))
total_SSB<-as.matrix(Function_readinpara("total_SSB"))
totalcatch_1_imple<-as.matrix(Function_readinpara("totalcatch_1_imple"))
totalcatch_2_imple<-as.matrix(Function_readinpara("totalcatch_2_imple"))
total_catch<-as.matrix(Function_readinpara("total_catch"))
comm_catch<-as.matrix(Function_readinpara("comm_catch"))
comm_catch_1<-as.matrix(Function_readinpara("comm_catch_1"))
comm_catch_2<-as.matrix(Function_readinpara("comm_catch_2"))
recr_catch<-as.matrix(Function_readinpara("recr_catch"))
Forhire_catch<-as.matrix(Function_readinpara("Forhire_catch"))
Forhire_catch_1<-as.matrix(Function_readinpara("Forhire_catch_1"))
Forhire_catch_2<-as.matrix(Function_readinpara("Forhire_catch_2"))
Private_catch<-as.matrix(Function_readinpara("Private_catch"))
Private_catch_1<-as.matrix(Function_readinpara("Private_catch_1"))
Private_catch_2<-as.matrix(Function_readinpara("Private_catch_2"))
comm_discards<-as.matrix(Function_readinpara("comm_discards"))
recr_discards<-as.matrix(Function_readinpara("recr_discards"))
comm_dislandratio<-as.matrix(Function_readinpara("comm_dislandratio"))
recr_dislandratio<-as.matrix(Function_readinpara("recr_dislandratio"))
true_fed_forhire_season_length<-as.matrix(Function_readinpara("true_fed_forhire_season_length"))
true_private_AL_season_length<-as.matrix(Function_readinpara("true_private_AL_season_length"))
true_private_FL_season_length<-as.matrix(Function_readinpara("true_private_FL_season_length"))
true_private_LA_season_length<-as.matrix(Function_readinpara("true_private_LA_season_length"))
true_private_MS_season_length<-as.matrix(Function_readinpara("true_private_MS_season_length"))
true_private_TX_season_length<-as.matrix(Function_readinpara("true_private_TX_season_length"))
true_private_AL_catch<-as.matrix(Function_readinpara("true_private_AL_catch"))
true_private_FL_catch<-as.matrix(Function_readinpara("true_private_FL_catch"))
true_private_LA_catch<-as.matrix(Function_readinpara("true_private_LA_catch"))
true_private_MS_catch<-as.matrix(Function_readinpara("true_private_MS_catch"))
true_private_TX_catch<-as.matrix(Function_readinpara("true_private_TX_catch"))
F_general<-as.matrix(Function_readinpara("F_general"))
AM_comm_overunderage<-as.matrix(Function_readinpara("AM_comm_overunderage"))
AM_forhire_overunderage<-as.matrix(Function_readinpara("AM_forhire_overunderage"))
AM_AL_overunderage<-as.matrix(Function_readinpara("AM_AL_overunderage"))
AM_FL_overunderage<-as.matrix(Function_readinpara("AM_FL_overunderage"))
AM_LA_overunderage<-as.matrix(Function_readinpara("AM_LA_overunderage"))
AM_MS_overunderage<-as.matrix(Function_readinpara("AM_MS_overunderage"))
AM_TX_overunderage<-as.matrix(Function_readinpara("AM_TX_overunderage"))
Year_green<-as.matrix(Function_readinpara("Year_green"))
F_ratio<-as.matrix(Function_readinpara("F_ratio"))
SSB_ratio<-as.matrix(Function_readinpara("SSB_ratio"))

Function_deletepara("R_1")
Function_deletepara("R_2")
Function_deletepara("SSB_1")
Function_deletepara("SSB_2")
Function_deletepara("total_SSB")
Function_deletepara("totalcatch_1_imple")
Function_deletepara("totalcatch_2_imple")
Function_deletepara("total_catch")
Function_deletepara("comm_catch")
Function_deletepara("comm_catch_1")
Function_deletepara("comm_catch_2")
Function_deletepara("recr_catch")
Function_deletepara("Forhire_catch")
Function_deletepara("Forhire_catch_1")
Function_deletepara("Forhire_catch_2")
Function_deletepara("Private_catch")
Function_deletepara("Private_catch_1")
Function_deletepara("Private_catch_2")
Function_deletepara("comm_discards")
Function_deletepara("recr_discards")
Function_deletepara("comm_dislandratio")
Function_deletepara("recr_dislandratio")
Function_deletepara("true_fed_forhire_season_length")
Function_deletepara("true_private_AL_season_length")
Function_deletepara("true_private_FL_season_length")
Function_deletepara("true_private_LA_season_length")
Function_deletepara("true_private_MS_season_length")
Function_deletepara("true_private_TX_season_length")
Function_deletepara("true_private_AL_catch")
Function_deletepara("true_private_FL_catch")
Function_deletepara("true_private_LA_catch")
Function_deletepara("true_private_MS_catch")
Function_deletepara("true_private_TX_catch")
Function_deletepara("F_general")
Function_deletepara("AM_comm_overunderage")
Function_deletepara("AM_forhire_overunderage")
Function_deletepara("AM_AL_overunderage")
Function_deletepara("AM_FL_overunderage")
Function_deletepara("AM_LA_overunderage")
Function_deletepara("AM_MS_overunderage")
Function_deletepara("AM_TX_overunderage")
Function_deletepara("Year_green")
Function_deletepara("F_ratio")
Function_deletepara("SSB_ratio")

Sys.time()

#save the following matrix for figures:
#Essential Figure
proj_year<-project_start_year:(project_start_year+Runtime_long-1)

#Figure 1 Line Gardient of total catch and total SSB, use median of 100 runs
total_catch_median<-colMedians(total_catch)*2.2046
total_SSB_median<-colMedians(total_SSB)

#plot(proj_year,total_catch_median, ylim=c(0,max(total_catch_median)*1.2))
#plot(proj_year,total_SSB_median, ylim=c(0,max(total_SSB_median)*1.2))

#Figure 2 Starcked Area Chart, use median of 100 runs
#total_catch_median<-colMedians(total_catch)
comm_catch_median<-colMedians(comm_catch)*2.2046
Forhire_catch_median<-colMedians(Forhire_catch)*2.2046
Private_catch_median<-colMedians(Private_catch)*2.2046

#plot(proj_year,comm_catch_median, ylim=c(0,max(comm_catch_median)*1.2))
#plot(proj_year,Forhire_catch_median, ylim=c(0,max(Forhire_catch_median)*1.2))
#plot(proj_year,Private_catch_median, ylim=c(0,max(Private_catch_median)*1.2))

#Figure 3 Starcked Area Chart, use median of 100 runs
#total_SSB_median<-colMedians(total_SSB)
SSB_1_median<-colMedians(SSB_1[,1:Runtime_long])
SSB_2_median<-colMedians(SSB_2[,1:Runtime_long])

#plot(proj_year,SSB_1_median, ylim=c(0,max(SSB_1_median)*1.2))
#plot(proj_year,SSB_2_median, ylim=c(0,max(SSB_2_median)*1.2))

#Figure 4 Starcked Area Chart, use median of 100 runs
#comm_catch_median<-colMedians(comm_catch)
comm_catch_1_median<-colMedians(comm_catch_1)*2.2046
comm_catch_2_median<-colMedians(comm_catch_2)*2.2046

#plot(proj_year,comm_catch_1_median, ylim=c(0,max(comm_catch_1_median)*1.2))
#plot(proj_year,comm_catch_2_median, ylim=c(0,max(comm_catch_2_median)*1.2))

#Figure 5 Starcked Area Chart, use median of 100 runs
#Forhire_catch_median<-colMedians(Forhire_catch)*2.2046
Forhire_catch_1_median<-colMedians(Forhire_catch_1)*2.2046
Forhire_catch_2_median<-colMedians(Forhire_catch_2)*2.2046
#plot(proj_year,Forhire_catch_1_median, ylim=c(0,max(Forhire_catch_1)*1.2))
#plot(proj_year,Forhire_catch_2_median, ylim=c(0,max(Forhire_catch_2)*1.2))

#Figure 6 Starcked Area Chart, use median of 100 runs
#Private_catch_median<-colMedians(Private_catch)*2.2046
Private_catch_1_median<-colMedians(Private_catch_1)*2.2046
Private_catch_2_median<-colMedians(Private_catch_2)*2.2046
#plot(proj_year,Private_catch_1_median, ylim=c(0,max(Private_catch_1_median)*1.2))
#plot(proj_year,Private_catch_2_median, ylim=c(0,max(Private_catch_2_median)*1.2))

#Figure 7 Bar Simple with variation, use median, and 95% CI of 100 runs
true_fed_forhire_season_length_median<-colMedians(true_fed_forhire_season_length)
true_fed_forhire_season_length_975<-apply(true_fed_forhire_season_length, 2, quantile, probs=0.975)
true_fed_forhire_season_length_025<-apply(true_fed_forhire_season_length, 2, quantile, probs=0.025)

#Figure 8 Share Dataset without the pie plot, use median of 100 runs
true_private_AL_season_length_median<-colMedians(true_private_AL_season_length)
true_private_FL_season_length_median<-colMedians(true_private_FL_season_length)  
true_private_LA_season_length_median<-colMedians(true_private_LA_season_length)  
true_private_MS_season_length_median<-colMedians(true_private_MS_season_length)  
true_private_TX_season_length_median<-colMedians(true_private_TX_season_length)  

#Figure 9 Kobe plot, use median of 100 runs
SSB_ratio_median<-colMedians(SSB_ratio)
F_ratio_median<-colMedians(F_ratio)
#plot(SSB_ratio_median,F_ratio_median,type="l")

#Figure 10 Starcked Area Chart, use median of 100 runs
R_1_median<-colMedians(R_1[,1:Runtime_long])
R_2_median<-colMedians(R_2[,1:Runtime_long])

#Other figures
#Figure 1 Confidence Band
comm_catch_median<-colMedians(comm_catch)*2.2046
comm_catch_975<-apply(comm_catch, 2, quantile, probs=0.975)*2.2046
comm_catch_025<-apply(comm_catch, 2, quantile, probs=0.025)*2.2046

#Figure 2 Confidence Band
recr_catch_median<-colMedians(recr_catch)*2.2046
recr_catch_975<-apply(recr_catch, 2, quantile, probs=0.975)*2.2046
recr_catch_025<-apply(recr_catch, 2, quantile, probs=0.025)*2.2046

#Figure 3 Confidence Band
Forhire_catch_median<-colMedians(Forhire_catch)*2.2046
Forhire_catch_975<-apply(Forhire_catch, 2, quantile, probs=0.975)*2.2046
Forhire_catch_025<-apply(Forhire_catch, 2, quantile, probs=0.025)*2.2046

#Figure 4 Confidence Band
Private_catch_median<-colMedians(Private_catch)*2.2046
Private_catch_975<-apply(Private_catch, 2, quantile, probs=0.975)*2.2046
Private_catch_025<-apply(Private_catch, 2, quantile, probs=0.025)*2.2046

#Figure 5 Confidence Band
F_general_median<-colMedians(F_general)
F_general_975<-apply(F_general, 2, quantile, probs=0.975)
F_general_025<-apply(F_general, 2, quantile, probs=0.025)

#Figure 6 Confidence Band
SSB_1_median<-colMedians(SSB_1[,1:Runtime_long])
SSB_1_975<-apply(SSB_1[,1:Runtime_long], 2, quantile, probs=0.975)
SSB_1_025<-apply(SSB_1[,1:Runtime_long], 2, quantile, probs=0.025)

#Figure 7 Confidence Band
SSB_2_median<-colMedians(SSB_2[,1:Runtime_long])
SSB_2_975<-apply(SSB_2[,1:Runtime_long], 2, quantile, probs=0.975)
SSB_2_025<-apply(SSB_2[,1:Runtime_long], 2, quantile, probs=0.025)

#Figure 8 Confidence Band
total_SSB_median<-colMedians(total_SSB)
total_SSB_975<-apply(total_SSB, 2, quantile, probs=0.975)
total_SSB_025<-apply(total_SSB, 2, quantile, probs=0.025)

#Figure 9 Confidence Band
true_private_AL_catch_median<-colMedians(true_private_AL_catch)*2.2046
true_private_AL_catch_975<-apply(true_private_AL_catch, 2, quantile, probs=0.975)*2.2046
true_private_AL_catch_025<-apply(true_private_AL_catch, 2, quantile, probs=0.025)*2.2046

#Figure 10 Confidence Band
true_private_AL_season_length_median<-colMedians(true_private_AL_season_length)
true_private_AL_season_length_975<-apply(true_private_AL_season_length, 2, quantile, probs=0.975)
true_private_AL_season_length_025<-apply(true_private_AL_season_length, 2, quantile, probs=0.025)

#Figure 11 Confidence Band
true_private_FL_catch_median<-colMedians(true_private_FL_catch)*2.2046
true_private_FL_catch_975<-apply(true_private_FL_catch, 2, quantile, probs=0.975)*2.2046
true_private_FL_catch_025<-apply(true_private_FL_catch, 2, quantile, probs=0.025)*2.2046

#Figure 12 Confidence Band
true_private_FL_season_length_median<-colMedians(true_private_FL_season_length)
true_private_FL_season_length_975<-apply(true_private_FL_season_length, 2, quantile, probs=0.975)
true_private_FL_season_length_025<-apply(true_private_FL_season_length, 2, quantile, probs=0.025)

#Figure 13 Confidence Band
true_private_LA_catch_median<-colMedians(true_private_LA_catch)*2.2046
true_private_LA_catch_975<-apply(true_private_LA_catch, 2, quantile, probs=0.975)*2.2046
true_private_LA_catch_025<-apply(true_private_LA_catch, 2, quantile, probs=0.025)*2.2046

#Figure 14 Confidence Band
true_private_LA_season_length_median<-colMedians(true_private_LA_season_length)
true_private_LA_season_length_975<-apply(true_private_LA_season_length, 2, quantile, probs=0.975)
true_private_LA_season_length_025<-apply(true_private_LA_season_length, 2, quantile, probs=0.025)

#Figure 15 Confidence Band
true_private_MS_catch_median<-colMedians(true_private_MS_catch)*2.2046
true_private_MS_catch_975<-apply(true_private_MS_catch, 2, quantile, probs=0.975)*2.2046
true_private_MS_catch_025<-apply(true_private_MS_catch, 2, quantile, probs=0.025)*2.2046

#Figure 16 Confidence Band
true_private_MS_season_length_median<-colMedians(true_private_MS_season_length)
true_private_MS_season_length_975<-apply(true_private_MS_season_length, 2, quantile, probs=0.975)
true_private_MS_season_length_025<-apply(true_private_MS_season_length, 2, quantile, probs=0.025)

#Figure 17 Confidence Band
true_private_TX_catch_median<-colMedians(true_private_TX_catch)*2.2046
true_private_TX_catch_975<-apply(true_private_TX_catch, 2, quantile, probs=0.975)*2.2046
true_private_TX_catch_025<-apply(true_private_TX_catch, 2, quantile, probs=0.025)*2.2046

#Figure 18 Confidence Band
true_private_TX_season_length_median<-colMedians(true_private_TX_season_length)
true_private_TX_season_length_975<-apply(true_private_TX_season_length, 2, quantile, probs=0.975)
true_private_TX_season_length_025<-apply(true_private_TX_season_length, 2, quantile, probs=0.025)



#### MSE comparison, need the same data for every scenario. Most of them have already exist in the previous section
#Essential Figures
#Figure 1 Basic Line Chart or The second half of the Shared Dataset
total_catch_median<-colMedians(total_catch)*2.2046

#Figure 2 Basic Line Chart or The second half of the Shared Dataset
total_SSB_median<-colMedians(total_SSB)

#Figure 3 Basic simple with variation
total_catch_MSEcomp<-rowSums(total_catch)*2.2046
total_catch_median_MSEcomp<-median(total_catch_MSEcomp)*2.2046
total_catch_upper_MSEcomp<-as.numeric(quantile(total_catch_MSEcomp,0.975))*2.2046
total_catch_lower_MSEcomp<-as.numeric(quantile(total_catch_MSEcomp,0.025))*2.2046
#calculate variation from Simrun_Num of results

#Figure 4 Basic simple with variation
catch_var_MSEcomp<-rowSds(total_catch)/rowMedians(total_catch)
catch_var_median_MSEcomp<-median(catch_var_MSEcomp)
catch_var_upper_MSEcomp<-as.numeric(quantile(catch_var_MSEcomp,0.975))
catch_var_lower_MSEcomp<-as.numeric(quantile(catch_var_MSEcomp,0.025))
#calculate variation from Simrun_Num of results

#Figure 5 Basic simple with variation
terminal_SSB_MSEcomp<-total_SSB[,ncol(total_SSB)]
terminal_SSB_median_MSEcomp<-median(terminal_SSB_MSEcomp)
terminal_SSB_upper_MSEcomp<-as.numeric(quantile(terminal_SSB_MSEcomp,0.975))
terminal_SSB_lower_MSEcomp<-as.numeric(quantile(terminal_SSB_MSEcomp,0.025))
#calculate variation from Simrun_Num of results

#Figure 6 Basic simple with variation
lowest_SSB_MSEcomp<-rowMins(total_SSB) 
lowest_SSB_median_MSEcomp<-median(lowest_SSB_MSEcomp)
lowest_SSB_upper_MSEcomp<-as.numeric(quantile(lowest_SSB_MSEcomp,0.975))
lowest_SSB_lower_MSEcomp<-as.numeric(quantile(lowest_SSB_MSEcomp,0.025))
#calculate variation from Simrun_Num of results

#Figure 7 Basic bar simple
percent_green_MSEcomp<-sum(Year_green>0)/(Simrun_Num*Runtime_long)
#calculate variation from Simrun_Num of results

#Figure 8 Basic Rader Chart median of the Figure 3-7
median(total_catch_MSEcomp)*2.2046
median(catch_var_MSEcomp)
median(terminal_SSB_MSEcomp)
median(lowest_SSB_MSEcomp)
percent_green_MSEcomp

#Figure 9 Basic simple with variation
total_discards_MSEcomp<-rowSums(comm_discards+recr_discards)*2.2046
total_discards_median_MSEcomp<-median(total_discards_MSEcomp)*2.2046
total_discards_upper_MSEcomp<-as.numeric(quantile(total_discards_MSEcomp,0.975))*2.2046
total_discards_lower_MSEcomp<-as.numeric(quantile(total_discards_MSEcomp,0.025))*2.2046

#Figure 10 Basic simple with variation
discards_var_MSEcomp<-rowSds(comm_discards+recr_discards)/rowMedians(comm_discards+recr_discards)
discards_var_median_MSEcomp<-median(discards_var_MSEcomp)
discards_var_upper_MSEcomp<-as.numeric(quantile(discards_var_MSEcomp,0.975))
discards_var_lower_MSEcomp<-as.numeric(quantile(discards_var_MSEcomp,0.025))

#Other Detailed Figures
#Section 2: Sector Comparison
#Figure 1 Basic Line chart or Shared Dataset Chart
comm_catch_median<-colMedians(comm_catch)*2.2046

#Figure 2 Basic Line chart or Shared Dataset Chart
Forhire_catch_median<-colMedians(Forhire_catch)*2.2046

#Figure 3 Basic Line chart or Shared Dataset Chart
Private_catch_median<-colMedians(Private_catch)*2.2046

#Figure 4 Basic bar simple with variation
total_SSB_first5median<-rowMedians(total_SSB[,1:5])
total_SSB_first5median_median<-median(total_SSB_first5median)
total_SSB_first5median_upper<-as.numeric(quantile(total_SSB_first5median,0.975))
total_SSB_first5median_lower<-as.numeric(quantile(total_SSB_first5median,0.025))

#Figure 5 Basic bar simple with variation
total_SSB_last5median<-rowMedians(total_SSB[,(Runtime_long-4):Runtime_long])
total_SSB_last5median_median<-median(total_SSB_last5median)
total_SSB_last5median_upper<-as.numeric(quantile(total_SSB_last5median,0.975))
total_SSB_last5median_lower<-as.numeric(quantile(total_SSB_last5median,0.025))

#Figure 6 Bar label rotation, maybe with variation?
Forhire_catch_first5median<-median(rowMedians(Forhire_catch[,1:5]))*2.2046
Private_catch_first5median<-median(rowMedians(Private_catch[,1:5]))*2.2046
Forhire_catch_last5median<-median(rowMedians(Forhire_catch[,(Runtime_long-4):Runtime_long]))*2.2046
Private_catch_last5median<-median(rowMedians(Private_catch[,(Runtime_long-4):Runtime_long]))*2.2046

#Figure 7 Bar label rotation, maybe with variation?
true_private_AL_catch_first5median<-median(rowMedians(true_private_AL_catch[,1:5]))*2.2046
true_private_FL_catch_first5median<-median(rowMedians(true_private_FL_catch[,1:5]))*2.2046
true_private_LA_catch_first5median<-median(rowMedians(true_private_LA_catch[,1:5]))*2.2046
true_private_MS_catch_first5median<-median(rowMedians(true_private_MS_catch[,1:5]))*2.2046
true_private_TX_catch_first5median<-median(rowMeans(true_private_TX_catch[,1:5]))*2.2046

#Figure 8 Bar label rotation, maybe with variation?
true_private_AL_catch_last5median<-median(rowMedians(true_private_AL_catch[,(Runtime_long-4):Runtime_long]))*2.2046
true_private_FL_catch_last5median<-median(rowMedians(true_private_FL_catch[,(Runtime_long-4):Runtime_long]))*2.2046
true_private_LA_catch_last5median<-median(rowMedians(true_private_LA_catch[,(Runtime_long-4):Runtime_long]))*2.2046
true_private_MS_catch_last5median<-median(rowMedians(true_private_MS_catch[,(Runtime_long-4):Runtime_long]))*2.2046
true_private_TX_catch_last5median<-median(rowMedians(true_private_TX_catch[,(Runtime_long-4):Runtime_long]))*2.2046

#Figure 9 Bar label rotation, maybe with variation?
true_private_AL_season_length_first5median<-median(rowMedians(true_private_AL_season_length[,1:5]))
true_private_FL_season_length_first5median<-median(rowMedians(true_private_FL_season_length[,1:5]))
true_private_LA_season_length_first5median<-median(rowMedians(true_private_LA_season_length[,1:5]))
true_private_MS_season_length_first5median<-median(rowMedians(true_private_MS_season_length[,1:5]))
true_private_TX_season_length_first5median<-median(rowMedians(true_private_TX_season_length[,1:5]))

#Figure 10 Bar label rotation, maybe with variation?
true_private_AL_season_length_last5median<-median(rowMedians(true_private_AL_season_length[,(Runtime_long-4):Runtime_long]))
true_private_FL_season_length_last5median<-median(rowMedians(true_private_FL_season_length[,(Runtime_long-4):Runtime_long]))
true_private_LA_season_length_last5median<-median(rowMedians(true_private_LA_season_length[,(Runtime_long-4):Runtime_long]))
true_private_MS_season_length_last5median<-median(rowMedians(true_private_MS_season_length[,(Runtime_long-4):Runtime_long]))
true_private_TX_season_length_last5median<-median(rowMedians(true_private_TX_season_length[,(Runtime_long-4):Runtime_long]))

#Figure 11 Basic bar simple with variation
total_catch_first5median<-rowMedians(total_catch[,1:5])*2.2046
total_catch_first5median_median<-median(total_catch_first5median)*2.2046
total_catch_first5median_upper<-as.numeric(quantile(total_catch_first5median,0.975))*2.2046
total_catch_first5median_lower<-as.numeric(quantile(total_catch_first5median,0.025))*2.2046

#Figure 12 Basic bar simple with variation
total_catch_last5median<-rowMedians(total_catch[,(Runtime_long-4):Runtime_long])*2.2046
total_catch_last5median_median<-median(total_catch_last5median)*2.2046
total_catch_last5median_upper<-as.numeric(quantile(total_catch_last5median,0.975))*2.2046
total_catch_last5median_lower<-as.numeric(quantile(total_catch_last5median,0.025))*2.2046

#Figure 13 Basic Line chart or Shared Dataset Chart
comm_discards_median<-colMedians(comm_discards)*2.2046

#Figure 14 Basic Line chart or Shared Dataset Chart
recr_discards_median<-colMedians(recr_discards)*2.2046

#Section 3: Within Sector Comparison
#Figure 1 Basic Rader Chart
comm_catch_var_MSEcomp<-median(rowSds(comm_catch)/rowMedians(comm_catch))
comm_catch_first5median_MSEcomp<-median(rowMedians(comm_catch[,1:5]))*2.2046
comm_catch_last5median_MSEcomp<-median(rowMedians(comm_catch[,(Runtime_long-4):Runtime_long]))*2.2046
total_SSB_last5median_MSEcomp<-median(rowMedians(total_SSB[,(Runtime_long-4):Runtime_long]))
percent_green_MSEcomp<-sum(Year_green>0)/(Simrun_Num*Runtime_long)
comm_dislandratio_ave20_MSEcomp<-median(rowMedians(comm_dislandratio))

#Figure 2 Basic Rader Chart
recr_catch_var_MSEcomp<-median(rowSds(recr_catch)/rowMedians(recr_catch))
Forhire_catch_first5median_MSEcomp<-median(rowMedians(Forhire_catch[,1:5]))*2.2046
Private_catch_first5median_MSEcomp<-median(rowMedians(Private_catch[,1:5]))*2.2046
Forhire_catch_last5median_MSEcomp<-median(rowMedians(Forhire_catch[,(Runtime_long-4):Runtime_long]))*2.2046
Private_catch_last5median_MSEcomp<-median(rowMedians(Private_catch[,(Runtime_long-4):Runtime_long]))*2.2046
total_SSB_last5median_MSEcomp<-median(rowMedians(total_SSB[,(Runtime_long-4):Runtime_long]))
percent_green_MSEcomp<-sum(Year_green>0)/(Simrun_Num*Runtime_long)
recr_dislandratio_ave20_MSEcomp<-median(rowMedians(recr_dislandratio))

#save.image("MSEspace3.RData")

## MSE advanced compare
# change target F and commercial and recreational ratio. make target SSB grey.
# call above simulation for multiple times...
# for every run record total catch, catch variation, terminal SSB, lowest SSB, and year to green.
