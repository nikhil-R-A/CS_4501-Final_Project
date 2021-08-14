#setwd("D:/GoogleDrive/Santiago/Grupo Bioinformatica Estructural y Ambiental URP/Equipo hCoV/Proyecto_hCoV/hCoV_24-04-21/")
#preparation##########################################
## to load libraries, establish hour, and establish theme of graphics
library(lubridate)
library(ggplot2)
library(ggpubr)
library(utils)
library(dplyr)
library(tidyr)
library(ggpmisc)
library(reshape2)
library(stringr)
library(cowplot)
theme_set(theme_bw(base_size = 15))
Sys.setlocale("LC_TIME", "English")
col_conts<-c("Africa"="#7b9e82","Asia"="#c46a6a","Europe"="#6a91c4","NorthAmerica"="#7d7c79","Oceania"="#896ac4","SouthAmerica"="#c4a96a")
#epiFilter##########################################
## Bayesian recursive filtering via EpiFilter
# From: Parag, KV, (2020) "Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves" BioRxiv.

# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion

# Inputs - grid on reproduction numbers (Rgrid), size of grid (m), diffusion noise (eta),
# prior on R (pR0), max time (nday), total infectiousness (Lday), incidence (Iday), confidence (a)

# Output - mean (Rmean), median (Rmed), 50% and 95% quantiles of estimates (Rhat),
# causal posterior over R (pR), pre-update (pRup) and state transition matrix (pstate)

epiFilter <- function(Rgrid, m, eta, pR0, nday, Lday, Iday, a){
  
  # Probability vector for R and prior
  pR = matrix(0, nday, m); pRup = pR
  pR[1, ] = pR0; pRup[1, ] = pR0
  
  # Mean and median estimates
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)
  
  # Initialise mean
  Rmean[1] = pR[1, ]%*%Rgrid
  # CDF of prior 
  Rcdf0 = cumsum(pR0)
  # Initialise quartiles
  idm = which(Rcdf0 >= 0.5, 1); Rmed[1] = Rgrid[idm[1]]
  id1 = which(Rcdf0 >= a, 1); id2 = which(Rcdf0 >= 1-a, 1)
  id3 = which(Rcdf0 >= 0.25, 1); id4 = which(Rcdf0 >= 0.75, 1)
  Rhat[1, 1] = Rgrid[id1[1]]; Rhat[2, 1] = Rgrid[id2[1]]
  Rhat[3, 1] = Rgrid[id3[1]]; Rhat[4, 1] = Rgrid[id4[1]]
  
  # Precompute state distributions for R transitions
  pstate = matrix(0, m, m);
  for(j in 1:m){
    pstate[j, ] = dnorm(Rgrid[j], Rgrid, sqrt(Rgrid)*eta)
  }
  
  # Update prior to posterior sequentially
  for(i in 2:nday){
    # Compute mean from Poisson renewal (observation model)
    rate = Lday[i]*Rgrid
    # Probabilities of observations
    pI = dpois(Iday[i], rate)
    
    # State predictions for R
    pRup[i, ]  = pR[i-1, ]%*%pstate
    # Update to posterior over R
    pR[i, ] = pRup[i, ]*pI
    pR[i, ] = pR[i, ]/sum(pR[i, ])
    
    # Posterior mean and CDF
    Rmean[i] = pR[i, ]%*%Rgrid
    Rcdf = cumsum(pR[i, ])
    
    # Quantiles for estimates
    idm = which(Rcdf >= 0.5, 1); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }
  
  # Main outputs: estimates of R and states
  epiFilter = list(Rmed, Rhat, Rmean, pR, pRup, pstate)
}
#recursPredict##########################################
## Bayesian recursive prediction via EpiFilter
# From: Parag, KV, (2020) "Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves" BioRxiv.
# Notes and assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - can apply causal or smoothing posteriors over R to predict incidence
# - need to set maximum on possible Igrid to interrogate (currently 800)

# Inputs - grid on reproduction numbers (Rgrid), posterior on R (pR), 
# total infectiousness (Lday), mean R esimate using pR (Rmean), confidence level (a and 50%)

# Output - mean prediction (pred) and confidence intervals (predInt)

recursPredict <- function(Rgrid, pR, Lday, Rmean, a){
  
  # Grid size and length of time series
  nday = nrow(pR); m = ncol(pR)
  # Test lengths of inputs
  if (length(Rgrid) != m | length(Lday) != nday){
    stop("Input vectors of incorrect dimension")
  }
  
  # Mean prediction: Lday[i] => Iday[i+1]
  pred = Lday*Rmean; pred = pred[1:length(pred)-1]
  
  # Discrete space of possible predictions
  Igrid = 0:800; lenI = length(Igrid);
  
  # Check if close to upper bound
  if (any(pred > 0.9*max(Igrid))){
    stop("Epidemic size too large")  
  }
  
  # Prediction cdf and quantiles (50% and 95%)
  Fpred = matrix(0, nday-1, lenI)
  predInt = matrix(0, 4, nday-1)
  
  # At every time construct CDF of predictions
  for(i in 1:(nday-1)){
    # Compute rate from Poisson renewal
    rate = Lday[i]*Rgrid
    # Prob of any I marginalised over Rgrid
    pI = rep(0, lenI)
    
    # Probabilities of observations 1 day ahead
    for(j in 1:lenI){
      # Raw probabilities of Igrid
      pIset = dpois(Igrid[j], rate)
      # Normalised by probs of R
      pI[j] = sum(pIset*pR[i, ])
    }
    
    # Quantile predictions and CDF at i+1
    Fpred[i, ] = cumsum(pI)/sum(pI)
    id1 = which(Fpred[i, ] >= a); id2 = which(Fpred[i, ] >= 1-a)
    id3 = which(Fpred[i, ] >= 0.25); id4 = which(Fpred[i, ] >= 0.75)
    
    # Assign prediction results
    predInt[1, i] = Igrid[id1[1]]; predInt[2, i] = Igrid[id2[1]]
    predInt[3, i] = Igrid[id3[1]]; predInt[4, i] = Igrid[id4[1]]
  }
  # Main outputs: mean and 95% predictions
  recursPredict = list(pred, predInt)
}
#epiSmoother##########################################
## Bayesian recursive smoothing via EpiFilter
# From: Parag, KV, (2020) "Improved real-time estimation of reproduction numbers
# at low case incidence and between epidemic waves" BioRxiv.
# Assumptions
# - observation model is Poisson renewal equation (as in EpiEstim)
# - reproduction number state space model is a simple diffusion
# - must have run epiFilter first to obtain forward distribution pR
# - method makes a backward pass to generate qR

# Inputs - grid on reproduction numbers (Rgrid), size of grid (m), filtered posterior (pR),
# update pre-filter (pRup), max time (nday), state transition matrix (pstate), confidence (a)

# Output - mean (Rmean), median (Rmed), lower (Rlow) amd upper (Rhigh) quantiles of estimates,
# smoothed posterior over R (qR) which is backwards and forwards

epiSmoother <- function(Rgrid, m, pR, pRup, nday, pstate, a){
  
  # Last smoothed distribution same as filtered
  qR = matrix(0, nday, m); qR[nday, ] = pR[nday, ]
  
  # Main smoothing equation iteratively computed
  for(i in seq(nday-1, 1)){
    # Remove zeros
    pRup[i+1, pRup[i+1, ] == 0] = 10^-8
    
    # Integral term in smoother
    integ = qR[i+1, ]/pRup[i+1, ]
    integ = integ%*%pstate
    
    # Smoothed posterior over Rgrid
    qR[i, ] = pR[i, ]*integ
    # Force a normalisation
    qR[i, ] = qR[i, ]/sum(qR[i, ]);
  }
  
  # Mean, median estimats of R
  Rmean = rep(0, nday); Rmed = Rmean
  # 50% and 95% (depends on a) confidence on R
  Rhat = matrix(0, 4, nday)
  
  # Compute at every time point
  for (i in 1:nday) {
    # Posterior mean and CDF
    Rmean[i] = qR[i, ]%*%Rgrid
    Rcdf = cumsum(qR[i, ])
    
    # Quantiles for estimates
    idm = which(Rcdf >= 0.5); Rmed[i] = Rgrid[idm[1]]
    id1 = which(Rcdf >= a, 1); id2 = which(Rcdf >= 1-a, 1)
    id3 = which(Rcdf >= 0.25, 1); id4 = which(Rcdf >= 0.75, 1)
    Rhat[1, i] = Rgrid[id1[1]]; Rhat[2, i] = Rgrid[id2[1]]
    Rhat[3, i] = Rgrid[id3[1]]; Rhat[4, i] = Rgrid[id4[1]]
  }
  
  # Main outputs: estimates of R and states
  epiSmoother = list(Rmed, Rhat, Rmean, qR)
}
#Rt_estimation function##########################################
##Preparation of the function to estimate the Rt
RT_estimate<-function(data,country,mut){
  library(dplyr)
  alldata = read.csv(data,sep="\t")
  type<-sapply(strsplit(data,"_"),"[",1)
  alldata_f<-alldata[which(alldata[[type]]==mut),]
  col_1<-paste(type,"confirm",sep="_")
  col_2<-paste("roll_mean",type,sep="_")
  coun_data<-alldata_f%>%dplyr::filter(country=={{country}})%>%
    select(Date,confirm=col_1,roll=col_2)
  for (i in 1:nrow(coun_data)){
    if (coun_data$confirm[i] > 2 & coun_data$confirm[i+1] > 2 & coun_data$confirm[i+2] > 2 ){
      begin<-i
      break
    }
  }
  coun_data_filt<-coun_data[-(1:i-1),]
  coun_data_filt<-coun_data_filt%>%
    dplyr::filter(row_number()<=n()-3)
  for (i in 1:nrow(coun_data_filt)){
    if(round(coun_data_filt$roll[i],0)==0){
      coun_data_filt$roll[i]<-1
    }
  }
  Iday = round(coun_data_filt$roll,0)
  dates  = coun_data_filt$Date
  nday = length(dates); tday = 1:nday
  wdist = dgamma(tday, shape = 2.3669, scale = 2.7463)
  Lday = rep(0, nday) 
  for(i in 2:nday){
    Lday[i] = sum(Iday[seq(i-1, 1, -1)]*wdist[1:(i-1)])    
  }
  Rmin = 0.01; Rmax = 10; eta = 0.1
  m = 200; pR0 = (1/m)*rep(1, m)
  Rgrid = seq(Rmin, Rmax, length.out = m)
  Rfilt = epiFilter(Rgrid, m, eta, pR0, nday, Lday[tday], Iday[tday], 0.025)
  Rsmooth = epiSmoother(Rgrid, m, Rfilt[[4]], Rfilt[[5]], nday, Rfilt[[6]], 0.025)
  final_data<-data.frame(Rt=Rsmooth[[3]],date=as.Date(dates,format="%Y-%m-%d"),l_975=Rsmooth[[2]][1,],u_975=Rsmooth[[2]][2,],l_75=Rsmooth[[2]][3,],u_75=Rsmooth[[2]][4,],country=country,type=type,mut=mut)
  return(final_data)
}

#Prepare metadata of genomes##########################################
##
#read the csv file with metadata
my_data<-read.csv(file = "metadata.tsv", sep = "\t", header = TRUE)
hq_ids<-as.list(read.csv(file="aligned_290Ns_minusbad_align_clear_ids.list",header=FALSE))
hq_ids_name<-sapply(strsplit(hq_ids[["V1"]],"\\|"),"[",1)
hq_data<-my_data%>%dplyr::filter(Virus.name%in%hq_ids_name&Host%in%c("Human","Environment"))%>%
  dplyr::select(Virus_name=Virus.name,EPI_ISL=Accession.ID,date=Collection.date,Location=Location,Host=Host,Submission_date=Submission.date)
#Fix the format of the dates, continent and country
hq_data$year<-sapply(strsplit(hq_data$date, "-"), "[", 1)
hq_data$month<-sapply(strsplit(hq_data$date, "-"), "[", 2)
hq_data$month<-as.numeric(hq_data$month)
hq_data$Month<-as.character(lubridate::month(ymd(010101) + months(hq_data$month-1),label=TRUE,abbr=TRUE))
hq_data$y_m<-paste(hq_data$year,hq_data$Month,sep="-")
hq_data$y_m<-as.factor(hq_data$y_m)
hq_data$y_m<-factor(hq_data$y_m, levels = c("2019-Dec","2020-Jan","2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar","2021-Apr"),ordered=TRUE)
hq_data$Location<-gsub(" ", "",hq_data$Location)
hq_data$continent<-sapply(strsplit(hq_data$Location, "/"), "[", 1)
hq_data$country<-sapply(strsplit(hq_data$Location, "/"), "[", 2)
#fix countries
hq_data$country<-gsub("HongKong", "China", hq_data$country)
hq_data$country<-gsub("Guadeloupe", "LeewardIslands", hq_data$country)
hq_data$country<-gsub("RepublicoftheCongo", "Congo", hq_data$country)
hq_data$country<-gsub("Mayotte", "Madagascar", hq_data$country)
hq_data$country<-gsub("Reunion", "Madagascar", hq_data$country)
#save the metadata here to use in other analyses
hq_data_analysis<-hq_data
hq_data_analysis$strain<-paste(hq_data_analysis$Virus_name,hq_data_analysis$date,hq_data_analysis$Submission_date,sep="|")
hq_data_analysis<-hq_data_analysis[,c(13,1,2,3,4,5,6,7,8,9,10,11,12)]
write.table(hq_data_analysis,"aligned_290Ns_minusbad_align_clear_analysis_metadata.tsv",sep="\t",row.names=FALSE,quote=FALSE)
writeLines(hq_data_analysis$strain,"aligned_290Ns_minusbad_align_clear_analysis_ids.txt")
#Generate the files to calculate the relative frequency of nucleotides in each genomic position by country-month##########################################
##
#Select country-months with 15 or more genomes
table_counts<-hq_data%>%
  dplyr::group_by(country,y_m)%>%
  dplyr::summarise(count=n())%>%
  dplyr::filter(count>=15)%>%
  tidyr::drop_na()
hq_data_filtered<-list()
for (x in levels(as.factor(table_counts$country))){
  test<-table_counts%>%dplyr::filter(country==x)
  for (i in test$y_m){
    name<-paste(x,i,sep="_")
    hq_data_filtered[[name]]<-hq_data%>%dplyr::filter(country==x&y_m==i)
  }
}
hq_data_filtered_final<-do.call("rbind",hq_data_filtered)
#generate the files to perform the NRFp by country
setwd("D:/GoogleDrive/Santiago/Grupo Bioinformatica Estructural y Ambiental URP/Equipo hCoV/Proyecto_hCoV/hCoV_24-04-21/Frequency_files/")
script_total<-list()
script_total<-append(script_total,"#!/bin/bash")
script_total<-append(script_total,"#SBATCH -p general")
script_total<-append(script_total,"#SBATCH --nodes=1")
script_total<-append(script_total,"#SBATCH --ntasks-per-node=5")
script_total<-append(script_total,"#SBATCH -o log_%j")
hq_<-list()
hq_ym_c_<-list()
countries<-list()
for (i in levels(as.factor(hq_data_filtered_final$country))){
  script_total<-append(script_total,paste("sbatch ",i,".slurm",sep=""))
  script<-list()
  script<-append(script,"#!/bin/bash")
  script<-append(script,"#SBATCH -p general")
  script<-append(script,"#SBATCH --nodes=1")
  script<-append(script,"#SBATCH --ntasks-per-node=5")
  script<-append(script,"#SBATCH -o log_%j")
  hq_[[i]]<-hq_data_filtered_final%>%dplyr::filter(country==i)
  for (x in levels(as.factor(hq_data_filtered_final$y_m))){
    hq_ym_c_[[x]]<-hq_[[i]]%>%dplyr::filter(hq_[[i]]$y_m==x)
    name<-paste(i,"_",x,"_ids.txt",sep="")
    name_2<-paste("python3 retrieve_sequences_from_txt_v2.py ",i,"_",x,"_ids.txt"," aligned_290Ns_minusbad_align_clear.fasta"," && python3 frequencies_v2.py ",i,"_",x,"_ids_seqs.fasta",sep="")
    if (length(hq_ym_c_[[x]]$Virus_name)>=15){
      hq_ym_c_[[x]]$new<-paste(hq_ym_c_[[x]]$Virus_name,hq_ym_c_[[x]]$date,hq_ym_c_[[x]]$Submission_date,sep="|")
      writeLines(hq_ym_c_[[x]]$new,name)
      script<-append(script,name_2)
      hq_ym_c_[[x]]$country<-as.factor(hq_ym_c_[[x]]$country)
      for (z in levels(hq_ym_c_[[x]]$country)){
        if (!(z %in% countries)){
          countries<-append(countries,z)
        }
      }
    }
    name_3<-paste(i,".slurm",sep="")
    writeLines(text=as.character(script),name_3)  
  }
}
writeLines(text=as.character(script_total),"NRFp.slurm") 
#those created files will be used in the next part
#Obtaining the number of cases by country-month##########################################
##
#download the .csv file from: https://www.ecdc.europa.eu/en/publications-data/data-national-14-day-notification-rate-covid-19
#the downloaded file is named "descarga", change the name to: cdc_europe_cases_29_04_20.csv and use the function "read.csv" to read it
setwd("D:/GoogleDrive/Santiago/Grupo Bioinformatica Estructural y Ambiental URP/Equipo hCoV/Proyecto_hCoV/hCoV_24-04-21/")
data <- read.csv("cdc_europe_cases_29_04_20.csv", sep = ",")
names(data)[1]<-c("country")
#fix names
data$country<-gsub("Bonaire, Saint Eustatius And Saba", "Bonaire", data$country)
data$country<-gsub("CuraÃ§ao", "Curacao", data$country)
data$country<-gsub("British Virgin Islands", "LeewardIslands", data$country)
data$country<-gsub("United States Virgin Islands", "LeewardIslands", data$country)
data$country<-gsub("Anguilla", "LeewardIslands", data$country)
data$country<-gsub("Antigua and Barbuda", "LeewardIslands", data$country)
data$country<-gsub("Montserrat", "LeewardIslands", data$country)
data$country<-gsub("Timor Leste", "Timor-Leste", data$country)
data$country<-gsub("United States Of America", "USA", data$country)
data$country<-gsub("United Kingdom", "UnitedKingdom", data$country)
data$country<-gsub("Democratic Republic Of The Congo", "DemocraticCongo", data$country)
#verify that names in GISAID metadata match with those from cases
data$country<-as.factor(data$country)
country_cases_levels<-as.list(levels(data$country))
country_data_my_data<-list()
for (i in countries){
  if (i %in% country_cases_levels){
    country_data_my_data<-append(country_data_my_data,i)
  }
}
#This part will print the countries that are not in the metadata of cases but are present in GISAID
for (i in countries){
  if (!(i %in% country_data_my_data)){
    print(i)
  }
}
#if list countries and country_data_my_data have the same length and the previous section did not print anything, then proceed to the next step
##
#take the number of cases of the country-month combinations selected from GISAID
data_cases_rdy<-data%>%dplyr::filter(country%in%country_data_my_data)%>%
  dplyr::select(country,continent,weekly_count,year_week)
#create lists with the countries belonging to each of the american subcontinents
#NA
countries_NA<-hq_data_filtered_final%>%dplyr::filter(continent=="NorthAmerica")
countries_list_NA<-as.list(levels(as.factor(countries_NA$country)))
#SA
countries_SA<-hq_data_filtered_final%>%dplyr::filter(continent=="SouthAmerica")
countries_list_SA<-as.list(levels(as.factor(countries_SA$country)))
#change the names of the continent America to North America, Central America or South America depending on the countries in the lists
data_cases_rdy<-within(data_cases_rdy,continent[continent=="America"&country%in%countries_list_NA]<-"NorthAmerica")
data_cases_rdy<-within(data_cases_rdy,continent[continent=="America"&country%in%countries_list_SA]<-"SouthAmerica")
#arrange the information of months and years
data_cases_rdy$year_week_new<-as.Date(paste("1", data_cases_rdy$year_week, sep = "-"), format = "%w-%Y-%W")
data_cases_rdy$year_week_new<-as.character(data_cases_rdy$year_week_new)
data_cases_rdy$year<-sapply(strsplit(data_cases_rdy$year_week_new, "-"), "[", 1)
data_cases_rdy$month<-sapply(strsplit(data_cases_rdy$year_week_new, "-"), "[", 2)
data_cases_rdy$month<-as.numeric(data_cases_rdy$month)
data_cases_rdy$Month<-as.character(lubridate::month(ymd(010101)+months(data_cases_rdy$month-1),label=TRUE,abbr=TRUE))
data_cases_rdy$y_m<-paste(data_cases_rdy$year,data_cases_rdy$Month,sep="-")
data_cases_rdy$y_m<-gsub("NA-NA", "2020-Dec",data_cases_rdy$y_m)
data_cases_rdy$y_m<-as.factor(data_cases_rdy$y_m)
data_cases_rdy$y_m<-factor(data_cases_rdy$y_m, levels = c("2019-Dec","2020-Jan","2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar","2021-Apr"))
#create a data frame with cases by country by month and bind it to genome counts
cases_c_m<-data_cases_rdy%>%
  dplyr::group_by(continent,country,y_m)%>%
  dplyr::summarise(count_cases=sum(weekly_count))
hq_data_c_m<-hq_data_filtered_final%>%
  dplyr::group_by(continent,country,y_m)%>%
  dplyr::summarise(count_genomes=(n()))
hq_data_c_m_rdy<-hq_data_c_m%>%
  dplyr::filter(count_genomes>=15)
cases_genomes_rdy<-dplyr::left_join(hq_data_c_m_rdy,cases_c_m)%>%
  drop_na()
cases_genomes_rdy_cont<-cases_genomes_rdy%>%
  dplyr::group_by(continent,y_m)%>%
  dplyr::summarise(count_genomes_cont=sum(count_genomes),count_cases_cont=sum(count_cases))
#plotear casos vs genomas
p_cases_genomes<-ggscatter(cases_genomes_rdy,x="count_genomes",y="count_cases",
                           add="reg.line",cor.method="spearman",conf.int=TRUE,
                           xlab="Number of genomes",ylab = "Number of cases")+stat_cor(aes(label=paste(..r.label..,..rr.label..,..p.label..,sep="~`,`~")),method="spearman")
ggsave("cas_gen_24-04-21.jpg",plot=p_cases_genomes,width=7,height=5,dpi=500)
#Estimate the relative frequencies of cases by mutation (NRFp)##########################################
##
#read the files with the frequencies calculated
freqs<-list()
tables<-list()
list_freqs<-dir(pattern="*._ids_seqs.fasta.freq")
for (i in list_freqs){
  freqs[[i]]<-read.csv(i,sep="\t",header=TRUE)
  if (length(freqs[[i]]$X.) == 0){
    freqs[[i]]$X.<-0
  }
  if (length(freqs[[i]]$N) == 0){
    freqs[[i]]$N<-0
  }
  table<-freqs[[i]]%>%dplyr::select("Position","A","C","G","T","N","X.")
  table$Position<-table$Position+203
  table_rdy<-reshape2::melt(table,id.var=c("Position"))
  tables[[i]]<-table_rdy
}
All_freqs<-table_rdy%>%dplyr::select(Position,variable)
cases_genomes_rdy$rdy_names<-paste(cases_genomes_rdy$continent,cases_genomes_rdy$country,cases_genomes_rdy$y_m,sep="_")
for (i in list_freqs){
  name<-paste(word(i,1,sep="_"),word(i,2,sep="_"),sep="_")
  for (k in cases_genomes_rdy$rdy_names){
    if (paste(word(k,2,sep="_"),word(k,3,sep="_"),sep="_") == name){
      tables[[i]]$cases_value<-tables[[i]]$value*cases_genomes_rdy$count_cases[cases_genomes_rdy$rdy_names==paste(word(k,1,sep="_"),name,sep="_")]
      new_table<-tables[[i]]%>%dplyr::select(Position,variable,cases_value)
      names(new_table)[3]<-paste(word(k,1,sep="_"),name,sep="_")
      All_freqs<-dplyr::full_join(All_freqs,new_table)
    }
  }
}  
All_freqs_total<-All_freqs%>%mutate(NRFp=rowSums(.[3:ncol(All_freqs)])/sum(cases_genomes_rdy$count_cases))
write.csv(All_freqs_total,"All_freqs_total.csv",row.names=FALSE)
#from the next line you can start reading the file created in the previous line
All_freqs_total<-read.csv("All_freqs_total.csv")
filename_2<-"D:/GoogleDrive/Santiago/Grupo Bioinformatica Estructural y Ambiental URP/Equipo hCoV/Proyecto_hCoV/8 PostAnalysis/Identities/Genomes/2_root_id/2_root_final_seqs_090620_align.fasta.id"
ref<-read.csv(filename_2, sep = "\t", header = TRUE)
ref$ref<-toupper(ref$ref)
ref<-ref%>%dplyr::select(Position=posit,ref)
ref$Position<-(ref$Position+202)
All_freqs_total<-left_join(All_freqs_total,ref,by=c("Position"))
#this part is to plot NRFp
setwd("D:/GoogleDrive/Santiago/Grupo Bioinformatica Estructural y Ambiental URP/Equipo hCoV/Proyecto_hCoV/hCoV_24-04-21/Frequencies")
All_freqs_total$variable<-gsub("X.","-",All_freqs_total$variable)
theme_set(theme_bw(base_size = 12))
All_freqs_total$variable<-factor(All_freqs_total$variable,levels=c("A","C","G","T","-","N"))
cols<-c("A"="chartreuse3","G"="gray80","C"="blue3","T"="red3","-"="black","N"="pink")
cols_int<-c("A"="chartreuse3","G"="gray80","C"="blue3","T"="red3","gap"="black","N"="pink")
#NRFp > 0.03
p<-ggplot(All_freqs_total,aes(Position,NRFp,fill=variable,color=variable))+ggtitle("World")+geom_point(data=subset(All_freqs_total,variable!=ref&NRFp>=0.03&variable!="N"),size=3,alpha=0.6,shape=21)+scale_x_continuous(breaks=seq(0,30000,1000))+scale_y_continuous(limits=c(0,1),breaks=seq(0,1.0,0.1))+scale_fill_manual(values=cols)+scale_color_manual(values=cols)+labs(x="Position",y=expression(paste(NRF[p],"")))+theme(axis.text.x=element_text(angle=45,hjust=0.5,vjust=0.5),legend.title=element_blank(),legend.position="bottom")+guides(fill=guide_legend(nrow=1))
ggsave("Figure_S1.jpg",plot=p,width=7.5,height=5,dpi=500)
#First observation of the dynamic of RM mutations (0.03 < NRFp)##########################################
##
#tabla de mutaciones con igual o más de 0.03 NRFp
RM_mutations<-All_freqs_total%>%
  dplyr::filter(variable!=ref&NRFp>=0.03&variable!="N")%>%
  dplyr::select(Position,variable,ref,NRFp)
RM_mutations<-RM_mutations%>%dplyr::mutate(ID=paste(RM_mutations$ref,RM_mutations$Position,RM_mutations$variable,sep=""))
write.table(RM_mutations,file="RM_table.txt",sep="\t",quote=FALSE,row.names=FALSE)          
#Other features of mutations were manually annotated and saved in a file named "AA_change.tsv"
RM_annotations<-read.csv("AA Change.tsv",sep="\t")
RM_annotations$code<-gsub("X.","-",RM_annotations$code)
#change ?? to triangle
RM_annotations$AA_Change<-gsub("del","??",RM_annotations$AA_Change)
RM_annotations$label<-paste(RM_annotations$code,RM_annotations$AA_Change,RM_annotations$Location,sep="_")
RM_annotations<-RM_annotations[order(RM_annotations$Position),]
RM_annotations$NRFp<-round(RM_annotations$NRFp,3)
RM_labels<-RM_annotations%>%
  dplyr::select(Position,code,label,AA_Change,Location)
#creation of table with metadata of RM mutations with NRFp data
RM_analysis<-All_freqs_total%>%
  dplyr::filter(NRFp>=0.03&variable!=ref&variable!="N")%>%
  dplyr::mutate(code=paste(ref,Position,variable,sep=""))%>%
  dplyr::select(!(c(NRFp,ref,Position,variable)))
t_RM_analysis<-reshape::melt(RM_analysis,id.vars="code")
t_RM_analysis$continent<-sapply(strsplit(as.character(t_RM_analysis$variable),"_"),"[",1)
t_RM_analysis$country<-sapply(strsplit(as.character(t_RM_analysis$variable),"_"),"[",2)
t_RM_analysis$y_m<-sapply(strsplit(as.character(t_RM_analysis$variable),"_"),"[",3)
t_RM_analysis$y_m<-gsub("\\.","-",t_RM_analysis$y_m)
t_RM_analysis$variable<-NULL
#analysis of the dynamic of RM mutations
t_RM_analysis_comp<-dplyr::full_join(t_RM_analysis,cases_genomes_rdy)%>%
  dplyr::mutate(freq=value/count_cases)%>%
  dplyr::mutate(CI=1.96*(sqrt(((1-freq)*(freq))/count_genomes)))%>%
  dplyr::mutate_all(~replace(.,is.nan(.),0))
t_RM_analysis_comp$y_m<-factor(t_RM_analysis_comp$y_m,levels=c("2020-Jan","2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar","2021-Apr"),ordered=TRUE)
t_RM_analysis_comp$code<-gsub("X.","-",t_RM_analysis_comp$code)
t_RM_analysis_comp$code<-as.factor(t_RM_analysis_comp$code)
#seleccionar paises para hacer los analisis por paises
countries_diff_freqs<-list()
#the next loop selects the countries that will be included in the analysis by countries
for (c in levels(as.factor(cases_genomes_rdy$country))){
  country_table<-cases_genomes_rdy%>%dplyr::filter(country==c)%>%dplyr::select(country,y_m)
  if (nrow(country_table)>=14){
    countries_diff_freqs<-append(countries_diff_freqs,c)
  }
}
#creation of data of RM mutations by seected countries
t_RM_count_count<-t_RM_analysis_comp%>%
  dplyr::filter(country%in%countries_diff_freqs)
t_RM_count_count$code<-gsub("X.","-",t_RM_count_count$code)
t_RM_count_count_labels<-left_join(t_RM_count_count,RM_labels)
#creation of data of RM mutations by continents
t_RM_cont_list<-list()
for (i in levels(as.factor(t_RM_analysis_comp$continent))){
  df<-t_RM_analysis_comp%>%
    dplyr::filter(continent==i)%>%
    dplyr::group_by(continent,y_m,code)%>%
    dplyr::mutate(count_genomes_cont=sum(count_genomes),count_cases_cont=sum(count_cases))%>%
    dplyr::summarise(NRFp=wtd.mean(freq,count_cases),NRFp_sd=wtd.var(freq,count_cases))
  t_RM_cont_list[[i]]<-df
}
t_RM_cont_v2<-do.call(rbind,t_RM_cont_list)
t_RM_cont_v2$code<-gsub("X.","-",t_RM_cont_v2$code)
t_RM_cont_labels_v2<-dplyr::left_join(t_RM_cont_v2,RM_labels)
#creation of data of RM mutations by world
t_RM_world_v2<-t_RM_analysis_comp%>%
  dplyr::group_by(y_m,code)%>%
  dplyr::mutate(count_genomes_cont=sum(count_genomes),count_cases_cont=sum(count_cases))%>%
  dplyr::summarise(NRFp=wtd.mean(freq,count_cases),NRFp_sd=wtd.var(freq,count_cases))%>%
  dplyr::mutate(group="Global")
t_RM_world_v2$code<-gsub("X.","-",t_RM_world_v2$code)
t_RM_world_labels_v2<-dplyr::left_join(t_RM_world_v2,RM_labels)
#Calculation of frequency change of RM mutations globally##########################################
##
#Calculating frequency changes
t_RM_world_v2$code<-as.factor(t_RM_world_v2$code)
dif_freqs_table_world_v2<-data.frame(y_m="test",NRFp_change=0,code="test")
for (z in levels(factor(t_RM_world_v2$code))){
  frame<-t_RM_world_v2%>%dplyr::filter(code==z)
  for (i in 1:length(levels(frame$y_m))){
    if (i != 1){
      for (y in 2:nrow(frame)){
        if (frame[y,"y_m"] == levels(frame$y_m)[i] && frame[y,"y_m"] != "2021-Apr"){
          l<-as.numeric(frame[y,"NRFp"])-as.numeric(frame[y-1,"NRFp"])
          query<-data.frame(y_m=levels(frame$y_m)[i-1],NRFp_change=l,code=z)
          dif_freqs_table_world_v2<-rbind(dif_freqs_table_world_v2,query)
        }
      }
    }
  }
}
for (i in 1:nrow(dif_freqs_table_world_v2)){
  if (dif_freqs_table_world_v2$code[i]=="test"){
    dif_freqs_table_world_v2<-dif_freqs_table_world_v2[-i,]
  }
}
dif_freqs_table_world_v2$y_m<-factor(dif_freqs_table_world_v2$y_m,levels=c("2020-Jan","2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar","2021-Apr"),ordered=TRUE)
dif_freqs_table_world_v2$code<-as.factor(dif_freqs_table_world_v2$code)
#Assignment of HF, MF and LF mutations##########################################
##
#creation of High Frequency (HF) tables and list
dif_freqs_list_world_HF_v2<-list()
for (i in 1:length(levels(dif_freqs_table_world_v2$code))){
  code_df<-dif_freqs_table_world_v2%>%
    dplyr::filter(code==levels(dif_freqs_table_world_v2$code)[i])
  if (max(code_df$NRFp_change)>0.10 && min(code_df$NRFp_change) > (-0.02)){
    dif_freqs_list_world_HF_v2[[levels(dif_freqs_table_world_v2$code)[i]]]<-code_df
  }
}
dif_freqs_world_HF_v2<-do.call(rbind,dif_freqs_list_world_HF_v2)
t_RM_world_HF_list_v2<-levels(as.factor(as.character(dif_freqs_world_HF_v2$code)))
t_RM_world_HF_v2<-t_RM_world_labels_v2%>%dplyr::filter(code%in%t_RM_world_HF_list_v2)
t_RM_world_HF_v2<-t_RM_world_HF_v2[order(t_RM_world_HF_v2$Position),]
t_RM_world_HF_v2$code<-as.factor(as.character(t_RM_world_HF_v2$code))
t_RM_world_notHF_v2<-t_RM_world_labels_v2%>%dplyr::filter(!(code%in%t_RM_world_HF_list_v2))
t_RM_world_notHF_v2$code<-as.factor(as.character(t_RM_world_notHF_v2$code))
#creation of low frequency advantageous (also HFs) tables and list
dif_freqs_list_world_LFa_v2<-list()
for (i in 1:length(levels(dif_freqs_table_world_v2$code))){
  code_df<-dif_freqs_table_world_v2%>%
    dplyr::filter(code==levels(dif_freqs_table_world_v2$code)[i])
  if (min(code_df$NRFp_change) > (-0.01) && max(code_df$NRFp_change) > 0.05){
    dif_freqs_list_world_LFa_v2[[levels(dif_freqs_table_world_v2$code)[i]]]<-code_df
  }
}
dif_freqs_world_LFa_v2<-do.call(rbind,dif_freqs_list_world_LFa_v2)
t_RM_world_LFa_v2<-t_RM_world_notHF_v2%>%dplyr::filter(code%in%as.character(dif_freqs_world_LFa_v2$code))
t_RM_world_LFa_list_v2<-levels(as.factor(as.character(t_RM_world_LFa_v2$code)))
t_RM_world_LFa_v2<-t_RM_world_LFa_v2[order(t_RM_world_LFa_v2$Position),]
t_RM_world_LFa_v2$code<-as.factor(as.character(t_RM_world_LFa_v2$code))
t_RM_world_notHF_notLFa_v2<-t_RM_world_notHF_v2%>%dplyr::filter(!(code%in%as.character(dif_freqs_world_LFa_v2$code)))
t_RM_world_notHF_notLFa_v2$code<-as.factor(as.character(t_RM_world_notHF_notLFa_v2$code))
#creation of Medium Frequency (MF) tables and list
t_RM_world_list_MF_v2<-list()
for (i in 1:length(levels(as.factor(as.character(t_RM_world_notHF_notLFa_v2$code))))){
  code_df<-t_RM_world_notHF_notLFa_v2%>%
    dplyr::filter(code==levels(t_RM_world_notHF_notLFa_v2$code)[i])
  if (max(code_df$NRFp)>0.15){
    t_RM_world_list_MF_v2[[levels(t_RM_world_notHF_notLFa_v2$code)[i]]]<-code_df
  }
}
t_RM_world_MF_v2<-do.call(rbind,t_RM_world_list_MF_v2)
t_RM_world_MF_v2<-t_RM_world_MF_v2[order(t_RM_world_MF_v2$Position),]
t_RM_world_MF_v2$code<-as.factor(as.character(t_RM_world_MF_v2$code))
t_RM_world_MF_v2<-t_RM_world_MF_v2%>%
  dplyr::filter(!(Position%in%c(210,23604)))
t_RM_world_MF_list_v2<-levels(as.factor(as.character(t_RM_world_MF_v2$code)))
#creation of Low Frequency (LF) tables and list
t_RM_world_list_LF_v2<-list()
for (i in 1:length(levels(as.factor(as.character(t_RM_world_notHF_notLFa_v2$code))))){
  code_df<-t_RM_world_notHF_notLFa_v2%>%
    dplyr::filter(code==levels(t_RM_world_notHF_notLFa_v2$code)[i])
  if (max(code_df$NRFp)<=0.15){
    t_RM_world_list_LF_v2[[levels(t_RM_world_notHF_notLFa_v2$code)[i]]]<-code_df
  }
}
t_RM_world_LF_v2<-do.call(rbind,t_RM_world_list_LF_v2)
t_RM_world_LF_v2<-t_RM_world_LF_v2[order(t_RM_world_LF_v2$Position),]
t_RM_world_LF_v2$code<-as.factor(as.character(t_RM_world_LF_v2$code))
t_RM_world_LF_v2<-t_RM_world_LF_v2%>%
  dplyr::filter(!(Position%in%c(203,204)))
t_RM_world_LF_list_v2<-levels(as.factor(as.character(t_RM_world_LF_v2$code)))
#test plots
p_int_RM_world_HF<-ggplot(subset(t_RM_world_HF_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=44))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_HF)
p_int_RM_world_notHF<-ggplot(subset(t_RM_world_notHF_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=73))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_notHF)
p_int_RM_world_LFa<-ggplot(subset(t_RM_world_LFa_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=23))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_LFa)
p_int_RM_world_notHF_notLFa<-ggplot(subset(t_RM_world_notHF_notLFa_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=50))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_notHF_notLFa)
p_int_RM_world_MF<-ggplot(subset(t_RM_world_MF_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=29))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_MF)
p_int_RM_world_LF<-ggplot(subset(t_RM_world_LF_v2,y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0),axis.title.x=element_blank())+labs(y="Normalized Relative Frequency")+scale_color_manual(values=rainbow(n=17))+theme(legend.position="none")+scale_y_continuous(limits=c(0,1))
ggplotly(p_int_RM_world_LF)
ggarrange(p_int_RM_world_HF,p_int_RM_world_LFa,p_int_RM_world_MF,p_int_RM_world_LF)
#add mutation type to the table of NRFps
t_RM_world_labels_v2$Type<-"any"
for (i in 1:nrow(t_RM_world_labels_v2)){
  if (t_RM_world_labels_v2$code[i] %in% t_RM_world_HF_list_v2){
    t_RM_world_labels_v2$Type[i]<-"HF"
  }else if(t_RM_world_labels_v2$code[i] %in% t_RM_world_LFa_list_v2){
    t_RM_world_labels_v2$Type[i]<-"LFa"
  }else if(t_RM_world_labels_v2$code[i] %in% t_RM_world_MF_list_v2){
    t_RM_world_labels_v2$Type[i]<-"MF"
  }else if(t_RM_world_labels_v2$code[i] %in% t_RM_world_LF_list_v2){
    t_RM_world_labels_v2$Type[i]<-"LF"
  }
}
#plot by type of mutations
pHF<-ggplot(data=subset(t_RM_world_labels_v2,Type%in%c("HF","LFa")&y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line(size=3,alpha=0.2)+scale_y_continuous(limits=c(0,1))+theme(legend.position = "none",axis.text=element_text(size=15),axis.text.x=element_text(angle=90),axis.title.x=element_blank(),title=element_text(size=20))+ggtitle("High-Frequency mutations")+labs(y="Normalized Relative Frequency (NRFp)")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1),alpha=0.2),show.legend=FALSE,width=.2,size=1)+annotate("text",x=2,y=0.875,label="Group 1",size=12)+annotate("text",x=12,y=0.5,label="Group 2",size=12)
pMF<-ggplot(data=subset(t_RM_world_labels_v2,Type=="MF"&y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line(size=3,alpha=0.2)+scale_y_continuous(limits=c(0,1))+theme(legend.position = "none",axis.text=element_text(size=15),axis.text.x=element_text(angle=90),axis.title.x=element_blank(),title=element_text(size=20))+ggtitle("Medium-Frequency mutations")+labs(y="Normalized Relative Frequency (NRFp)")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1),alpha=0.2),show.legend=FALSE,width=.2,size=1)
pLF<-ggplot(data=subset(t_RM_world_labels_v2,Type=="LF"&y_m!="2021-Apr"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line(size=3,alpha=0.2)+scale_y_continuous(limits=c(0,1))+theme(legend.position = "none",axis.text=element_text(size=15),axis.text.x=element_text(angle=90),axis.title.x=element_blank(),title=element_text(size=20))+ggtitle("Low-Frequency mutations")+labs(y="Normalized Relative Frequency (NRFp)")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1),alpha=0.2),show.legend=FALSE,width=.2,size=1)
p_temp_dyn<-ggarrange(pHF,pMF,pLF,ncol=1,nrow=3)+draw_plot_label(label=c("a","b","c"),size=40,x=0,y=c(1.01,0.67666,0.34333))
ggsave(filename="Figure_1.jpg",plot=p_temp_dyn,dpi=500,width=10.5,height=20)
#Plots expanded for supplementary figures
t_RM_world_labels_v2$label<-factor(t_RM_world_labels_v2$label,levels=unique(t_RM_world_labels_v2$label[order(t_RM_world_labels_v2$Position,t_RM_world_labels_v2$label)]),ordered=TRUE)

p_RM_HF_nonSyn_nsp<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("HF")&grepl("nsp",Location)&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_HF_nonSyn_nsp_fin<-p_RM_HF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_nonSyn_nsp_last<-p_RM_HF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_nonSyn_Str<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("HF")&!(grepl("nsp",Location))&!(grepl("UTR",Location))&!(grepl("Intergenic",Location))&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_HF_nonSyn_Str_fin<-p_RM_HF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_nonSyn_Str_last<-p_RM_HF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_Syn<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("HF")&AA_Change=="Syn"|Position%in%c(241,28271)&y_m!="2021-Apr"&Type%in%c("HF")),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_HF_Syn_fin<-p_RM_HF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_Syn_last<-p_RM_HF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_HF_nonSyn_nsp_legend<-get_legend(p_RM_HF_nonSyn_nsp_fin)
p_RM_HF_nonSyn_Str_legend<-get_legend(p_RM_HF_nonSyn_Str_fin)
p_RM_HF_Syn_legend<-get_legend(p_RM_HF_Syn_fin)
p_RM_HF<-ggarrange(p_RM_HF_nonSyn_nsp_last,p_RM_HF_nonSyn_Str_last,p_RM_HF_Syn_last,p_RM_HF_nonSyn_nsp_legend,p_RM_HF_nonSyn_Str_legend,p_RM_HF_Syn_legend,heights=c(2,1),ncol=3,nrow=2)+draw_plot_label(label=c("a","b","c"),size=20,x=c(0,0.3333,0.6666),y=c(1.05,1.05,1.05))

p_RM_LFa_nonSyn_nsp<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LFa")&grepl("nsp",Location)&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_LFa_nonSyn_nsp_fin<-p_RM_LFa_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_nonSyn_nsp_last<-p_RM_LFa_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_nonSyn_Str<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LFa")&!(grepl("nsp",Location))&!(grepl("UTR",Location))&!(grepl("Intergenic",Location))&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_LFa_nonSyn_Str_fin<-p_RM_LFa_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_nonSyn_Str_last<-p_RM_LFa_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_Syn<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LFa")&AA_Change=="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_LFa_Syn_fin<-p_RM_LFa_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_Syn_last<-p_RM_LFa_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LFa_nonSyn_nsp_legend<-get_legend(p_RM_LFa_nonSyn_nsp_fin)
p_RM_LFa_nonSyn_Str_legend<-get_legend(p_RM_LFa_nonSyn_Str_fin)
p_RM_LFa_Syn_legend<-get_legend(p_RM_LFa_Syn_fin)
p_RM_LFa<-ggarrange(p_RM_LFa_nonSyn_nsp_last,p_RM_LFa_nonSyn_Str_last,p_RM_LFa_Syn_last,p_RM_LFa_nonSyn_nsp_legend,p_RM_LFa_nonSyn_Str_legend,p_RM_LFa_Syn_legend,heights=c(2,1),ncol=3,nrow=2)+draw_plot_label(label=c("d","e","f"),size=20,x=c(0,0.3333,0.6666),y=c(1.05,1.05,1.05))

text_1<-ggplot()+ggplot2::annotate("text",x=4, y = 25, size=8, label = "Non-synonymous nsp")+theme_void()
text_2<-ggplot()+ggplot2::annotate("text",x=4, y = 25, size=8, label = "Non-synonymous Structural")+theme_void()
text_3<-ggplot()+ggplot2::annotate("text",x=4, y = 25, size=8, label = "Synonymous")+theme_void()
titles<-ggarrange(text_1,text_2,text_3,ncol=3,nrow=1)

p_RM_HF_LFa<-ggarrange(titles,p_RM_HF,p_RM_LFa,nrow=3,heights=c(0.7,6,6))
ggsave("Figure_S2.jpg",plot=p_RM_HF_LFa,dpi=500,width=30,height=15)

p_RM_MF_nonSyn_nsp<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("MF")&grepl("nsp",Location)&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_MF_nonSyn_nsp_fin<-p_RM_MF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_nonSyn_nsp_last<-p_RM_MF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_nonSyn_Str<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("MF")&!(grepl("nsp",Location))&!(grepl("UTR",Location))&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_MF_nonSyn_Str_fin<-p_RM_MF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_nonSyn_Str_last<-p_RM_MF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_Syn<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("MF")&AA_Change=="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,1))+theme(axis.text=element_text(size=15))
p_RM_MF_Syn_fin<-p_RM_MF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_Syn_last<-p_RM_MF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_MF_nonSyn_nsp_legend<-get_legend(p_RM_MF_nonSyn_nsp_fin)
p_RM_MF_nonSyn_Str_legend<-get_legend(p_RM_MF_nonSyn_Str_fin)
p_RM_MF_Syn_legend<-get_legend(p_RM_MF_Syn_fin)
p_RM_MF<-ggarrange(p_RM_MF_nonSyn_nsp_last,p_RM_MF_nonSyn_Str_last,p_RM_MF_Syn_last,p_RM_MF_nonSyn_nsp_legend,p_RM_MF_nonSyn_Str_legend,p_RM_MF_Syn_legend,heights=c(2,1),ncol=3,nrow=2)+draw_plot_label(label=c("a","b","c"),size=20,x=c(0,0.3333,0.6666),y=c(1.05,1.05,1.05))

p_RM_LF_nonSyn_nsp<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LF")&grepl("nsp",Location)&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,0.2))+theme(axis.text=element_text(size=15))
p_RM_LF_nonSyn_nsp_fin<-p_RM_LF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_nonSyn_nsp_last<-p_RM_LF_nonSyn_nsp+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_nonSyn_Str<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LF")&!(grepl("nsp",Location))&!(grepl("UTR",Location))&AA_Change!="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,0.2))+theme(axis.text=element_text(size=15))
p_RM_LF_nonSyn_Str_fin<-p_RM_LF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_nonSyn_Str_last<-p_RM_LF_nonSyn_Str+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=4,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_Syn<-ggplot(data=subset(t_RM_world_labels_v2,y_m!="2021-Apr"&Type%in%c("LF")&AA_Change=="Syn"),aes(x=y_m,y=NRFp,color=label,group=label))+geom_line()+scale_color_manual(values=cols)+scale_y_continuous(limits=c(0,0.2))+theme(axis.text=element_text(size=15))
p_RM_LF_Syn_fin<-p_RM_LF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="bottom",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_Syn_last<-p_RM_LF_Syn+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+theme(legend.position="none",axis.text.x=element_text(angle=90,hjust=0.5,vjust=0.5),axis.title.x=element_blank(),legend.title=element_blank(),legend.text=element_text(size=13))+guides(color=guide_legend(ncol=3,byrow=TRUE))+labs(y="Normalized Relative Frequency (NRFp)")
p_RM_LF_nonSyn_nsp_legend<-get_legend(p_RM_LF_nonSyn_nsp_fin)
p_RM_LF_nonSyn_Str_legend<-get_legend(p_RM_LF_nonSyn_Str_fin)
p_RM_LF_Syn_legend<-get_legend(p_RM_LF_Syn_fin)
p_RM_LF<-ggarrange(p_RM_LF_nonSyn_nsp_last,p_RM_LF_nonSyn_Str_last,p_RM_LF_Syn_last,p_RM_LF_nonSyn_nsp_legend,p_RM_LF_nonSyn_Str_legend,p_RM_LF_Syn_legend,heights=c(2,1),ncol=3,nrow=2)+draw_plot_label(label=c("d","e","f"),size=20,x=c(0,0.3333,0.6666),y=c(1.05,1.05,1.05))

p_RM_MF_LF<-ggarrange(titles,p_RM_MF,p_RM_LF,nrow=3,heights=c(0.7,6,6))
ggsave("Figure_S3.jpg",plot=p_RM_MF_LF,dpi=500,width=30,height=15)
#Analysis of RM mutations by type by continent##########################################
##
#Add type of mutations to the table by continents
t_RM_cont_labels_v2$Type<-"any"
for (i in 1:nrow(t_RM_cont_labels_v2)){
  if (t_RM_cont_labels_v2$code[i] %in% t_RM_world_HF_list_v2){
    t_RM_cont_labels_v2$Type[i]<-"HF"
  }else if(t_RM_cont_labels_v2$code[i] %in% t_RM_world_LFa_list_v2){
    t_RM_cont_labels_v2$Type[i]<-"LFa"
  }else if(t_RM_cont_labels_v2$code[i] %in% t_RM_world_MF_list_v2){
    t_RM_cont_labels_v2$Type[i]<-"MF"
  }else if(t_RM_cont_labels_v2$code[i] %in% t_RM_world_LF_list_v2){
    t_RM_cont_labels_v2$Type[i]<-"LF"
  }
}
t_RM_cont_labels_v2$label<-factor(t_RM_cont_labels_v2$label,levels=unique(t_RM_cont_labels_v2$label[order(t_RM_cont_labels_v2$Position,t_RM_cont_labels_v2$label)]),ordered=TRUE)
#plot of mutations by type by continent
p_RM_cont_HF_group_1<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="HF"&y_m!="2021-Apr"&code%in%c("C241T","C3037T","C14408T","A23403G")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,nrow=1)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Group 1")
p_RM_cont_HF_group_2<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="HF"&y_m!="2021-Apr"&!(code%in%c("C241T","C3037T","C14408T","A23403G"))),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=7)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Group 2 - subtype 1")
p_RM_cont_LFa_group_3<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LFa"&y_m!="2021-Apr"&code!="A28095T"),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=7)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Group 2 - subtype 2")
p_RM_cont_LFa_group_4<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LFa"&y_m!="2021-Apr"&code=="A28095T"),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=7)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Group 2 - subtype 3")
p_RM_cont_HF_12<-ggarrange(p_RM_cont_HF_group_1,p_RM_cont_HF_group_2,ncol=1,heights=c(1,5))
p_RM_cont_HF_34<-ggarrange(p_RM_cont_LFa_group_3,p_RM_cont_LFa_group_4,ncol=1,heights=c(2,1))
ggsave("Figure_S4.jpg",plot=p_RM_cont_HF_12,dpi=500,width=15,height=24)
ggsave("Figure_S5.jpg",plot=p_RM_cont_HF_34,dpi=500,width=15,height=17.15)

p_RM_cont_MF_Europe<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("T445C","C6286T","G21255C","C22227T","C26801G","C28932T","G29645T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=5)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Europe - subtype 3")
t_RM_cont_labels_v2_t<-t_RM_cont_labels_v2
t_RM_cont_labels_v2_t$label<-factor(t_RM_cont_labels_v2_t$label,levels=c("C10319T_L89F_nsp5","A18424G_N129D_nsp14","C21304T_Syn_nsp16","G25907T_G172V_ORF3a","C27964T_S24L_ORF8","C28472T_P67S_N","C28869T_P199L_N","C1059T_T85I_nsp2"),ordered=TRUE)
p_RM_cont_MF_NorthAmerica<-ggplot(data=subset(t_RM_cont_labels_v2_t,Type=="MF"&y_m!="2021-Apr"&code%in%c("C1059T","C10319T","A18424G","C21304T","G25907T","C27964T","C28472T","C28869T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=5)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("North America - subtypes 4-5")
p_RM_cont_MF_Africa<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("C14120T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=5)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Africa - subtype 6")
p_RM_cont_MF_SouthAmerica<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("G25088T","T27299C","T29148C")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=5)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("South America - subtypes 1-2")
t_RM_cont_labels_v2_t<-t_RM_cont_labels_v2
t_RM_cont_labels_v2_t$label<-factor(t_RM_cont_labels_v2_t$label,levels=c("T22917G_L452R_S","C26681T_F53F_M","G29402T_D377Y_N","C18877T_Syn_nsp14","C26735T_Y71Y_M","C28854T_S194L_N","G25563T_Q57H_ORF3a","G28881A_R203K_N","G28882A_R203R_N","G28883C_G204R_N"),ordered=TRUE)
p_RM_cont_MF_morethanone<-ggplot(data=subset(t_RM_cont_labels_v2_t,Type=="MF"&y_m!="2021-Apr"&code%in%c("C18877T","C26735T","C28854T","T22917G","C26681T","G29402T","G25563T","G28881A","G28882A","G28883C")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label,ncol=5)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("More than one region - subtypes 7-10")
p_RM_cont_MF<-ggarrange(p_RM_cont_MF_SouthAmerica,p_RM_cont_MF_Europe,p_RM_cont_MF_NorthAmerica,p_RM_cont_MF_Africa,p_RM_cont_MF_morethanone,ncol=1,heights=c(1,2,2,1,2))
ggsave("Figure_S6.jpg",plot=p_RM_cont_MF,dpi=500,width=15,height=24)

p_RM_cont_MF_group_1<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("T445C","C6286T","G21255C","C22227T","C26801G","C28932T","G29645T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_2<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("C1059T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_3<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("C10319T","A18424G","C21304T","G25907T","C27964T","C28472T","C28869T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_4<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("C14120T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_5<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("G25088T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_6<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("T27299C","T29148C")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_7<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("C18877T","C26735T","C28854T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_8<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("T22917G","G29402T","C26681T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_9<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("G25563T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))
p_RM_cont_MF_group_10<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="MF"&y_m!="2021-Apr"&code%in%c("G28881A","G28882A","G28883C")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))

p_RM_cont_LF_Europe<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LF"&y_m!="2021-Apr"&code%in%c("C27944T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Europe - subtype 3")
p_RM_cont_LF_Asia<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LF"&y_m!="2021-Apr"&code%in%c("C313T","C22444T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Asia - subtypes 5-6")
p_RM_cont_LF_NorthAmerica<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LF"&y_m!="2021-Apr"&code%in%c("G8083A","C14805T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("North America - subtype 4")
p_RM_cont_LF_SouthAmerica<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LF"&y_m!="2021-Apr"&code%in%c("C12053T","C28253T","G28975T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("South America - subtype 1-2")
p_RM_cont_LF_Oceania<-ggplot(data=subset(t_RM_cont_labels_v2,Type=="LF"&y_m!="2021-Apr"&code%in%c("G22992A")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("Oceania - subtype 7")
t_RM_cont_labels_v2_t<-t_RM_cont_labels_v2
t_RM_cont_labels_v2_t$label<-factor(t_RM_cont_labels_v2_t$label,levels=c("C4543T_Syn_nsp3","G15766T_V776L_nsp12","G11083T_L37F_nsp6","A17615G_K460R_nsp13","T19839C_Syn_nsp15","A20268G_Syn_nsp15","C28887T_T205I_N"),ordered=TRUE)
p_RM_cont_LF_morethanone<-ggplot(data=subset(t_RM_cont_labels_v2_t,Type=="LF"&y_m!="2021-Apr"&code%in%c("C4543T","G11083T","G15766T","A17615G","T19839C","A20268G","C25710","C28887T")),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+facet_wrap(~label)+theme(strip.background=element_rect(color="white",fill="white"),legend.position="bottom",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0,size=10),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle("More than one region - subtypes 8-13")
p_RM_cont_LF<-ggarrange(p_RM_cont_LF_SouthAmerica,p_RM_cont_LF_Europe,p_RM_cont_LF_NorthAmerica,p_RM_cont_LF_Asia,p_RM_cont_LF_Oceania,p_RM_cont_LF_morethanone,ncol=1,heights=c(1,1,1,1,1,3))
ggsave("Figure_S13.jpg",plot=p_RM_cont_LF,dpi=500,width=15,height=24)
#Function to determine the presence of type of mutations in the continents
ver_pres_muts<-function(type="HF"){
  data_filt<-t_RM_cont_labels_v2%>%dplyr::filter(Type==type)
  pre_muts<-list()
  for (z in levels(as.factor(data_filt$continent))){
    pre_muts[[z]]<-list()
    for (i in 1:length(levels(as.factor(as.character(data_filt$code))))){
      df<-data_filt%>%dplyr::filter(code==data_filt$code[i],continent==z)
      if (max(df$NRFp)>0){
        pre_muts[[z]][[data_filt$code[i]]]<-data.frame(mutation=data_filt$code[i],Presence="Yes",continent=z)
      }else{
        pre_muts[[z]][[data_filt$code[i]]]<-data.frame(mutation=data_filt$code[i],Presence="No",continent=z)
      }  
    }
  }
  pre_muts_total<-do.call(rbind,pre_muts)
  pre_muts_total_total<-do.call(rbind,pre_muts_total)
  return(pre_muts_total_total)
}
pre_muts_HF<-ver_pres_muts(type="HF")
pre_muts_LFa<-ver_pres_muts(type="LFa")
pre_muts_MF<-ver_pres_muts(type="MF")
pre_muts_LF<-ver_pres_muts(type="LF")
#add type of mutations to the table of RM mutations
RM_labels$Type<-"any"
for (i in 1:nrow(RM_labels)){
  if (RM_labels$code[i] %in% t_RM_world_HF_list_v2){
    RM_labels$Type[i]<-"HF"
  }else if(RM_labels$code[i] %in% t_RM_world_LFa_list_v2){
    RM_labels$Type[i]<-"LFa"
  }else if(RM_labels$code[i] %in% t_RM_world_MF_list_v2){
    RM_labels$Type[i]<-"MF"
  }else if(RM_labels$code[i] %in% t_RM_world_LF_list_v2){
    RM_labels$Type[i]<-"LF"
  }
}
#Calculate the cumulative cases (total number of cases) of RM mutations by continent
t_RM_cont_total<-t_RM_analysis_comp%>%
  dplyr::filter(y_m!="2021-Apr")%>%
  dplyr::group_by(continent,code)%>%
  dplyr::summarise(NRFp=wtd.mean(freq,count_cases),NRFp_sd=wtd.var(freq,count_cases),total_cases=sum(count_cases))%>%
  dplyr::left_join(RM_labels)
t_RM_cont_total_2<-t_RM_cont_total%>%
  dplyr::mutate(p=1-NRFp)
t_RM_cont_total_3<-t_RM_cont_total_2%>%
  dplyr::mutate(p_cases=p*total_cases,NRFp_cases=NRFp*total_cases)
t_RM_cont_total_3$label<-factor(t_RM_cont_total_3$label,levels=unique(t_RM_cont_total_3$label[order(t_RM_cont_total_3$Position,t_RM_cont_total_3$label)]),ordered=TRUE)
t_RM_cont_total_3$code<-factor(t_RM_cont_total_3$code,levels=unique(t_RM_cont_total_3$code[order(t_RM_cont_total_3$Position,t_RM_cont_total_3$code)]),ordered=TRUE)
cont_tables<-list()
for (i in levels(t_RM_cont_total_3$code)){
  print(i)
  df<-t_RM_cont_total_3%>%dplyr::filter(code==i)
  cont_tables[[i]]<-data.frame(mutation=c(i,substr(i,1,nchar(i)-1)))
  row.names(cont_tables[[i]])<-c(i,substr(i,1,nchar(i)-1))
  for (z in 1:nrow(df)){
    cont_tables[[i]][[df$continent[z]]]<-c(round(df$NRFp_cases[z],0),round(df$p_cases[z],0))
  }
  cont_tables[[i]]$mutation<-NULL
}
#chi and pearson residual analysis and plot
#MF mutations
MF_sample<-c("G25088T","T27299C","C22227T","C28472T","C1059T","C14120T","T22917G","C18877T","G25563T","G28883C")
plot_chi_MF_t<-list()
for(z in MF_sample){
  for(i in 1:length(cont_tables)){
    if(names(cont_tables[i])==z&&names(cont_tables[i])!="G28883C"){
      p1<-ggplot(data=subset(t_RM_cont_labels_v2,y_m!="2021-Apr"&code==names(cont_tables[i])),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+theme(legend.position="none",legend.title=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle(t_RM_cont_labels_v2$label[t_RM_cont_labels_v2$code==names(cont_tables[i])][1])
      p2<-ggplot(data=subset(t_RM_cont_total_3,code==names(cont_tables[i])),aes(x=continent,y=NRFp,fill=continent))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2,alpha=0.5)+scale_fill_manual(values=col_conts)+theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+geom_text(aes(label=formatC(round(NRFp_cases,0),format="e",digits=1)),position=position_dodge(width=0.9),vjust="inward")
      chisq<-chisq.test(cont_tables[[i]])
      contrib <- 100*chisq$residuals^2/chisq$statistic
      contrib_df<-as.data.frame(chisq$residuals)
      contrib_df$code<-row.names(contrib_df)
      t_contrib_df<-melt(contrib_df,id.vars="code")
      p3<-ggplot(t_contrib_df, aes(x=variable,y=code,fill=value))+geom_tile()+theme(axis.title=element_blank())+scale_fill_gradient2(low = "grey", mid = "white", high = "brown")+theme(axis.text.x=element_blank())+labs(fill="Residuals")+ggtitle(ifelse(chisq$p.value==0,"p < 0.01",chisq$p.value))
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_1",sep="")]]<-p1
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_2",sep="")]]<-p2
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_3",sep="")]]<-p3
    }else if(names(cont_tables[i])==z&&names(cont_tables[i])=="G28883C"){
      p1<-ggplot(data=subset(t_RM_cont_labels_v2,y_m!="2021-Apr"&code==names(cont_tables[i])),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+theme(legend.position="none",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle(t_RM_cont_labels_v2$label[t_RM_cont_labels_v2$code==names(cont_tables[i])][1])
      p2<-ggplot(data=subset(t_RM_cont_total_3,code==names(cont_tables[i])),aes(x=continent,y=NRFp,fill=continent))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2,alpha=0.5)+scale_fill_manual(values=col_conts)+theme(legend.position="none",axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5),axis.title.x=element_blank())+geom_text(aes(label=formatC(round(NRFp_cases,0),format="e",digits=1)),position=position_dodge(width=0.9),vjust="inward")
      chisq<-chisq.test(cont_tables[[i]])
      contrib <- 100*chisq$residuals^2/chisq$statistic
      contrib_df<-as.data.frame(chisq$residuals)
      contrib_df$code<-row.names(contrib_df)
      t_contrib_df<-melt(contrib_df,id.vars="code")
      p3<-ggplot(t_contrib_df, aes(x=variable,y=code,fill=value))+geom_tile()+theme(axis.title=element_blank())+scale_fill_gradient2(low = "grey", mid = "white", high = "brown")+theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5))+labs(fill="Residuals")+ggtitle(ifelse(chisq$p.value==0,"p < 0.01",chisq$p.value))
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_1",sep="")]]<-p1
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_2",sep="")]]<-p2
      plot_chi_MF_t[[paste(names(cont_tables[i]),"_3",sep="")]]<-p3
    }
  }
}
ptotal_MF<-egg::ggarrange(plots=plot_chi_MF_t,nrow=10,labels=c("a","","","b","","","c","","","d","","","e","","","f","","","g","","","h","","","i","","","j","",""),label.args=list(gp=grid::gpar(fontsize=20,fontface="bold"),hjust=-0.5))
ggsave("Figure_2.jpg",plot=ptotal_MF,dpi=500,width=20,height=20)
#LF mutations
LF_sample<-c("C12053T","C28253T","C27944T","G8083A","C313T","C22444T","G22992A","C4543T","G11083T","A17615G","T19839C","A20268G","C28887T")
plot_chi_LF_t<-list()
for(z in LF_sample){
  for(i in 1:length(cont_tables)){
    if(names(cont_tables[i])==z&&names(cont_tables[i])!="C28887T"){
      p1<-ggplot(data=subset(t_RM_cont_labels_v2,y_m!="2021-Apr"&code==names(cont_tables[i])),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+theme(legend.position="none",legend.title=element_blank(),axis.text.x=element_blank(),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle(t_RM_cont_labels_v2$label[t_RM_cont_labels_v2$code==names(cont_tables[i])][1])
      p2<-ggplot(data=subset(t_RM_cont_total_3,code==names(cont_tables[i])),aes(x=continent,y=NRFp,fill=continent))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2,alpha=0.5)+scale_fill_manual(values=col_conts)+theme(legend.position="none",axis.text.x=element_blank(),axis.title.x=element_blank())+geom_text(aes(label=formatC(round(NRFp_cases,0),format="e",digits=1)),position=position_dodge(width=0.9),vjust="inward")
      chisq<-chisq.test(cont_tables[[i]])
      contrib <- 100*chisq$residuals^2/chisq$statistic
      contrib_df<-as.data.frame(chisq$residuals)
      contrib_df$code<-row.names(contrib_df)
      t_contrib_df<-melt(contrib_df,id.vars="code")
      p3<-ggplot(t_contrib_df, aes(x=variable,y=code,fill=value))+geom_tile()+theme(axis.title=element_blank())+scale_fill_gradient2(low = "grey", mid = "white", high = "brown")+theme(axis.text.x=element_blank())+labs(fill="Residuals")+ggtitle(ifelse(chisq$p.value==0,"p < 0.01",chisq$p.value))
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_1",sep="")]]<-p1
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_2",sep="")]]<-p2
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_3",sep="")]]<-p3
    }else if(names(cont_tables[i])==z&&names(cont_tables[i])=="C28887T"){
      p1<-ggplot(data=subset(t_RM_cont_labels_v2,y_m!="2021-Apr"&code==names(cont_tables[i])),aes(x=y_m,y=NRFp,color=continent,group=continent))+geom_line(size=1)+theme(legend.position="none",legend.title=element_blank(),axis.text.x=element_text(angle=90,hjust=1,vjust=0),axis.title.x=element_blank())+labs(y="NRFp")+scale_color_manual(values=col_conts)+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2)+guides(color=guide_legend(nrow=1))+ggtitle(t_RM_cont_labels_v2$label[t_RM_cont_labels_v2$code==names(cont_tables[i])][1])
      p2<-ggplot(data=subset(t_RM_cont_total_3,code==names(cont_tables[i])),aes(x=continent,y=NRFp,fill=continent))+geom_bar(stat="identity")+geom_errorbar(aes(ymin=pmax(0,NRFp-NRFp_sd),ymax=pmin(NRFp+NRFp_sd,1)),show.legend=FALSE,width=.2,alpha=0.5)+scale_fill_manual(values=col_conts)+theme(legend.position="none",axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5),axis.title.x=element_blank())+geom_text(aes(label=formatC(round(NRFp_cases,0),format="e",digits=1)),position=position_dodge(width=0.9),vjust="inward")
      chisq<-chisq.test(cont_tables[[i]])
      contrib <- 100*chisq$residuals^2/chisq$statistic
      contrib_df<-as.data.frame(chisq$residuals)
      contrib_df$code<-row.names(contrib_df)
      t_contrib_df<-melt(contrib_df,id.vars="code")
      p3<-ggplot(t_contrib_df, aes(x=variable,y=code,fill=value))+geom_tile()+theme(axis.title=element_blank())+scale_fill_gradient2(low = "grey", mid = "white", high = "brown")+theme(axis.text.x=element_text(angle=45,vjust=0.5,hjust=0.5))+labs(fill="Residuals")+ggtitle(ifelse(chisq$p.value==0,"p < 0.01",chisq$p.value))
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_1",sep="")]]<-p1
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_2",sep="")]]<-p2
      plot_chi_LF_t[[paste(names(cont_tables[i]),"_3",sep="")]]<-p3
    }
  }
}
ptotal_LF<-egg::ggarrange(plots=plot_chi_LF_t,nrow=13,labels=c("a","","","b","","","c","","","d","","","e","","","f","","","g","","","h","","","i","","","j","","","k","","","l","","","m","",""),label.args=list(gp=grid::gpar(fontsize=20,fontface="bold"),hjust=-0.5))
ggsave("Figure_S14.jpg",plot=ptotal_LF,dpi=500,width=20,height=20)
#analysis of the MF mutations specific of regions by country
t_RM_count_sample_MF<-t_RM_analysis_comp%>%
  dplyr::filter(code%in%MF_sample)%>%
  dplyr::left_join(RM_labels)%>%
  dplyr::ungroup()%>%
  dplyr::group_by(code,country)%>%
  dplyr::mutate(cum_sum=cumsum(value))
t_RM_count_sample_MF_cum<-t_RM_count_sample_MF%>%
  dplyr::ungroup()%>%
  tidyr::complete(y_m,nesting(code,continent,country,label,Position,AA_Change,Location,Type),fill=list(value=0,count_genomes=0,count_cases=0,freq=0,CI=0))%>%
  dplyr::group_by(code,country)%>%
  dplyr::mutate(cum_sum=cumsum(value))
p_fr_coun<-function(code_t,continent_t){
  df<-t_RM_count_sample_MF%>%
    dplyr::filter(code==code_t,continent==continent_t)
  p_count_MF_fr<-ggplot(data=subset(df,y_m!="2021-Apr"&cum_sum>0),aes(x=y_m,y=freq,color=continent))+geom_point()+geom_errorbar(aes(ymin=pmax(0,freq-CI),ymax=pmin(freq+CI,1)),width=0.2)+facet_wrap(~country,ncol=1)+theme(strip.background=element_rect(fill="white",color="white"),legend.position="none",axis.title=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1))+ggtitle(df$label[1])+scale_color_manual(values=col_conts)
  return(p_count_MF_fr)
}
p_cs_coun<-function(code_t,continent_t){
  df<-t_RM_count_sample_MF_cum%>%
  dplyr::filter(code==code_t,continent==continent_t)
  p_count_MF_cs<-ggplot(data=subset(df,y_m!="2021-Apr"&cum_sum>0),aes(x=y_m,y=cum_sum))+geom_bar(aes(fill=continent),stat="identity")+facet_wrap(~country,scales="free_y",ncol=1)+theme(strip.background=element_rect(fill="white",color="white"),legend.position="none",axis.title=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1))+ggtitle(df$label[1])+scale_fill_manual(values=col_conts)
  return(p_count_MF_cs)
}
p_fr_list<-list()
p_cs_list<-list()
p_fr_list[["G25088T"]]<-p_fr_coun("G25088T","SouthAmerica")
p_fr_list[["T27299C"]]<-p_fr_coun("T27299C","SouthAmerica")
p_fr_list[["C28472T"]]<-p_fr_coun("C28472T","NorthAmerica")
p_fr_list[["C1059T"]]<-p_fr_coun("C1059T","NorthAmerica")
p_cs_list[["G25088T"]]<-p_cs_coun("G25088T","SouthAmerica")
p_cs_list[["T27299C"]]<-p_cs_coun("T27299C","SouthAmerica")
p_cs_list[["C28472T"]]<-p_cs_coun("C28472T","NorthAmerica")
p_cs_list[["C1059T"]]<-p_cs_coun("C1059T","NorthAmerica")

p_fr<-ggarrange(plotlist=p_fr_list,labels=letters,ncol=4)  
p_fr_f<-annotate_figure(p_fr,left=text_grob("Relative Frequencies",rot=90,size=15))
ggsave("Figure_S7.jpg",plot=p_fr_f,dpi=500,height=10,width=15)
p_cs<-ggarrange(plotlist=p_cs_list,labels=letters,ncol=4)  
p_cs_f<-annotate_figure(p_cs,left=text_grob("Cumulative Number of Cases",rot=90,size=15))
ggsave("Figure_S9.jpg",plot=p_cs_f,dpi=500,height=10,width=15)

df<-t_RM_count_sample_MF%>%
  dplyr::filter(code=="C22227T",continent=="Europe")
p_count_MF_fr<-ggplot(data=subset(df,y_m!="2021-Apr"&cum_sum>0),aes(x=y_m,y=freq,color=continent))+geom_point()+geom_errorbar(aes(ymin=pmax(0,freq-CI),ymax=pmin(freq+CI,1)),width=0.2)+facet_wrap(~country)+theme(strip.background=element_rect(fill="white",color="white"),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=0,vjust=0))+ggtitle(df$label[1])+scale_color_manual(values=col_conts)+labs(y="Relative Frequencies")
ggsave("Figure_S8.jpg",plot=p_count_MF_fr,dpi=500,height=10,width=15)

df<-t_RM_count_sample_MF_cum%>%
  dplyr::filter(code=="C22227T",continent=="Europe")
p_count_MF_cs<-ggplot(data=subset(df,y_m!="2021-Apr"&cum_sum>0),aes(x=y_m,y=cum_sum))+geom_bar(aes(fill=continent),stat="identity")+facet_wrap(~country,scales="free_y")+theme(strip.background=element_rect(fill="white",color="white"),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_text(angle=90,hjust=0,vjust=0,size=10))+ggtitle(df$label[1])+scale_fill_manual(values=col_conts)+labs(y="Cumulative Number of Cases")
ggsave("Figure_S10.jpg",plot=p_count_MF_cs,dpi=500,height=10,width=15)

#Generating table 1 of the paper##########################################
##
#preparation
RM_annotations$Type<-"test"
for (i in 1:nrow(RM_annotations)){
  if (RM_annotations$code[i] %in% t_RM_world_HF_list_v2){
    RM_annotations$Type[i]<-"HF"
  }else if(RM_annotations$code[i] %in% t_RM_world_LFa_list_v2){
    RM_annotations$Type[i]<-"HF"
  }else if(RM_annotations$code[i] %in% t_RM_world_MF_list_v2){
    RM_annotations$Type[i]<-"MF"
  }else if(RM_labels$code[i] %in% t_RM_world_LF_list_v2){
    RM_annotations$Type[i]<-"LF"
  }
}
RM_table_paper<-RM_annotations%>%
  dplyr::filter(!grepl("-_---_5-UTR",label))%>%
  dplyr::select(ID=label,NRFp=NRFp,Type,nt_change=code,AA_change=AA_Change,Location,wt_codon=ref_codon,mut_codon)
write.xlsx(RM_table_paper,"Table_1_paper.xlsx")
#Creating the table of stringency index information##########################################
##
#information of stringency index by country
#obtained from: https://ourworldindata.org/covid-stringency-index
read_pols<-function(file){
  str_ind<-read.csv(file)
  str_ind$Entity<-gsub("United States","USA",str_ind$Entity)
  str_ind$Entity<-gsub("United Kingdom","UnitedKingdom",str_ind$Entity)
  str_ind$Entity<-gsub("Democratic Republic of Congo","DemocraticCongo",str_ind$Entity)
  str_ind_R<-str_ind%>%
    dplyr::filter(Entity%in%sel_countries)%>%
    dplyr::select(country=Entity,date=Day,4)
  str_ind_R$date<-as.Date(str_ind_R$date)
  names(str_ind_R)[3]<-paste(word(sapply(strsplit(file,"-"),"[",1),start=1,end=4,sep="|"),word(sapply(strsplit(file,"-"),"[",2),start=1,end=4,sep="|"),sep="_")
  return(str_ind_R)
}
str_ind_R<-read_pols("covid-stringency-index.csv")
#Creating the table of cases by day to estimate the Rt per mutation##########################################
#rawdata of cases by day
#downloaded from here https://ourworldindata.org/covid-cases
cas_wor<-read.csv("owid-covid-data.csv")
cas_wor$location<-gsub(" ","",cas_wor$location)
cas_wor$location<-gsub("UnitedStates","USA",cas_wor$location)
cas_filt<-cas_wor%>%
  dplyr::filter(location%in%countries_diff_freqs)%>%
  dplyr::select(country=location,date,total_cases,new_cases)
cas_filt$year<-sapply(strsplit(cas_filt$date,"-"),"[",1)
cas_filt$month<-sapply(strsplit(cas_filt$date,"-"),"[",2)
cas_filt$month<-as.numeric(cas_filt$month)
cas_filt$Month<-as.character(lubridate::month(ymd(010101) + months(cas_filt$month-1),label=TRUE,abbr=TRUE))
cas_filt$y_m<-paste(cas_filt$year,cas_filt$Month,sep="-")
cas_count<-cas_filt%>%
  dplyr::filter(y_m%in%c("2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar"))%>%
  dplyr::select(date=date,confirm=new_cases,country)
#creating the column of week
cas_count$Date<-as.Date(cas_count$date,format="%Y-%m-%d")
rep_week<-cas_count%>%
  dplyr::mutate(week=cut.Date(Date,breaks="1 week",labels=FALSE))%>%
  dplyr::select(country,Date,confirm,week)
rep_week$confirm[is.na(rep_week$confirm)] <- 0
#Calculating the relative frequencies of mutations of interest by week of countries with sufficient data##########################################
##
#function to verify which countries has more than 14 genomes in the months of interest
test_genomes<-function(pais){
  df_r<-hq_data_analysis%>%
    dplyr::filter(country==pais & y_m%in%c("2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar"))
  df_r$colect<-as.Date(df_r$date,format="%Y-%m-%d")
  df_new<-df_r%>%
    dplyr::mutate(week=cut.Date(colect,breaks="1 week",labels=FALSE))%>%
    dplyr::arrange(colect)
  df_2<-df_new%>%  
    dplyr::group_by(week)%>%
    dplyr::summarise(count=n())
  df_3<-df_2%>%
    dplyr::filter(count<10)
  p_aus<-ggplot(df_3,aes(x=week,y=count))+geom_bar(stat="identity")
  ggplotly(p_aus)
}
test_genomes("USA")
#after test if of the sixteen countries those in "sel_countries" were the selected
sel_countries<-c("Australia","Canada","Germany","India","Japan","Netherlands","Switzerland","UnitedKingdom","USA")
#HF
#creation of lists (DG_NY) of the genomes that contain the mutations D614G or N501Y
#these lists were obtained using python scripts and are described on the folder preparations
D614G<-scan(file="aligned_290Ns_minusbad_align_clear_analysis_ids_seqs_D614G.list",what="",sep=",")
N501Y<-scan(file="aligned_290Ns_minusbad_align_clear_analysis_ids_seqs_N501Y.list",what="",sep=",")
DG_NY<-c(intersect(D614G,N501Y),setdiff(D614G,N501Y),setdiff(N501Y,D614G))
HFs<-data.frame(strain=DG_NY,HF="yes")
DG<-data.frame(strain=D614G,DG="yes")
NY<-data.frame(strain=N501Y,NY="yes")
#MF
#creation of lists (DG_NY) of the genomes that contain the mutations D614G or N501Y
#this list was obtained using python scripts and are described on the folder preparations
R203K<-scan(file="aligned_290Ns_minusbad_align_clear_analysis_ids_seqs_R203K.list",what="",sep=",")
MFs<-data.frame(strain=R203K,MF="yes")
#creating the table that contain information of what mutation contains each genome
df_r_MHF<-hq_data_analysis%>%
  dplyr::filter(country%in%countries_diff_freqs & y_m%in%c("2020-Feb","2020-Mar","2020-Apr","2020-May","2020-Jun","2020-Jul","2020-Aug","2020-Sep","2020-Oct","2020-Nov","2020-Dec","2021-Jan","2021-Feb","2021-Mar"))%>%
  dplyr::select(strain,date,country,continent)
df_r_MHF$Date<-as.Date(df_r_MHF$date,format="%Y-%m-%d")
df_r_HF_id<-dplyr::left_join(df_r_MHF,HFs)
df_r_HF_id$HF[is.na(df_r_HF_id$HF)] <- "no"
df_r_MHF_id<-dplyr::left_join(df_r_HF_id,MFs)
df_r_MHF_id$MF[is.na(df_r_MHF_id$MF)] <- "no"
df_r_MHFD_id<-dplyr::left_join(df_r_MHF_id,DG)
df_r_MHFD_id$DG[is.na(df_r_MHFD_id$DG)] <- "no"
df_r_MHFDN_id<-dplyr::left_join(df_r_MHFD_id,NY)
df_r_MHFDN_id$NY[is.na(df_r_MHFDN_id$NY)] <- "no"

df_r_MHFDN_week<-df_r_MHFDN_id%>%
  dplyr::mutate(week=cut.Date(Date,breaks="1 week",labels=FALSE))
#calculating the relative frequencies of genomes by week that contains D614G or N501Y
df_HF_tables<-list()
for (i in sel_countries){
  tab<-df_r_MHFDN_week%>%
    dplyr::filter(country==i)%>%
    tidyr::drop_na(week)%>%
    dplyr::group_by(week,HF)%>%
    dplyr::summarise(c_HF=n())%>%
    dplyr::mutate(p_HF=c_HF/sum(c_HF))
  tab$country<-i
  df_HF_tables[[i]]<-tab
}
df_HF_grouped<-do.call(rbind,df_HF_tables)
df_HF_group_com<-df_HF_grouped%>%
  tidyr::complete(HF,nesting(country,week),fill=list(c_HF=0,p_HF=0))
df_HF_group_com_f<-df_HF_group_com%>%
  dplyr::ungroup()%>%
  tidyr::complete(week,nesting(HF,country),fill=list(p_HF=0,c_HF=0))
df_HF_group_com_fin<-df_HF_group_com_f%>%
  dplyr::group_by(country,week)%>%
  dplyr::mutate(tot_gen_week=sum(c_HF))
df_HF_group_com_fin$CI<-1.96*(sqrt(((1-df_HF_group_com_fin$p_HF)*(df_HF_group_com_fin$p_HF))/df_HF_group_com_fin$tot_gen_week))
df_HF_group_com_fin$CI[is.nan(df_HF_group_com_fin$CI)]<-0
#calculating the relative frequencies of genomes by week that contains R203K
df_MF_tables<-list()
for (i in sel_countries){
  tab<-df_r_MHFDN_week%>%
    dplyr::filter(country==i)%>%
    tidyr::drop_na(week)%>%
    dplyr::group_by(week,MF)%>%
    dplyr::summarise(c_MF=n())%>%
    dplyr::mutate(p_MF=c_MF/sum(c_MF))
  tab$country<-i
  df_MF_tables[[i]]<-tab
}
df_MF_grouped<-do.call(rbind,df_MF_tables)
df_MF_group_com<-df_MF_grouped%>%
  tidyr::complete(MF,nesting(country,week),fill=list(c_MF=0,p_MF=0))
df_MF_group_com_f<-df_MF_group_com%>%
  dplyr::ungroup()%>%
  tidyr::complete(week,nesting(MF,country),fill=list(p_MF=0,c_MF=0))
df_MF_group_com_fin<-df_MF_group_com_f%>%
  dplyr::group_by(country,week)%>%
  dplyr::mutate(tot_gen_week=sum(c_MF))
df_MF_group_com_fin$CI<-1.96*(sqrt(((1-df_MF_group_com_fin$p_MF)*(df_MF_group_com_fin$p_MF))/df_MF_group_com_fin$tot_gen_week))
df_MF_group_com_fin$CI[is.nan(df_MF_group_com_fin$CI)]<-0
#calculating the relative frequencies of genomes by week that contains D614G
df_DG_tables<-list()
for (i in sel_countries){
  tab<-df_r_MHFDN_week%>%
    dplyr::filter(country==i)%>%
    tidyr::drop_na(week)%>%
    dplyr::group_by(week,DG)%>%
    dplyr::summarise(c_DG=n())%>%
    dplyr::mutate(p_DG=c_DG/sum(c_DG))
  tab$country<-i
  df_DG_tables[[i]]<-tab
}
df_DG_grouped<-do.call(rbind,df_DG_tables)
df_DG_group_com<-df_DG_grouped%>%
  tidyr::complete(DG,nesting(country,week),fill=list(c_DG=0,p_DG=0))
df_DG_group_com_f<-df_DG_group_com%>%
  dplyr::ungroup()%>%
  tidyr::complete(week,nesting(DG,country),fill=list(p_DG=0,c_DG=0))
df_DG_group_com_fin<-df_DG_group_com_f%>%
  dplyr::group_by(country,week)%>%
  dplyr::mutate(tot_gen_week=sum(c_DG))
df_DG_group_com_fin$CI<-1.96*(sqrt(((1-df_DG_group_com_fin$p_DG)*(df_DG_group_com_fin$p_DG))/df_DG_group_com_fin$tot_gen_week))
df_DG_group_com_fin$CI[is.nan(df_DG_group_com_fin$CI)]<-0
#calculating the relative frequencies of genomes by week that contains N501Y
df_NY_tables<-list()
for (i in sel_countries){
  tab<-df_r_MHFDN_week%>%
    dplyr::filter(country==i)%>%
    tidyr::drop_na(week)%>%
    dplyr::group_by(week,NY)%>%
    dplyr::summarise(c_NY=n())%>%
    dplyr::mutate(p_NY=c_NY/sum(c_NY))
  tab$country<-i
  df_NY_tables[[i]]<-tab
}
df_NY_grouped<-do.call(rbind,df_NY_tables)
df_NY_group_com<-df_NY_grouped%>%
  tidyr::complete(NY,nesting(country,week),fill=list(c_NY=0,p_NY=0))
df_NY_group_com_f<-df_NY_group_com%>%
  dplyr::ungroup()%>%
  tidyr::complete(week,nesting(NY,country),fill=list(p_NY=0,c_NY=0))
df_NY_group_com_fin<-df_NY_group_com_f%>%
  dplyr::group_by(country,week)%>%
  dplyr::mutate(tot_gen_week=sum(c_NY))
df_NY_group_com_fin$CI<-1.96*(sqrt(((1-df_NY_group_com_fin$p_NY)*(df_NY_group_com_fin$p_NY))/df_NY_group_com_fin$tot_gen_week))
df_NY_group_com_fin$CI[is.nan(df_NY_group_com_fin$CI)]<-0
#Creating the table of frequencies of mutations by week and cases by day##########################################
#merge cases with genomes
rep_week_count<-rep_week%>%
  dplyr::filter(country%in%sel_countries)
rep_week_count_com<-rep_week_count%>%
  complete(Date,nesting(country),fill=list(confirm=0))%>%
  tidyr::fill(week)
rep_week_count_com$week<-as.numeric(rep_week_count_com$week)

#Estimation of number of cases (corrected logistically when necessary)##########################################
##
#merge the tables of frequencies of mutations with the table of number of cases
#next we calculated the estimated number of cases with a mutation present
HF_cas<-dplyr::full_join(df_HF_group_com_fin,rep_week_count_com)
HF_cas$HF_confirm<-round(HF_cas$p_HF*HF_cas$confirm,0)
HF_cas_final_roll<-HF_cas%>%
  dplyr::group_by(country,HF)%>%
  dplyr::mutate(roll_mean=zoo::rollmean(confirm,k=7,fill=0),roll_mean_HF=zoo::rollmean(HF_confirm,k=7,fill=0))

MF_cas<-dplyr::full_join(df_MF_group_com_fin,rep_week_count_com)
MF_cas$MF_confirm<-round(MF_cas$p_MF*MF_cas$confirm,0)
MF_cas_final_roll<-MF_cas%>%
  dplyr::group_by(country,MF)%>%
  dplyr::mutate(roll_mean=zoo::rollmean(confirm,k=7,fill=0),roll_mean_MF=zoo::rollmean(MF_confirm,k=7,fill=0))

DG_cas<-dplyr::full_join(df_DG_group_com_fin,rep_week_count_com)
pred_prob_DG_list<-list()
p_pred_prob_DG_list<-list()
for(i in sel_countries){
  test_df<-DG_cas%>%
    dplyr::filter(DG=="yes",country==i)
  model<-glm(p_DG~week,data=test_df,family="binomial")
  pred_DG<-data.frame(p_DG_c=model$fitted.values,week=test_df$week,DG="yes",country=i)
  pred_DG_u<-unique(pred_DG)
  pred_notDG<-data.frame(p_DG_c=1-model$fitted.values,week=test_df$week,DG="no",country=i)
  pred_notDG_u<-unique(pred_notDG)
  pred_prob_DG_list[[i]]<-rbind(pred_DG_u,pred_notDG_u)
  p_pred_prob_DG_list[[i]]<-ggplot(test_df,aes(x=as.POSIXct(as.Date(week*7,origin="2020-02-01")),y=p_DG))+geom_point()+geom_errorbar(aes(ymin=pmax(0,p_DG-CI),ymax=pmin(p_DG+CI,1)))+geom_smooth(method="glm",method.args=list(family="binomial"),color="#1857de",fill="#1857de")+scale_x_datetime(date_breaks="4 weeks",date_labels="%Y-%b")+theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1))+labs(y="Proportion of HF mutation (D614G)")+ggtitle(i)
}
pred_prob_DG<-do.call(rbind,pred_prob_DG_list)
DG_cas_cor<-left_join(DG_cas,pred_prob_DG)
DG_cas_cor$DG_confirm<-round(DG_cas_cor$p_DG_c*DG_cas_cor$confirm,0)
DG_cas_final_roll<-DG_cas_cor%>%
  dplyr::group_by(country,DG)%>%
  dplyr::mutate(roll_mean=zoo::rollmean(confirm,k=7,fill=0),roll_mean_DG=zoo::rollmean(DG_confirm,k=7,fill=0))

NY_cas<-dplyr::full_join(df_NY_group_com_fin,rep_week_count_com)
pred_prob_NY_list<-list()
p_pred_prob_NY_list<-list()
for(i in sel_countries){
  test_df<-NY_cas%>%
    dplyr::filter(NY=="yes",country==i)
  model<-glm(p_NY~week,data=test_df,family="binomial")
  pred_NY<-data.frame(p_NY_c=model$fitted.values,week=test_df$week,NY="yes",country=i)
  pred_NY_u<-unique(pred_NY)
  pred_notNY<-data.frame(p_NY_c=1-model$fitted.values,week=test_df$week,NY="no",country=i)
  pred_notNY_u<-unique(pred_notNY)
  pred_prob_NY_list[[i]]<-rbind(pred_NY_u,pred_notNY_u)
  p_pred_prob_NY_list[[i]]<-ggplot(test_df,aes(x=as.POSIXct(as.Date(week*7,origin="2020-02-01")),y=p_NY))+geom_point()+geom_errorbar(aes(ymin=pmax(0,p_NY-CI),ymax=pmin(p_NY+CI,1)))+geom_smooth(method="glm",method.args=list(family="binomial"),color="#1857de",fill="#1857de")+scale_x_datetime(date_breaks="4 weeks",date_labels="%Y-%b")+theme(axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1))+labs(y="Proportion of HF mutation (N501Y)")+ggtitle(i)
}
pred_prob_NY<-do.call(rbind,pred_prob_NY_list)
NY_cas_cor<-left_join(NY_cas,pred_prob_NY)
NY_cas_cor$NY_confirm<-round(NY_cas_cor$p_NY_c*NY_cas_cor$confirm,0)
NY_cas_final_roll<-NY_cas_cor%>%
  dplyr::group_by(country,NY)%>%
  dplyr::mutate(roll_mean=zoo::rollmean(confirm,k=7,fill=0),roll_mean_NY=zoo::rollmean(NY_confirm,k=7,fill=0))
p_log_count<-ggarrange(plotlist=p_pred_prob_NY_list,labels=letters,ncol=3,nrow=3)
ggsave("Figure_S16.jpg",plot=p_log_count,dpi=500,height=15,width=15)
#create the tables that will be used to estimate the Rt
write.table(HF_cas_final_roll,"HF_estimations.tsv",row.names=FALSE,sep="\t",quote=FALSE)
write.table(MF_cas_final_roll,"MF_estimations.tsv",row.names=FALSE,sep="\t",quote=FALSE)
write.table(DG_cas_final_roll,"DG_estimations.tsv",row.names=FALSE,sep="\t",quote=FALSE)
write.table(NY_cas_final_roll,"NY_estimations.tsv",row.names=FALSE,sep="\t",quote=FALSE)
#Estimation of the Rt by mutations and bind with policies data##########################################
##
#RT_estimate function was used
f_mat_list<-list()
for(i in sel_countries){
  name_HF<-paste("HF",i,"yes",sep="_")
  f_mat_list[[name_HF]]<-RT_estimate("HF_estimations.tsv",i,"yes")
  name_HF<-paste("HF",i,"no",sep="_")
  f_mat_list[[name_HF]]<-RT_estimate("HF_estimations.tsv",i,"no")
  name_MF<-paste("MF",i,"yes",sep="_")
  f_mat_list[[name_MF]]<-RT_estimate("MF_estimations.tsv",i,"yes")
  name_MF<-paste("MF",i,"no",sep="_")
  f_mat_list[[name_MF]]<-RT_estimate("MF_estimations.tsv",i,"no")
  name_DG<-paste("DG",i,"yes",sep="_")
  f_mat_list[[name_DG]]<-RT_estimate("DG_estimations.tsv",i,"yes")
  name_DG<-paste("DG",i,"no",sep="_")
  f_mat_list[[name_DG]]<-RT_estimate("DG_estimations.tsv",i,"no")
  name_NY<-paste("NY",i,"yes",sep="_")
  f_mat_list[[name_NY]]<-RT_estimate("NY_estimations.tsv",i,"yes")
  name_NY<-paste("NY",i,"no",sep="_")
  f_mat_list[[name_NY]]<-RT_estimate("NY_estimations.tsv",i,"no")
}
f_mat<-do.call(rbind,f_mat_list)
#bind Rt data of RM mutations with policies data
R_matrix_final<-Reduce(function(x,y) dplyr::left_join(x,y,all=TRUE),list(f_mat,str_ind_R))
R_matrix_final$Date<-as.character(R_matrix_final$date)
R_matrix_final$year<-sapply(strsplit(R_matrix_final$Date, "-"), "[", 1)
R_matrix_final$month<-sapply(strsplit(R_matrix_final$Date, "-"), "[", 2)
R_matrix_final$month<-as.numeric(R_matrix_final$month)
R_matrix_final$Month<-as.character(lubridate::month(ymd(010101) + months(R_matrix_final$month-1),label=TRUE,abbr=TRUE))
R_matrix_final$y_m<-paste(R_matrix_final$year,R_matrix_final$Month,sep="-")
R_matrix_final$str_fac<-cut(R_matrix_final$cov_str,breaks=seq(0,100,10),labels=c("0","1","2","3","4","5","6","7","8","9"),include.lowest=TRUE)
R_matrix_final$str_fac<-as.factor(R_matrix_final$str_fac)
R_matrix_final$country<-as.factor(R_matrix_final$country)
R_matrix_final$type<-as.factor(R_matrix_final$type)
R_matrix_final$mut<-as.factor(R_matrix_final$mut)
#calculating boot_mean by level of stringency by country##########################################
##
R_matrix_for_boot_1<-R_matrix_final%>%
  dplyr::filter(!(y_m%in%c("2020-Feb","2020-Mar")))
R_matrix_for_boot_list<-list()
matrix_boot_list<-list()
n_iter<-1000
for (c in levels(R_matrix_for_boot_1$country)){
  for (l in levels(R_matrix_for_boot_1$str_fac)){
    for (m in levels(R_matrix_for_boot_1$mut)){
      for (t in levels(R_matrix_for_boot_1$type)){
        dat_yes<-R_matrix_for_boot_1%>%
          dplyr::filter(country==c&str_fac==l&mut=="yes"&type==t)
        dat_no<-R_matrix_for_boot_1%>%
          dplyr::filter(country==c&str_fac==l&mut=="no"&type==t)
        dat<-R_matrix_for_boot_1%>%
          dplyr::filter(country==c&str_fac==l&mut==m&type==t)
        if(nrow(dat_yes)>=25 & nrow(dat_no)>=25){
          dat_new<-dat%>%
            dplyr::arrange(date)%>%
            dplyr::mutate(Rt_14=0)
          for (i in 1:nrow(dat_new)){
            if (i %in% (nrow(dat_new)-13):nrow(dat_new)){
              dat_new$Rt_14[i]<-NA
            }else{
              dat_new$Rt_14[i]<-dat_new$Rt[i+13]
            }
          }
          dat_fin<-dat_new%>%tidyr::drop_na()
          n<-nrow(dat_fin)
          mean<-numeric(n_iter)
          for (iter in 1:n_iter){
            inc<-sample(1:n,replace=TRUE)
            mean[iter]<-mean(dat_fin$Rt_14[inc])
          }
          name<-paste(c,l,m,t,sep="_")
          R_matrix_for_boot_list[[name]]<-dat_fin
          matrix_boot_list[[name]]<-data.frame(country=c,str_fac=l,type=t,mut=m,boot_mean=mean)
        }
      }
    }
  }
}
matrix_boot<-do.call(rbind,matrix_boot_list)
matrix_boot$country<-as.factor(matrix_boot$country)
matrix_boot$str_fac<-as.factor(matrix_boot$str_fac)
matrix_boot$type<-as.factor(matrix_boot$type)
matrix_boot$mut<-as.factor(matrix_boot$mut)
R_matrix_for_boot<-do.call(rbind,R_matrix_for_boot_list)
R_matrix_for_boot$country<-as.factor(R_matrix_for_boot$country)
R_matrix_for_boot$str_fac<-as.factor(R_matrix_for_boot$str_fac)
R_matrix_for_boot$type<-as.factor(R_matrix_for_boot$type)
R_matrix_for_boot$mut<-as.factor(R_matrix_for_boot$mut)

#Calculating significant differences between Rt_14 of different levels of stringency by country##########################################
#MF (R203K)
#Calculating significant difference of MF
matrix_boot_MF<-R_matrix_for_boot%>%
  dplyr::filter(type=="MF")
str_df_list_MF<-list()
p_str_df_list_MF<-list()
matrix_boot_MF_fin_list<-list()
for(c in levels(matrix_boot_MF$country)){
  for(l in levels(matrix_boot_MF$str_fac)){
    yes<-matrix_boot_MF%>%
      dplyr::filter(country==c,str_fac==l,mut=="yes")
    no<-matrix_boot_MF%>%
      dplyr::filter(country==c,str_fac==l,mut=="no")
    if(nrow(yes)>=10 & nrow(no)>=10){
      Tobs<-mean(yes$Rt_14)-mean(no$Rt_14)
      B<-1000
      Tboot<-rep(0,B)
      for(i in 1:B){
        yes_X<-sample(c(yes$Rt,no$Rt_14),length(yes$Rt_14),replace=TRUE)
        no_Y<-sample(c(yes$Rt,no$Rt_14),length(no$Rt_14),replace=TRUE)
        Tboot[i]<-mean(yes_X)-mean(no_Y)
      }
      p<-sum(abs(Tboot)>=abs(Tobs))/B
      name<-paste(c,l,sep="_")
      str_df_list_MF[[name]]<-data.frame(str_fac=l,group1="no",group2="yes",p.value=p,p.signif=ifelse(p<0.05,"**","ns"),country=c)
      p_str_df_list_MF[[name]]<-ggplot(data.frame(x=Tboot),aes(x=x))+geom_histogram()+geom_vline(xintercept=Tobs,color="#039103")+geom_vline(xintercept=-Tobs,color="#039103")+ggtitle(paste("p =",p,c,l,sep=" "))
      matrix_boot_MF_fin_list[[name]]<-rbind(yes,no)
    }
  }
}
matrix_boot_MF_fin<-do.call(rbind,matrix_boot_MF_fin_list)
str_df_MF<-do.call(rbind,str_df_list_MF)
p_dist_MF<-ggarrange(plotlist=p_str_df_list_MF)
ggsave("Figure_S17.jpg",plot=p_dist_MF,dpi=500,height=25,width=25)
#plot of the Rt of MF
p_Rt_MF<-ggplot(data=subset(R_matrix_for_boot),aes(x=as.POSIXct(date),y=Rt,ymin=l_975,ymax=u_975,group=mut,color=mut))+facet_wrap(~country)+
  geom_bar(data=subset(R_matrix_for_boot,mut=="yes"&type%in%c("MF")),aes(y=(as.numeric(str_fac)-1)/2),stat="identity",color="gray84",fill="gray",alpha=0.5)+
  geom_line(data=subset(R_matrix_for_boot,mut=="yes"&type%in%c("MF")),color="#039103",size=1)+
  geom_line(data=subset(R_matrix_for_boot,mut=="no"&type%in%c("MF")),color="#5e7a5e",size=1)+
  geom_ribbon(data=subset(R_matrix_for_boot,mut=="yes"&type%in%c("MF")),fill="#039103",alpha=0.4,color=NA)+
  geom_ribbon(data=subset(R_matrix_for_boot,mut=="no"&type%in%c("MF")),fill="#5e7a5e",alpha=0.4,color=NA)+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b")+scale_y_continuous("Effective Reproduction Number (Rt)",limits=c(0,4.5),sec.axis=sec_axis(~.*2,name="Level of Stringency",breaks=seq(0,9,1)))+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8))+
  geom_hline(yintercept=1,linetype="dashed")

#plot of the significant differences
R_matrix_for_boot_MF<-R_matrix_for_boot%>%
  dplyr::filter(type=="MF")
R_matrix_for_boot_MF$mut<-gsub("yes","MF",R_matrix_for_boot_MF$mut)
R_matrix_for_boot_MF$mut<-gsub("no","notMF",R_matrix_for_boot_MF$mut)
leg_Rt_MF<-ggplot(R_matrix_for_boot_MF,aes(x=date,y=Rt,color=mut))+geom_line(size=3)+scale_color_manual(values=c("#039103","#5e7a5e"))+theme(legend.position="bottom")+labs(color="Mutation Type")
legend_Rt_MF<-get_legend(leg_Rt_MF)
p_Rt_MF_fin<-ggarrange(p_Rt_MF,common.legend=TRUE,legend.grob=legend_Rt_MF,legend="bottom")

p_effect_MF<-ggerrorplot(data=subset(matrix_boot,type=="MF"),size=0.2,width=0.5,error.plot="errorbar",x="str_fac",y="boot_mean",color="mut",desc_stat="median_q1q3",add="mean",add.params=list(group="mut",size=0.2,color="mut"),ggtheme=theme_bw(base_size=15))+
  facet_wrap(~country,scales="free_x")+scale_color_manual(values=c("#5e7a5e","#039103"))+labs(y="Effective Reproduction Number (Rt)",x="Level of Stringency")+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none")+
  scale_y_continuous(limits=c(min(subset(matrix_boot,type=="NY")$boot_mean),2.1))+
  stat_pvalue_manual(str_df_MF,x="str_fac",y.position=2.05,label="p.signif")

leg_effect_MF<-ggplot(R_matrix_for_boot_MF,aes(x=date,y=Rt,color=mut))+geom_point(size=3)+scale_color_manual(values=c("#039103","#5e7a5e"))+theme(legend.position="bottom")+labs(color="Mutation Type")
legend_effect_MF<-get_legend(leg_effect_MF)
p_effect_MF_fin<-ggarrange(p_effect_MF,common.legend=TRUE,legend.grob=legend_effect_MF,legend="bottom")

p_ef_MG_final<-ggarrange(p_Rt_MF_fin,p_effect_MF_fin,labels=letters,ncol=1)
ggsave(filename="Figure_S15.jpg",plot=p_ef_MG_final,dpi=500,width=10,height=15)
#HF(N501Y)
#Calculating significant difference of NY
matrix_boot_NY<-R_matrix_for_boot%>%
  dplyr::filter(type=="NY")
str_df_list_NY<-list()
p_str_df_list_NY<-list()
matrix_boot_NY_fin_list<-list()
for(c in levels(matrix_boot_NY$country)){
  for(l in levels(matrix_boot_NY$str_fac)){
    yes<-matrix_boot_NY%>%
      dplyr::filter(country==c,str_fac==l,mut=="yes")
    no<-matrix_boot_NY%>%
      dplyr::filter(country==c,str_fac==l,mut=="no")
    if(nrow(yes)>=10 & nrow(no)>=10){
      Tobs<-mean(yes$Rt_14)-mean(no$Rt_14)
      B<-1000
      Tboot<-rep(0,B)
      for(i in 1:B){
        yes_X<-sample(c(yes$Rt_14,no$Rt_14),length(yes$Rt_14),replace=TRUE)
        no_Y<-sample(c(yes$Rt_14,no$Rt_14),length(no$Rt_14),replace=TRUE)
        Tboot[i]<-mean(yes_X)-mean(no_Y)
      }
      p<-sum(abs(Tboot)>=abs(Tobs))/B
      name<-paste(c,l,sep="_")
      str_df_list_NY[[name]]<-data.frame(str_fac=l,group1="no",group2="yes",p.value=p,p.signif=ifelse(p<0.05,"**","ns"),country=c)
      p_str_df_list_NY[[name]]<-ggplot(data.frame(x=Tboot),aes(x=x))+geom_histogram()+geom_vline(xintercept=Tobs,color="#1857de")+geom_vline(xintercept=-Tobs,color="#1857de")+ggtitle(paste("p =",p,c,l,sep=" "))
      matrix_boot_NY_fin_list[[name]]<-rbind(yes,no)
    }
  }
}
matrix_boot_NY_fin<-do.call(rbind,matrix_boot_NY_fin_list)
str_df_NY<-do.call(rbind,str_df_list_NY)
p_dist_NY<-ggarrange(plotlist=p_str_df_list_NY)
ggsave("Figure_S18.jpg",plot=p_dist_NY,dpi=500,height=20.8,width=25)
#plot of the Rt of MF
p_Rt_NY<-ggplot(data=subset(R_matrix_for_boot),aes(x=as.POSIXct(date),y=Rt,ymin=l_975,ymax=u_975,group=mut,color=mut))+facet_wrap(~country)+
  geom_bar(data=subset(R_matrix_for_boot,mut=="no"&type%in%c("NY")),aes(y=(as.numeric(str_fac)-1)/2),stat="identity",color="gray84",fill="gray",alpha=0.5)+
  geom_line(data=subset(R_matrix_for_boot,mut=="yes"&type%in%c("NY")),color="#1857de",size=1)+
  geom_line(data=subset(R_matrix_for_boot,mut=="no"&type%in%c("NY")),color="#8697bd",size=1)+
  geom_ribbon(data=subset(R_matrix_for_boot,mut=="yes"&type%in%c("NY")),fill="#1857de",alpha=0.4,color=NA)+
  geom_ribbon(data=subset(R_matrix_for_boot,mut=="no"&type%in%c("NY")),fill="#8697bd",alpha=0.4,color=NA)+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b")+scale_y_continuous("Effective Reproduction Number (Rt)",limits=c(0,4.5),sec.axis=sec_axis(~.*2,name="Level of Stringency",breaks=seq(0,9,1)))+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8))+
  geom_hline(yintercept=1,linetype="dashed")
#plot of the significant differences
R_matrix_for_boot_NY<-R_matrix_for_boot%>%
  dplyr::filter(type=="NY")
R_matrix_for_boot_NY$mut<-gsub("yes","HF",R_matrix_for_boot_NY$mut)
R_matrix_for_boot_NY$mut<-gsub("no","notHF",R_matrix_for_boot_NY$mut)
leg_Rt_NY<-ggplot(R_matrix_for_boot_NY,aes(x=date,y=Rt,color=mut))+geom_line(size=3)+scale_color_manual(values=c("#1857de","#8697bd"))+theme(legend.position="bottom")+labs(color="Mutation Type")
legend_Rt_NY<-get_legend(leg_Rt_NY)
p_Rt_NY_fin<-ggarrange(p_Rt_NY,common.legend=TRUE,legend.grob=legend_Rt_NY,legend="bottom")

p_effect_NY<-ggerrorplot(data=subset(matrix_boot,type=="NY"),size=0.2,width=0.5,error.plot="errorbar",x="str_fac",y="boot_mean",color="mut",desc_stat="median_q1q3",add="mean",add.params=list(group="mut",size=0.2,color="mut"),ggtheme=theme_bw(base_size=15))+
  facet_wrap(~country,scales="free_x")+scale_color_manual(values=c("#8697bd","#1857de"))+labs(y="Effective Reproduction Number (Rt)",x="Level of Stringency")+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none")+
  scale_y_continuous(limits=c(min(subset(matrix_boot,type=="NY")$boot_mean),2.1))+
  stat_pvalue_manual(str_df_NY,x="str_fac",y.position=2.05,label="p.signif")

leg_effect_NY<-ggplot(R_matrix_for_boot_NY,aes(x=date,y=Rt,color=mut))+geom_point(size=3)+scale_color_manual(values=c("#1857de","#8697bd"))+theme(legend.position="bottom")+labs(color="Mutation Type")
legend_effect_NY<-get_legend(leg_effect_NY)
p_effect_NY_fin<-ggarrange(p_effect_NY,common.legend=TRUE,legend.grob=legend_effect_NY,legend="bottom")

#Calculations of change of Rt based on the change 14 days after##########################################
##
#Calculating change of Rt based on 14 days after
R_matrix_GR_delay_list<-list()
for(c in levels(R_matrix_final$country)){
  for(t in levels(R_matrix_final$type)){
    for(m in levels(R_matrix_final$mut)){
      GR_df<-R_matrix_final%>%
        dplyr::filter(country==c,type==t,mut==m)%>%
        dplyr::arrange(date)%>%
        dplyr::mutate(gr_R=0)
      for (i in 1:nrow(GR_df)){
        if (i %in% 1:14){
          GR_df$gr_R[i]<-NA
        }else if (i %in% (nrow(GR_df)-14):nrow(GR_df)){
          GR_df$gr_R[i]<-NA
          GR_df$gr_R[(i-14)]<-GR_df$Rt[i]-mean(GR_df$Rt[(i-14):(i-1)])
        }else{
          GR_df$gr_R[(i-14)]<-GR_df$Rt[i]-mean(GR_df$Rt[(i-14):(i-1)])
        }
        name<-paste(c,t,m,sep="_")
        R_matrix_GR_delay_list[[name]]<-GR_df
      }
    }
  }
}
R_matrix_GR_delay<-do.call(rbind,R_matrix_GR_delay_list)
R_matrix_GR_delay_for_boot<-R_matrix_GR_delay%>%
  dplyr::filter(!(y_m%in%c("2020-Feb","2020-Mar")))%>%
  tidyr::drop_na()
#calculating the bootstrap distribution of the mean of change of Rt by stringency level
matrix_GR_delay_boot_list<-list()
n_iter<-1000
for (c in levels(R_matrix_GR_delay_for_boot$country)){
  for (l in levels(R_matrix_GR_delay_for_boot$str_fac)){
    for (m in levels(R_matrix_GR_delay_for_boot$mut)){
      for (t in levels(R_matrix_GR_delay_for_boot$type)){
        dat_yes<-R_matrix_GR_delay_for_boot%>%
          dplyr::filter(country==c&str_fac==l&mut=="yes"&type==t)
        dat_no<-R_matrix_GR_delay_for_boot%>%
          dplyr::filter(country==c&str_fac==l&mut=="no"&type==t)
        dat<-R_matrix_GR_delay_for_boot%>%
          dplyr::filter(country==c&str_fac==l&mut==m&type==t)
        if(nrow(dat_yes)>=10 & nrow(dat_no)>=10){
          n<-nrow(dat)
          mean<-numeric(n_iter)
          for (iter in 1:n_iter){
            inc<-sample(1:n,replace=TRUE)
            mean[iter]<-mean(dat$gr_R[inc])
          }
          name<-paste(c,l,m,t,sep="_")
          matrix_GR_delay_boot_list[[name]]<-data.frame(country=c,str_fac=l,type=t,mut=m,gr_boot_mean=mean)
        }
      }
    }
  }
}
matrix_GR_delay_boot<-do.call(rbind,matrix_GR_delay_boot_list)
matrix_GR_delay_boot$country<-as.factor(matrix_GR_delay_boot$country)
matrix_GR_delay_boot$str_fac<-as.factor(matrix_GR_delay_boot$str_fac)
matrix_GR_delay_boot$type<-as.factor(matrix_GR_delay_boot$type)
matrix_GR_delay_boot$mut<-as.factor(matrix_GR_delay_boot$mut)

#Plot of change of Rt##########################################
##
#Differences between change of Rt in different stringency levels
matrix_boot_GR_delay_NY<-R_matrix_GR_delay_for_boot%>%
  dplyr::filter(type=="NY")
str_df_list_GR_delay_NY<-list()
p_str_df_list_GR_delay_NY<-list()
matrix_boot_GR_delay_NY_fin_list<-list()
for(c in levels(matrix_boot_GR_delay_NY$country)){
  for(l in levels(matrix_boot_GR_delay_NY$str_fac)){
    yes<-matrix_boot_GR_delay_NY%>%
      dplyr::filter(country==c,str_fac==l,mut=="yes")
    no<-matrix_boot_GR_delay_NY%>%
      dplyr::filter(country==c,str_fac==l,mut=="no")
    if(nrow(yes)>=10 & nrow(no)>=10){
      Tobs<-mean(yes$gr_R)-mean(no$gr_R)
      B<-1000
      Tboot<-rep(0,B)
      for(i in 1:B){
        yes_X<-sample(c(yes$gr_R,no$gr_R),length(yes$gr_R),replace=TRUE)
        no_Y<-sample(c(yes$gr_R,no$gr_R),length(no$gr_R),replace=TRUE)
        Tboot[i]<-mean(yes_X)-mean(no_Y)
      }
      p<-sum(abs(Tboot)>=abs(Tobs))/B
      name<-paste(c,l,sep="_")
      str_df_list_GR_delay_NY[[name]]<-data.frame(str_fac=l,group1="no",group2="yes",p.value=p,p.signif=ifelse(p<0.05,"**","ns"),country=c)
      p_str_df_list_GR_delay_NY[[name]]<-ggplot(data.frame(x=Tboot),aes(x=x))+geom_histogram()+geom_vline(xintercept=Tobs,color="#1857de")+geom_vline(xintercept=-Tobs,color="#1857de")+ggtitle(paste("p =",p,c,l,sep=" "))
      matrix_boot_GR_delay_NY_fin_list[[name]]<-rbind(yes,no)
    }
  }
}
matrix_boot_GR_delay_NY_fin<-do.call(rbind,matrix_boot_GR_delay_NY_fin_list)
str_df_GR_delay_NY<-do.call(rbind,str_df_list_GR_delay_NY)
p_dist_GR_delay_NY<-ggarrange(plotlist=p_str_df_list_GR_delay_NY)
ggsave("Figure_S19.jpg",plot=p_dist_GR_delay_NY,dpi=500,height=20.8,width=25)
#plot of change of Rt in time
p_gr_Rt_delay_NY<-ggplot(data=subset(R_matrix_GR_delay_for_boot),aes(x=as.POSIXct(date),y=gr_R,group=mut,color=mut))+facet_wrap(~country)+
  geom_segment(data=subset(R_matrix_GR_delay_for_boot,mut=="no"&type%in%c("NY")),aes(y=((as.numeric(str_fac)-1)/2.25)-2,xend=as.POSIXct(date),yend=-2),color="gray84",alpha=0.5)+
  geom_line(data=subset(R_matrix_GR_delay_for_boot,mut=="yes"&type%in%c("NY")),color="#1857de",size=1)+
  geom_line(data=subset(R_matrix_GR_delay_for_boot,mut=="no"&type%in%c("NY")),color="#8697bd",size=1)+
  scale_x_datetime(date_breaks="1 month",date_labels="%Y-%b")+scale_y_continuous("Change of Rt",limits=c(-2,2),sec.axis=sec_axis(trans=~(.+2)*2.25,name="Level of Stringency",breaks=seq(0,9,1)))+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none",axis.title.x=element_blank(),axis.text.x=element_text(angle=45,hjust=1,vjust=1,size=8))+
  geom_hline(yintercept=0,linetype="dashed")
p_gr_Rt_delay_NY_fin<-ggarrange(p_gr_Rt_delay_NY,common.legend=TRUE,legend.grob=legend_Rt_NY,legend="bottom")
#plot of differences between change of Rt in different stringency levels
p_gr_delay_effect_NY<-ggerrorplot(data=subset(matrix_GR_delay_boot,type=="NY"),size=0.2,width=0.5,error.plot="errorbar",x="str_fac",y="gr_boot_mean",color="mut",desc_stat="median_q1q3",add="mean",add.params=list(group="mut",size=0.2,color="mut"),ggtheme=theme_bw(base_size=15))+
  facet_wrap(~country,scales="free_x")+scale_color_manual(values=c("#8697bd","#1857de"))+labs(y="Change of Rt",x="Level of Stringency")+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none")+
  scale_y_continuous(limits=c(-0.45,0.45))+
  stat_pvalue_manual(str_df_GR_delay_NY,x="str_fac",y.position=0.4,label="p.signif")
p_gr_delay_effect_NY_fin<-ggarrange(p_gr_delay_effect_NY,common.legend=TRUE,legend.grob=legend_effect_NY,legend="bottom")
ggsave("Figure_S20.jpg",plot=p_gr_delay_effect_NY_fin,dpi=500,width=10,height=7.5)
#Plot of difference between Rt in different stringency levels and plot of change of Rt in time
p_gr_delay_ef_NY_final<-ggarrange(p_effect_NY_fin,p_gr_Rt_delay_NY_fin,labels=letters,ncol=1)
ggsave("Figure_4.jpg",plot=p_gr_delay_ef_NY_final,dpi=500,width=10,height=15)

#Correlation between stringency levels and Rt 14 days after##########################################
##
#HF(N501Y)
matrix_boot_GR_delay_NY_fin_14days_list<-list()
for (c in levels(matrix_boot_GR_delay_NY_fin$country)){
  for (m in levels(matrix_boot_GR_delay_NY_fin$mut)){
    df<-matrix_boot_GR_delay_NY_fin%>%
      dplyr::filter(country==c&mut==m)%>%
      dplyr::arrange(date)%>%
      dplyr::mutate(Rt_14=0)
    for (i in 1:nrow(df)){
      if (i %in% (nrow(df)-13):nrow(df)){
        df$Rt_14[i]<-NA
      }else{
        df$Rt_14[i]<-df$Rt[i+13]
      }
    name<-paste(c,m,sep="_")
    matrix_boot_GR_delay_NY_fin_14days_list[[name]]<-df
    }
  }
}
matrix_boot_GR_delay_NY_fin_14days<-do.call(rbind,matrix_boot_GR_delay_NY_fin_14days_list)
p_glm_rt_str_NY_14<-ggplot(data=subset(matrix_boot_GR_delay_NY_fin_14days,type=="NY"&!is.na(Rt_14)),aes(x=cov_str/10,y=Rt_14,color=mut))+
  geom_point(alpha=0.5)+facet_wrap(~country)+stat_smooth(method="glm",aes(fill=mut))+
  theme(strip.background=element_rect(color=NA,fill="white"),strip.text=element_text(hjust=0),legend.position="none")+
  scale_y_continuous(limits=c(0,4))+scale_x_continuous(limits=c(3,9),breaks=seq(0,9,1))+
  scale_color_manual(values=c("#8697bd","#1857de"))+scale_fill_manual(values=c("#8697bd","#1857de"))+
  labs(y="Effective Reproduction Number (Rt)",x="Level of Stringency")+
  geom_hline(yintercept=1,linetype="dashed")+
  stat_cor(aes(label=paste(..r.label..,..rr.label..,..p.label..,sep="~`,`~")),method="spearman",p.accuracy=0.001,r.accuracy=0.01)
p_glm_rt_str_NY_fin_14<-ggarrange(p_glm_rt_str_NY_14,common.legend=TRUE,legend.grob=legend_effect_NY,legend="bottom")

p_ef_NY_final<-ggarrange(p_Rt_NY_fin,p_glm_rt_str_NY_fin_14,labels=letters,ncol=1)
ggsave("Figure_3.jpg",plot=p_ef_NY_final,dpi=500,width=10,height=15)

