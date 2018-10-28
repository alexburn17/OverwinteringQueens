###########################################################################################
# Data Analysis for 2016 Bombus Overwintering Queen Study
# Kim Mack-Nair 
# April 12, 2018
###########################################################################################

#Preliminaries:
# Clear memory of characters
ls()
rm(list=ls())
#install.packages("ggplot2")
#install.packages("dplyr")
#install.packages("plyr")
#install.packages("lme4")
#install.packages("car")


# Call Required Packages
library("ggplot2")
library("dplyr")
library("plyr")
library("lme4")
library("car")

# Set Working Directory: 
# setwd("~/Desktop/OverwinteringQueens")

# read in data:
Bee_Data <- read.table("Bee_Data.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

OverwinteringQueens_Results <- read.table("OverwinteringQueens_Results.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE)

# merge data sets
OverwinteringQueens_Results <- merge(Bee_Data, OverwinteringQueens_Results, by=c("Sample.Name"), all.y=TRUE )

###########################################################################################
# Functions

###########################################################################
# function name: PrelimClean (requires dplyr package)
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned 
###########################################################################

PrelimClean <- function(data=data){
  
  # take only columns that we want:
  data <- dplyr::select(data, Sample.Name, Target.Name, Cq.Mean, Quantity.Mean, Treatment, Nosema, dil.factor)
  
  # remove duplicate rows
  data<-data[!duplicated(data),]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$Sample.Name=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$Sample.Name=="Gblock"),]
  # remove Gblock rows from dataframe:
  data<-data[!(data$Sample.Name=="gBlocks"),]
  
  return(data)
}

###########################################################################
# END OF FUNCITON
###########################################################################



###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr<- 100
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  total_extr_vol<- 600
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  

  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$Quantity.Mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################






###########################################################################
# function name: actinNormal
# description: normalizes virus data with actin values 
# parameters: a data frame with actin values
# returns a dataframe with normalized virus values 
###########################################################################

actinNormal <- function(data=data){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$Target.Name=="ACTIN"),]
  
  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$Sample.Name, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("Sample.Name", "ACT_genomeCopy")
  
  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("Sample.Name"), all.x=TRUE)
  
  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)
  
  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean
  
  return(data)
}


###########################################################################
# END OF FUNCITON
###########################################################################




###########################################################################
# function name: CT_Threash
# description: creates binary data and makes genome copy 0 if below Ct threash
# parameters: dataframe
###########################################################################

CT_Threash <- function(data=data){
  
  splitDF <- split(data, data$Target.Name)
  
  # make norm_genome_copbee 0 if Ct value is > 32.918
  splitDF$DWV$NormGenomeCopy[which(splitDF$DWV$Cq.Mean > 32.918)] <- 0
  splitDF$BQCV$NormGenomeCopy[which(splitDF$BQCV$Cq.Mean > 32.525)] <- 0
  
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Cq.Mean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Cq.Mean > 32.525, 0, 1)
  
  # merge split dataframe back into "BombSurv" dataframe:
  data <- rbind(splitDF$DWV, splitDF$BQCV)
  
  # norm_genomeCopy is 0 if NA
  data$virusBINY[is.na(data$virusBINY)] <- 0
  
  return(data)
  
}



###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# running functions to clean data set:
OverwinteringQueens_Results <- PrelimClean(OverwinteringQueens_Results)
OverwinteringQueens_Results<- VirusNorm(number_bees = 1, OverwinteringQueens_Results)
OverwinteringQueens_Results<- actinNormal(OverwinteringQueens_Results)
OverwinteringQueens_Results<- CT_Threash(OverwinteringQueens_Results)



#write.csv(OverwinteringQueens_Results, "overwintering.csv")











##############################################################################
# PRELIM ANALYSIS:


# read in data:
dat <- read.table("overwintering.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

dat$Binary <- ifelse(dat$NormGenomeCopy>0, 1, 0)


# split data set by virus
OWQ_split <- split(dat, dat$Target.Name)



# create BQCV data set
OWQ_BQCV <- OWQ_split$BQCV

# create DWV data set
OWQ_DWV <- OWQ_split$DWV


#Nosema load by treatment plot
boxplot(OWQ_DWV$Nosema~OWQ_DWV$Treatment)
#Nosema load by treatment plot
summary(aov(OWQ_DWV$Nosema~OWQ_DWV$Treatment))
#Nosema load by treatment plot
kruskal.test(OWQ_DWV$Nosema~as.factor(OWQ_DWV$Treatment))

boxplot(log10(OWQ_DWV$NormGenomeCopy)~OWQ_DWV$Treatment)
summary(aov(OWQ_DWV$NormGenomeCopy~OWQ_DWV$Treatment))
kruskal.test(OWQ_DWV$NormGenomeCopy~as.factor(OWQ_DWV$Treatment))

boxplot(log(OWQ_BQCV$NormGenomeCopy)~OWQ_BQCV$Treatment)
summary(aov(OWQ_BQCV$NormGenomeCopy~OWQ_BQCV$Treatment))
kruskal.test(OWQ_BQCV$NormGenomeCopy~as.factor(OWQ_BQCV$Treatment))


fisher.test(OWQ_BQCV$Binary, OWQ_BQCV$Treatment)




