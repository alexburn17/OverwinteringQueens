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
setwd("~/Desktop/GitHub_air_05042022/PA_Burnham/OverwinteringQueens")

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
beeDist <- read.table("beeDist.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 
beeWeight <- read.table("beeWeight.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 
multPlot <- read.table("overwinteringMultiPlot.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

multPlot$Binary <- ifelse(multPlot$load>0, 1, 0)
multPlot$logLoad <- log10(multPlot$load+1)

multplotNo0 <- multPlot[!multPlot$logLoad==0,]
multplotNo0 <- multplotNo0[!is.na(multplotNo0$logLoad),]


dat$Binary <- ifelse(dat$NormGenomeCopy>0, 1, 0)
dat$logVirus <- log10(dat$NormGenomeCopy+1)

datNo0 <- dat[!dat$logVirus==0,]

# split data set by virus
dat_split <- split(dat, dat$Target.Name)
datno0_split <- split(datNo0, datNo0$Target.Name)


# BQCV Stats

chisq.test(dat_split$BQCV$Treatment, dat_split$BQCV$Binary)
t.test(log10(datno0_split$BQCV$NormGenomeCopy+1)~datno0_split$BQCV$Treatment)


# Nosema Stats
multplotSplit <- split(multPlot, multPlot$pathogen)


chisq.test(multplotSplit$DWV$Binary, multplotSplit$DWV$Treatment)
t.test(log10(datno0_split$BQCV$Nosema +1)~datno0_split$BQCV$Treatment)


# depth under the soil distribution 
ggplot(data=beeDist, aes(x=depth_cm)) +
  geom_histogram( binwidth=2, fill="steelblue", color="#e9ecef", alpha=0.9) +
  labs(y="Frequency", x="Depth Under the Soil (cm)") +
  theme_minimal(base_size = 17) 



ggplot(beeWeight, aes(x = Group, y = Bee_Mass_Est_g., fill=Group)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha=1)+
  labs(y="Bee Mass (g)", x="Season") +
  theme_minimal(base_size = 17) +
  scale_fill_manual(values = c("goldenrod", "chartreuse4"))






dat <- aov(data=beeWeight, Bee_Mass_Est_g.~Group)
summary(dat)
kruskal.test(beeWeight$Bee_Mass_Est_g.~beeWeight$Group)



library(MASS)
library(scales)

ggplot(multplotNo0, aes(x = pathogen, y = load, fill=Treatment)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(jitter.width = 0.1), alpha=1)+
  labs(y="Pathogen Load log10(per/bee)", x="Pathogen") +
  theme_minimal(base_size = 17) +
  scale_fill_manual(values = c("goldenrod", "chartreuse4"))+
  theme(legend.position = 'top') +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) 


#Checking out by plant species
data <- ddply(multPlot, c("Treatment", "pathogen"), summarise, 
              n = length(Binary),
              a = sum(Binary, na.rm = T)+1,
              b = n - a + 1,
              lower = qbeta(.025, shape1 = a, shape2 = b),
              upper = qbeta(.975, shape1 = a, shape2 = b),
              mean = mean(Binary, na.rm=TRUE))




ggplot(data, aes(x=pathogen, y=mean, group=Treatment, fill=Treatment)) +
  geom_point(stat="identity", position=position_dodge(width = .9), shape=23, size=13, color = "black") + 
  labs(y="Pathogen Prevalence", x="Pathogen") + 
  theme_minimal(base_size = 17) + 
  theme(legend.position="top") +
  scale_fill_manual(values=c("goldenrod", "chartreuse4")) +
  geom_errorbar(aes(ymin = lower, ymax = upper, width = 0.5), 
                position=position_dodge(width = .9)) +
  coord_cartesian(ylim=c(0,1.3))

  

data






#Nosema load by treatment plot
boxplot(OWQ_DWV$Nosema~OWQ_DWV$Treatment,
        ylab="Nosema Load", 
        xlab="Season")
#Nosema load by treatment plot
summary(aov(OWQ_DWV$Nosema~OWQ_DWV$Treatment))
#Nosema load by treatment plot
kruskal.test(OWQ_DWV$Nosema~as.factor(OWQ_DWV$Treatment))

boxplot(log10(OWQ_DWV$NormGenomeCopy)~OWQ_DWV$Treatment)
summary(aov(OWQ_DWV$NormGenomeCopy~OWQ_DWV$Treatment))
kruskal.test(OWQ_DWV$NormGenomeCopy~as.factor(OWQ_DWV$Treatment))


boxplot(log10(OWQ_BQCV$NormGenomeCopy)~OWQ_BQCV$Treatment, 
        ylab="BQCV Load log10(genome copies/bee)", 
        xlab="Season")


summary(aov(OWQ_BQCV$NormGenomeCopy~OWQ_BQCV$Treatment))
kruskal.test(OWQ_BQCV$NormGenomeCopy~as.factor(OWQ_BQCV$Treatment))


fisher.test(OWQ_BQCV$Binary, OWQ_BQCV$Treatment)




