# here we will merge data into one unified data set
# P. Alexander Burnham
# 8 April 2026

# Read in data
spatial_dist <- read.csv("data/beeDist.csv", header = TRUE)
weight <- read.csv("data/beeWeight.csv", header = TRUE)
pathogens <- read.csv("data/overwinteringPathogen_data.csv", header = TRUE)


# Read in libraries
library(ggplot2)
library(lme4)
library(dplyr)

# Merge weight and pathogens
df <- merge(x = pathogens, y = weight, by = "ID", all.x = T)

# remove group
df$Group <- NULL

# rename treatment to season
names(df)[2] <- "Season"

# Export Data
#write.csv(df, "data/overwintering_data_clean.csv")
