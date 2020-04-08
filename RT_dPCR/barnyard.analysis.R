#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# code to analyze barnyard digital PCR data
###############################################

#clear environment
rm(list=ls())


library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)

#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# code to analyze barnyard ddPCR data github version
###############################################

#clear environment
rm(list=ls())


library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)

# load dPCR and meta data info
A <- read.table(file = "path.to.file/20190221.barnyard.fulltable.tab", stringsAsFactors = F, sep = "\t")
colnames(A) <- A[1,] # create column names
normalize.A <- A[2:NROW(A),] # copy dataframe for making calculations
# change to numeric. NAs cause it to be loaded as a character 
normalize.A$Concentration.in.CY5 <- as.numeric(normalize.A$Concentration.in.CY5)
normalize.A$Concentration.in.VIC <- as.numeric(normalize.A$Concentration.in.VIC)
#subtract background
normalize.A$Concentration.in.VIC <- normalize.A$Concentration.in.VIC - mean(normalize.A$Concentration.in.VIC[16:17]) #subtract mean of negative control
normalize.A$Concentration.in.CY5 <- normalize.A$Concentration.in.CY5 - mean(normalize.A$Concentration.in.CY5[13:15]) #subtract mean of negative control
normalize.A <- normalize.A[-c(13,14,15,16,17),] # remove control rows

# change to long format
test <- melt(normalize.A, id.vars = c("IPantibody", "SampleType"), measure.vars = c("Concentration.in.CY5", "Concentration.in.VIC"))
test$SampleType <- factor(test$SampleType,levels = c("Input Control", "HA IP 1", "V5 IP 2", "V5 IP 1", "HA IP 2")) # put in the order I want them plotted
test$variable <- gsub("Concentration.in.CY5", "mCherry", test$variable) # change to gene name
test$variable <- gsub("Concentration.in.VIC", "GFP", test$variable)

# summarize data for plotting
gd <- normalize.A %>% 
  group_by(SampleType) %>% 
  summarise(Concentration.in.VIC = mean(Concentration.in.VIC),
            Concentration.in.CY5 = mean(Concentration.in.CY5))
colnames(gd) <- c("Sample", "GFP", "mCherry")
# calculate ratios
gd$mCherry_TO_GFP <- gd$mCherry/gd$GFP 
gd$GFP_TO_mCherry <- gd$GFP/gd$mCherry
gd$Sample <- factor(gd$Sample,levels = c("Input Control", "HA IP 1", "V5 IP 2", "V5 IP 1", "HA IP 2")) # put in the order I want them plotted

ggplot(gd, aes(x=Sample, y=mCherry_TO_GFP, fill=Sample)) +
  ylab(expression(paste("[", italic("mCherry"), "]/[", italic("GFP"), "]", sep = ""))) +
  scale_y_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) + 
  geom_col(position = "dodge") +
  scale_fill_manual(name="Sample",
                    values=c("Input Control" = "gold", 
                             "HA IP 1" = "deepskyblue",
                             "V5 IP 2" = "orchid3",
                             "V5 IP 1" = "orchid3",
                             "HA IP 2" = "deepskyblue")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
