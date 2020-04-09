#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# code to analyze dPCR data iGlut/Gaba experiments
# 
###############################################

#clear environment
rm(list=ls())

#load libraries
library(readxl)
library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(stringr)

#functions
`%notin%` <- Negate(`%in%`)

#read in dPCR data and meta data table
full.table <- read.table(file = "path.to.file/ddPCR/20200129.full.table.iGaN.iGABA.tab", stringsAsFactors = F, sep = "\t")
# force these variables to be numeric
full.table$Concentration_in_CY5 <- as.numeric(full.table$Concentration_in_CY5)
full.table$Concentration_in_VIC <- as.numeric(full.table$Concentration_in_VIC)

gd <- full.table %>% 
  group_by(EH_sample_ID, MH_Sample_ID, Sample_Type, HA_cell_type, V5_cell_type, IP_type, HA_donor, V5_donor) %>% 
  summarise(GFP = mean(Concentration_in_VIC),
            mCherry = mean(Concentration_in_CY5), 
            GFP.sd = sd(Concentration_in_VIC), 
            mCherry.sd = sd(Concentration_in_CY5))

gd$Culture_Type <- c("coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "coculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture", "monoculture")

# calculate ratios
gd$mCherry_TO_GFP <- gd$mCherry/gd$GFP 
gd$GFP_TO_mCherry <- gd$GFP/gd$mCherry 

# factor sample numbers
EH.samples <- as.character(seq(from = 1, to = 34, by = 1))
gd$EH_sample_ID <- factor(gd$EH_sample_ID, levels = EH.samples)

#define sample types and factor
sample.types <- c("IP HA", "IC", "IP V5")
gd$Sample_Type <- factor(gd$Sample_Type, levels = sample.types)

# make a combined EH and MH label so replicates are more easily seen
gd$combined.label <- paste(gd$EH_sample_ID, gd$MH_Sample_ID, sep = "_")

# make a combined donor column
gd$combined.donor <- paste(gd$HA_donor, gd$V5_donor, sep = "")

cocultures <- gd %>% filter(Culture_Type == "coculture")
cocultures$IP_celltype <- c("none", "none", "none", "none", "none", "none", "iGluN", "iGluN", "iGaN", "iGaN", "iGluN", "iGaN", "iGluN", "iGaN", "iGluN", "iGaN", "iGluN", "iGaN")
ip.celltypes <- c("iGluN", "none", "iGaN")
cocultures$IP_celltype <- factor(cocultures$IP_celltype, levels = ip.celltypes)
sample.types <- c("IP HA", "IC", "IP V5")
cocultures$Sample_Type <- factor(cocultures$Sample_Type, levels = sample.types)

coculture.subset <- cocultures[-c(1:2,7:10),] # gets rid of cocultures composed of the identical cell types with opposite tags (e.g. iGaN RPL22-HA with iGaN RPL22-V5)
ggplot(coculture.subset, aes(x=IP_celltype, y=mCherry_TO_GFP, color=Sample_Type)) +
  geom_point(position = "identity", size = 6) +
  scale_y_continuous(trans='log2', 
                     breaks = c(0.015625,0.0625,0.25,1,4,16,64),
                     limits = c(0.015625,64),
                     labels = trans_format("log2", math_format(2^.x))) +
  #scale_y_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
  geom_line(aes(group = MH_Sample_ID), color="black") +
  scale_color_manual(name="Antibody",
                     values=c("IP HA" = "orchid3", 
                              "IC" = "gold", 
                              "IP V5" = "deepskyblue"),
                     labels = c("IP HA" = "HA",
                                "IC" = "none",
                                "IP V5" = "V5")) +
  ylab(expression(paste("[", italic("mCherry"), "]/[", italic("GFP"), "]", sep = ""))) +
  theme_minimal() +
  scale_x_discrete(labels = c("iGluN IP", "co-culture", "iGaN IP"))

