###############################################
#
# James Gregory (jgregory@nygenome.org), ehoelzli@nygenome.org
# code to analyze digital PCR from organoid/microglia co-cultures
# 
###############################################

#clear environment
rm(list=ls())

#functions
`%notin%` <- Negate(`%in%`) # another way for a 'not' in function

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)

dPCR.data <- read.table(file = "/path.to.file/20191203organoid.microglia.tab",
                        sep = "\t", 
                        stringsAsFactors = F)

           
# put means and standard deviations into a table
dPCR.summary <- dPCR.data %>% 
  group_by(matched_samples, sample_names, IPantibody, cell_type, culture_type) %>% 
  summarise(GFP = mean(Concentration.in.VIC),
            MTCO2 = mean(Concentration.in.CY5), 
            GFP.sd = sd(Concentration.in.VIC), 
            MTCO2.sd = sd(Concentration.in.CY5))

dPCR.summary$ratio <- dPCR.summary$GFP/dPCR.summary$MTCO2 # calculate ratio
dPCR.summary$IPantibody <- gsub("none", "input", dPCR.summary$IPantibody)
dPCR.summary$IPantibody <- factor(dPCR.summary$IPantibody, levels=c("input", "HA"))

# select organoid cell types
dPCR.summary.microglia.cocultures <- dPCR.summary %>% filter(cell_type %notin% "organoid")
ggplot(dPCR.summary.microglia.cocultures, aes(x=IPantibody, y=ratio, color=cell_type)) +
  geom_point(size=3) +
  geom_line(aes(group = matched_samples), color="black") +
  scale_y_continuous(trans='log2', 
                     breaks = c(0.015625,0.0625,0.25,1,4,16,64),
                     labels = trans_format("log2", math_format(2^.x))) +
  scale_color_manual(name="Culture Type",
                     values=c("microglia" = "#FF67A4", 
                              "microglia_organoid" = "#619CFF"),
                     labels = c("microglia" = "microglia monoculture", 
                                "microglia_organoid" = "microglia/organoid co-culture")) +
  ylab(expression(paste("[", italic("GFP"), "]/[", italic("MTCO2"), "]", sep = ""))) +
  theme_minimal() 

# melt data into long format for plotting with ggplot
test <- melt(dPCR.summary.microglia.cocultures, id.vars = c("IPantibody", "matched_samples", "cell_type"),
             measure.vars = c("GFP", "MTCO2"))
ggplot(test, aes(x=IPantibody, y=value, color=cell_type, shape=variable)) +
  geom_jitter(size=4, alpha=0.6, width=0.025) +
  scale_color_manual(name="Culture Type",
                     values=c("microglia" = "#FF67A4", 
                              "microglia_organoid" = "#619CFF"),
                     labels = c("microglia" = "microglia monoculture", 
                                "microglia_organoid" = "microglia/organoid co-culture")) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x)),
                breaks = c(1,10,100,1000,10000),
                limits = c(10,10000)) +
  scale_shape_manual(name="gene",
                     labels = c("GFP" = "GFP",
                                "MTCO2"="MTCO2"),
                     values = c("GFP"=16,
                                "MTCO2"=17)) +
  geom_line(data=subset(test, variable == "GFP"), aes(group=matched_samples), color="black") +
  geom_line(data=subset(test, variable == "MTCO2"), aes(group=matched_samples), color="black") +
  theme_minimal() +
  ylab(expression(paste("[", italic("RNA"), "] copies/ul", sep = ""))) 
