#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# calculate depletion efficiency for MN/astro cocultures
###############################################

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(stringr)
library(limma)
library(lme4)


#clear environment
rm(list=ls())
project <- "path.to.working.dir" 

#functions
'%!in%' <- function(x,y)!('%in%'(x,y))  # create a function that is the opposite of %in%; see this URL https://www.r-bloggers.com/the-notin-operator/
`%notin%` <- Negate(`%in%`) # another way for a 'not' in function

#load meta data
meta.data <- read.table(file = "path.to.file/20190918/20200422_MNastro_meta.data.tab", 
                        sep = "\t", 
                        stringsAsFactors = F,
                        header = T)

#load gene count matrix (per gene)
count.matrix <- read.table(file = file.path(project, "20200422.count.matrix.hg38.mm20.tab"), 
                           header = T, 
                           sep = "\t", 
                           stringsAsFactors = F)

# load samtools count matrix (per chromosome)
samtools.count.matrix <- read.table(file = file.path(project, "20200423.samtools.counts.MNastro.tab"), 
                                    sep = "\t",
                                    stringsAsFactors = F,
                                    header = T)

# vector of matched libraries by row (co.culture, MN IP, mAstro IP)
matched.cocultures <- c("Sample_lib2", "Sample_lib5", "Sample_lib8",
                        "Sample_lib10", "Sample_lib12", "Sample_lib14",
                        "Sample_lib17", "Sample_lib20", "Sample_lib18", 
                        "Sample_lib21", "Sample_lib28",
                        "Sample_lib23", "Sample_lib30", 
                        "Sample_lib25", "Sample_lib31")

# select only libraries from matched cocultures
meta.data.subset <- meta.data %>% filter(library_name %in% matched.cocultures)
chromosomes <- samtools.count.matrix$chromosome
samtools.count.matrix.subset <- samtools.count.matrix %>% select(all_of(matched.cocultures))
count.matrix <- count.matrix %>% select(all_of(c(matched.cocultures, "ensembl_gene_id", "species")))


###############################
## calculate depletion efficiency for all possible samples
###############################

# generate a list of all the matched samples; these are indicated by color in the paired.samples column
sample.lib <- "Sample_lib"
coculture1 <- c(2,8,5)
coculture2 <- c(10,12,14)
coculture3 <- c(17,18,20)
coculture4 <- c(21,28)
coculture5 <- c(23,30)
coculture6 <- c(25,31)
coculture.list <- list(coculture1, coculture2, coculture3, coculture4, coculture5, coculture6)

# function for putting together library names
paste.names <- function(x) {  
  y <- paste(sample.lib, x, sep = "") 
  return(y)
}
coculture.matched.list <- lapply(coculture.list, paste.names)

# this function calculates the change in uniquely mapping reads (IP - input)
normalize.to.input <- function(x) {
  # select desired columns, each list has the matched samples
  subset.dataframe <- samtools.count.matrix %>% select(all_of(x))
  column.number <- NCOL(subset.dataframe) # get column number; need this later
  for (i in 1:column.number) {
    # first calculate the percentage of the library for each chromosome
    subset.dataframe[,paste(colnames(subset.dataframe[i]), ".normalized", sep = "")] <- subset.dataframe[i]/sum(subset.dataframe[i])
  }
  #subest just the normalized data
  normalized.coculture <- subset.dataframe[,c((column.number+1):NCOL(subset.dataframe))]
  cols.for.normalized <- NCOL(normalized.coculture) # get column number; need this later
  # calculate difference IP - Input
  for (n in 2:cols.for.normalized) {
    normalized.coculture[,paste(colnames(normalized.coculture[n]), ".difference", sep = "")] <- (normalized.coculture[n] - normalized.coculture[1])
  }
  return(normalized.coculture)
}

#samtools.count.matrix.subset$chromosome <- samtools.count.matrix$chromosome # add back the chromosome
MN.astro.normalized.difference <- lapply(coculture.matched.list, normalize.to.input)

# the resulting data from do.call cbind is a dataframe with the number of reads for each chromosome divided
# by the total reads in the library (these columns are called normalized). The columns labelled difference are the IPs
# minus the input co-culture
combined.dataframe <- do.call(cbind, MN.astro.normalized.difference) # change list of lists to a dataframe

combined.dataframe <- combined.dataframe*100 # change to percent from fraction
combined.dataframe$chromosome <- samtools.count.matrix$chromosome # add chromosome labels back
# add column for species
combined.dataframe$species <- "mouse" # set mouse for all of them
combined.dataframe$species[23:nrow(samtools.count.matrix)] <- "human" # replace appropriate rows with human
#combined.dataframe$species <- inputcontrol.data.frame$cell.species # add species info
combined.dataframe$species[48:50] <- "reporter"
combined.dataframe <- combined.dataframe[c(1:47,49,50),] # remove Cre, which was not used in this work 
combined.dataframe$mito <- "no" # generically label all of them
combined.dataframe$mito[c(22,47)] <- "MT" #specifically label the mito chromosome
canonical.chromosomes <- combined.dataframe %>% 
  filter(chromosome %notin% c("mCherry", "GFP")) %>%
  filter(mito %notin% "MT")


######################
## plot 
#########################
# HA IP for 20190603 MN differentiation co-culture - mAstros
ggplot(combined.dataframe, aes(x=Sample_lib2.normalized, 
                               y=Sample_lib5.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

# calculate linear regression
HA.20190603.mAstro <- lmList(Sample_lib5.normalized.difference ~ Sample_lib2.normalized|species, data=canonical.chromosomes)
summary(HA.20190603.mAstro[[1]])


# V5 IP for 20190603 MN differentiation co-culture - MNs (<1ng yield on IP)
ggplot(combined.dataframe, aes(x=Sample_lib2.normalized, 
                               y=Sample_lib8.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

V5.20190603.MN <- lmList(Sample_lib8.normalized.difference ~ Sample_lib2.normalized|species, data=canonical.chromosomes)
summary(V5.20190603.MN)
summary(V5.20190603.MN[[2]])

# HA IP for 20190528 MN differentiation - mAstros
ggplot(combined.dataframe, aes(x=Sample_lib10.normalized, 
                               y=Sample_lib12.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

HA.20190528.mAstro <- lmList(Sample_lib12.normalized.difference ~ Sample_lib10.normalized|species, data=canonical.chromosomes)
summary(HA.20190528.mAstro[[1]])

# V5 IP for 20190528 MN differentiation - MNs
ggplot(combined.dataframe, aes(x=Sample_lib10.normalized, 
                               y=Sample_lib14.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

V5.20190528.MN <- lmList(Sample_lib14.normalized.difference ~ Sample_lib10.normalized|species, data=canonical.chromosomes)
summary(V5.20190528.MN)
summary(V5.20190528.MN[[2]])

#HA IP from MN differentiation 20190417 (mAstros)
ggplot(combined.dataframe, aes(x=Sample_lib17.normalized, 
                               y=Sample_lib20.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

HA.20190417.mAstro <- lmList(Sample_lib20.normalized.difference ~ Sample_lib17.normalized|species, data=canonical.chromosomes)
summary(HA.20190417.mAstro)
summary(HA.20190417.mAstro[[1]])

#V5 IP from MN differentiation 20190417 (MNs) (<1ng yield on IP)
ggplot(combined.dataframe, aes(x=Sample_lib17.normalized, 
                               y=Sample_lib18.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

V5.20190417.MN <- lmList(Sample_lib18.normalized.difference ~ Sample_lib17.normalized|species, data=canonical.chromosomes)
summary(V5.20190417.MN)
summary(V5.20190417.MN[[2]])

#V5 IP from MN differentiation 20190610 (MNs) (<1ng yield)
ggplot(combined.dataframe, aes(x=Sample_lib21.normalized, 
                               y=Sample_lib28.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

V5.20190610.MN <- lmList(Sample_lib28.normalized.difference ~ Sample_lib21.normalized|species, data=canonical.chromosomes)
summary(V5.20190610.MN)
summary(V5.20190610.MN[[2]])

#HA IP from MN differentiation 20190610 (MNs) 
ggplot(combined.dataframe, aes(x=Sample_lib23.normalized, 
                               y=Sample_lib30.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

HA.20190610.MN <- lmList(Sample_lib30.normalized.difference ~ Sample_lib23.normalized|species, data=canonical.chromosomes)
summary(HA.20190610.MN)
summary(HA.20190610.MN[[2]])

#V5 IP from MN differentiation 20190617 (MNs) (<1ng yield)
ggplot(combined.dataframe, aes(x=Sample_lib25.normalized, 
                               y=Sample_lib31.normalized.difference,
                               color=factor(species))) + 
  geom_point(alpha=0.6, size=3) +
  geom_point(data = subset(combined.dataframe, mito == "MT"), alpha = 0.6, color='black', size = 2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue',
                                'reporter'='coral'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38',
                                 'reporter'='GFP/mCherry')) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

V5.20190617.MN <- lmList(Sample_lib31.normalized.difference ~ Sample_lib25.normalized|species, data=canonical.chromosomes)
summary(V5.20190617.MN)
summary(V5.20190617.MN[[2]])


#function for returning values for each linear regression
get.lmlist.info <- function(x) {  
  regression.summary <- summary(x)
  slope <- regression.summary$coefficients[2,1]
  r.squared <- regression.summary$adj.r.squared
  p.value <- regression.summary$coefficients[2,4]
  #y <- list(slope, r.squared, p.value)
  y <- c(slope, r.squared, p.value)
  return(y)
}
regression.list <- list(HA.20190603.mAstro[[1]], 
                        V5.20190603.MN[[2]], 
                        HA.20190528.mAstro[[1]],
                        V5.20190528.MN[[2]],
                        HA.20190417.mAstro[[1]],
                        V5.20190417.MN[[2]],
                        V5.20190610.MN[[2]],
                        HA.20190610.MN[[2]],
                        V5.20190617.MN[[2]])

regression.values <- lapply(regression.list, get.lmlist.info) # call function to get desired values
dd  <-  data.frame(do.call(rbind, regression.values)) # convert list to dataframe
colnames(dd) <- c("slope", "r_squared", "pvalue")

# add the rest of the info from the meta data to the dataframe
libraries <- c(5,8,12,14,20,18,28,30,31)
library.names <- lapply(libraries, paste.names)
dd$library_name <- as.character(library.names)
ref.library <- c(2,2,10,10,17,17,21,23,25)
ref.library.names <- lapply(ref.library, paste.names) 
dd$input.library <- as.character(ref.library.names)
dd$mn.differentiation <- c(20190603,20190603,20190528,20190528,20190417,20190417,20190610,20190610,20190617)
dd$IP.cell.type <- c("mAstro","MN",	"mAstro",	"MN",	"mAstro",	"MN",	"MN",	"MN",	"MN")
dd$RiboTag <- c("HA",	"V5",	"HA",	"V5",	"HA",	"V5",	"V5",	"HA",	"V5")
dd$species <- c("mouse",	"human",	"mouse",	"human",	"mouse",	"human",	"human",	"human",	"human")
dd$depletion <- dd$slope*(-1)

# plot depletion efficiencies against cell type
ggplot(dd, aes(x=IP.cell.type, y=depletion, color=r_squared)) +
  geom_boxplot() +
  geom_jitter(alpha=0.6, size=4, width=0.3, aes(shape=RiboTag)) +
  scale_color_gradient(low = "blue", high = "red") +
  ylim(c(0,1)) +
  labs(y = "depletion efficiency", 
       x = "cell type") 

# plot depletion efficiencies across differentations
ggplot(dd, aes(x=as.factor(mn.differentiation), y=depletion, color=IP.cell.type)) +
  geom_point(alpha=0.6, size=4, aes(shape=RiboTag)) +
  scale_colour_manual(name = 'Cell Type', 
                      values =c('mAstro'='darkorchid1',
                                'MN'='deepskyblue'),
                      labels = c('mAstro'='mAstro',
                                 'MN'='motor neuron')) +
  ylim(c(0,1)) +
  labs(y = "depletion efficiency (%)", 
       x = "Replicate") +
  theme(axis.text.x = element_text(angle = 90)) 


