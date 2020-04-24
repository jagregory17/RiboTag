#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# enrichment analysis of MN & mouse astro co-cultures
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
library(cowplot)
library(lme4)
library(edgeR)

#clear environment
rm(list=ls())
project <- "path.to.working.dir" 

#functions
'%!in%' <- function(x,y)!('%in%'(x,y))  # create a function that is the opposite of %in%; see this URL https://www.r-bloggers.com/the-notin-operator/
`%notin%` <- Negate(`%in%`) # another way for a 'not' in function

#load meta data
meta.data <- read.table(file = "path.to.file/20200422_MNastro_meta.data.tab", 
                        sep = "\t", 
                        stringsAsFactors = F,
                        header = T)
meta.data.subset <- meta.data %>% filter(Cell.Type %notin% c("GABA", "NPC")) # keep only the motor neuron/mouse astro samples
#load gene count matrix
count.matrix <- read.table(file = file.path(project, "20200422.count.matrix.hg38.mm20.tab"), 
                           header = T, 
                           sep = "\t", 
                           stringsAsFactors = F)

# vector of matched libraries by row (co.culture, MN IP, mAstro IP)
matched.cocultures <- c("Sample_lib2", "Sample_lib5", "Sample_lib8",
                        "Sample_lib10", "Sample_lib12", "Sample_lib14",
                        "Sample_lib17", "Sample_lib20", "Sample_lib18", 
                        "Sample_lib21", "Sample_lib28",
                        "Sample_lib23", "Sample_lib30", 
                        "Sample_lib25", "Sample_lib31")
#################################
# plot expression in RPKM for mCherry/GFP
###################################

# calculate RPKM of reporter genes based on cell type, not total culture
# so we can compare expression in each cell type across samples
# most mouse cells are GFP (HA) and most human cells are mCherry (V5); 
# calculate RPKM for both reporter genes for all samples in both cell types,
# then define the on-target cell type

human.reads <- count.matrix %>% filter(species == "TRUE")
mouse.reads <- count.matrix %>% filter(species == "FALSE")
GFP.counts <- count.matrix %>% filter(ensembl_gene_id == "GFP")
mouse.gfp.counts <- rbind(mouse.reads[,c(1:40)], GFP.counts) # add to mouse genes
mouse.gfp.rpkm <- data.frame(rpkm(mouse.gfp.counts[,c(1:36)], gene.length = mouse.gfp.counts$Length, log = F))
rownames(mouse.gfp.rpkm) <- mouse.gfp.counts$ensembl_gene_id

mouse.gfp.rpkm <- mouse.gfp.rpkm %>% select(meta.data.subset$library_name) 
mouse.gfp.rpkm <- mouse.gfp.rpkm[c(NROW(mouse.gfp.rpkm)),] # get just the rpkm for gfp

mCherry.counts <- count.matrix %>% filter(ensembl_gene_id == "mCherry")
human.mcherry.counts <- rbind(human.reads[,c(1:40)], mCherry.counts) # add to human genes
human.mcherry.rpkm <- data.frame(rpkm(human.mcherry.counts[,c(1:36)], gene.length = human.mcherry.counts$Length, log = F))
human.mcherry.rpkm <- human.mcherry.rpkm %>% select(meta.data.subset$library_name)
human.mcherry.rpkm <- human.mcherry.rpkm[c(NROW(human.mcherry.rpkm)),]

# there is also a single MN samples that used RPL22-HA, so need to calculate GFP rpkm also
human.gfp.counts <- rbind(human.reads[,c(1:40)], GFP.counts)
human.gfp.rpkm <- data.frame(rpkm(human.gfp.counts[,c(1:36)], gene.length = human.gfp.counts$Length, log = F))
human.gfp.rpkm <- human.gfp.rpkm %>% select(meta.data.subset$library_name)
human.gfp.rpkm <- human.gfp.rpkm[c(NROW(human.gfp.rpkm)),]

# subset metadata to add new rpkm values for gfp/mcherry that are species specific
meta.data.subset$library_name == names(human.gfp.rpkm) #make sure they are in the right order
meta.data.subset$human.gfp.rpkm <- as.numeric(t(human.gfp.rpkm))
meta.data.subset$human.mcherry.rpkm <- as.numeric(t(human.mcherry.rpkm))
meta.data.subset$mouse.gfp.rpkm <- as.numeric(t(mouse.gfp.rpkm))
meta.data.subset$celltype.antibody <- paste(meta.data.subset$Cell.Type, 
                                            meta.data.subset$Antibody,
                                            sep = ".")


melt.all.libraries <- melt(meta.data.subset, id.vars = c("library_name", "Group", "Antibody", "celltype.antibody", "paired.samples"),
                           measure.vars = c("human.gfp.rpkm", "human.mcherry.rpkm", "mouse.gfp.rpkm"))

HA.motor.neuron <- "Motor Neuron.HA" # for subsetting the motor neurons that expressed RPL22-HA instead of RPL22-V5
MN.HA.libs <- c("Sample_lib21", "Sample_lib22", "Sample_lib28", "Sample_lib29", "Sample_lib36") # these are the specific libraries of MN's expressing RPL22-HA

# retrive RPL22-HA MN samples
MN.HA.melt <- melt.all.libraries %>% filter(library_name %in% MN.HA.libs & variable == "human.gfp.rpkm")
MN.HA.melt$RiboTag <- "RPL22-HA" # add an extra column for which RiboTag is being expressed

# get rid of all human.gfp (e.g. RPL22-HA) and libs from the MN.HA samples. leaving us with the remaining on target measurements
remaining.samples <- melt.all.libraries %>% 
  filter(variable != "human.gfp.rpkm") %>%
  filter(library_name %notin% MN.HA.libs)

# from these, get all on target mCherry samples (e.g. remove astro only samples)
# so few human (i.e. off-target reads) that even a single read in mcherry looks high in the rpkm
mcherry.on.target <- remaining.samples %>% 
  filter(variable == "human.mcherry.rpkm") %>%
  filter(Group %notin% c("astro.input", "astro.IP"))
mcherry.on.target$RiboTag <- "RPL22-V5"

# do the analogous for on target GFP reads
gfp.on.target <- remaining.samples %>% 
  filter(variable == "mouse.gfp.rpkm") %>%
  filter(Group %notin% c("motor.neuron.input", "motor.neuron.IP"))
gfp.on.target$RiboTag <- "RPL22-HA"
merge.on.target <- rbind(MN.HA.melt, mcherry.on.target, gfp.on.target)
merge.on.target <- merge.on.target %>% filter(library_name != "Sample_lib27") # removing an astro input that was not transduced
merge.on.target <- merge.on.target %>% filter(value != 0) # remaining samples with 0's are where there is no gfp.rpkm because the astros were not transduced

# harmonize naming for on target reporter gene in MNs
merge.on.target$variable <- gsub("human.gfp.rpkm", "human.on.target", merge.on.target$variable)
merge.on.target$variable <- gsub("human.mcherry.rpkm", "human.on.target", merge.on.target$variable)
merge.on.target$variable <- gsub("mouse.gfp.rpkm", "mouse.on.target", merge.on.target$variable)

merge.on.target$Group <- factor(merge.on.target$Group, levels = c("astro.input",
                                                                  "astro.IP",
                                                                  "co.culture.astro.IP",
                                                                  "co.culture.input",
                                                                  "co.culture.motor.neuron.IP",
                                                                  "motor.neuron.IP",
                                                                  "motor.neuron.input"))
ggplot(merge.on.target, aes(x=Group, y=(value))) +
  #geom_boxplot(color="black") +
  geom_boxplot(aes(color=factor(variable))) +
  geom_jitter(width = 0.2, size=3, alpha=0.6, aes(shape=factor(RiboTag))) +
  scale_color_manual(values= c("deepskyblue", "orchid")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_y_log10(limits=c(1,100000)) +
  labs(y="RPKM")

ggplot(merge.on.target, aes(x=RiboTag, y=value, color=variable)) +
  geom_boxplot(color = "black") +
  geom_jitter(width = 0.2, size =3, alpha=0.6) +
  scale_color_manual(name="Cell Type",
                     values= c("deepskyblue", "orchid"),
                     labels=c("hiPSC-MN","mAstro")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_y_log10(limits=c(1,100000)) +
  labs(y="RPKM (GFP or mCherry)")

ggplot(merge.on.target, aes(x=variable, y=value, color=RiboTag)) +
  geom_boxplot(color = "black") +
  geom_jitter(width = 0.2, size =3, alpha=0.6) +
  scale_x_discrete(name="cell type",
                   labels=c("hiPSC-MN", "mAstro")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_y_log10(limits=c(1,100000)) +
  labs(y="RPKM (GFP or mCherry)")





