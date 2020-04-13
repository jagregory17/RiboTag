#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# Code to determine species distribution of RiboSeq from 3T3/HEK293 co-culture
###############################################

library(dplyr)
library(reshape2)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(gridExtra)
library(stringr)
library(plyr)
library(nlme) 
library(edgeR)
library(gsubfn)

#clear environment
rm(list=ls())
`%notin%` <- Negate(`%in%`)

project <- "path.to.working.dir"
setwd(project)

##################
# load data from samtools count data, which calculates the number of uniquely mapping reads at the chromosome level
################

inputcontrol.data.frame <- read.table(file = "path.to.file/barnyard.samtools.count.tab", sep = "\t", stringsAsFactors = F)
#normalize to the number of uniquely mapping reads

#normalize each sample to the number of uniquely mapping reads
for (x in 2:ncol(inputcontrol.data.frame)) {
  inputcontrol.data.frame[,paste(colnames(inputcontrol.data.frame[x]), ".normalized", sep = "")] <- inputcontrol.data.frame[x]/sum(inputcontrol.data.frame[x])
}

inputcontrol.data.frame$cell.species <- "mouse" # set mouse for all of them
inputcontrol.data.frame$cell.species[23:nrow(inputcontrol.data.frame)] <- "human" # replace appropriate rows with human
inputcontrol.data.frame$cell.species[48:nrow(inputcontrol.data.frame)] <- "Exo" # GFP and mCherry rows
normalized.sample.names <- c("Input.Control", "HA.IP2_human", "V5.IP2_mouse", "Mouse", "human", "HA.IP1_human", "V5.IP1_mouse")
colnames(inputcontrol.data.frame)[9:15] <- normalized.sample.names # use appropriate sample names

#Look at GFP and mCherry expression
exo.genes <- inputcontrol.data.frame[49:50,]
melt.exo.genes <- melt(exo.genes, id.vars = c("chromosome", "cell.species"), measure.vars = colnames(inputcontrol.data.frame[,c(12,13,9,14,11,15,10)])) #select them in the right order
#plot gfp/mcherry reads as a percentage of the uniquely mapping reads
ggplot(melt.exo.genes, aes(x=variable, y=value, fill = chromosome)) + 
  geom_col(position = "dodge") + 
  theme(axis.text.x = element_text(angle = 90)) + 
  scale_fill_manual(values=c("deepskyblue","orchid3"))

# extract gfp/mcherry reads, normalize to gene length, and calculate ratios
exo.genes <- exo.genes[,9:15] 
exo.genes <- data.frame(t(exo.genes))
colnames(exo.genes) <- c("mCherry", "GFP")
exo.genes$sample <- rownames(exo.genes)
exo.genes$mCherry <- exo.genes$mCherry/708 # normalize for transcript length
exo.genes$GFP <- exo.genes$GFP/717 # normalize for transcript legnth
exo.genes$ratio <- exo.genes$mCherry/exo.genes$GFP

# plot ratios for each sample
new.order <- c("human", "Mouse", "Input.Control", "HA.IP1_human", "V5.IP2_mouse", "V5.IP1_mouse", "HA.IP2_human") 
exo.genes$sample <- factor(exo.genes$sample ,levels = new.order)
exo.genes <- exo.genes[c(1:3,6:7),]
ggplot(exo.genes, aes(x=sample, y=ratio)) + 
  geom_col(position = "dodge", fill=(values=c("gold", "deepskyblue", "orchid3", "orchid3", "deepskyblue"))) + 
  theme_minimal() + 
  scale_y_continuous(trans='log2', breaks = trans_breaks("log2", function(x) 2^x), labels = trans_format("log2", math_format(2^.x))) +
  ylab(expression(paste("[", italic("mCherry"), "]/[", italic("GFP"), "]", sep = "")))

#calculate summary data
table.sum <- aggregate(test$value, by=list(species=test$cell.species, sampleIP=test$variable), FUN=sum)
table.sum$x <- as.integer(table.sum$x * 100) # reduce to integer
table.sum$sampleIP <- factor(table.sum$sampleIP ,levels = new.order) # put in the order I want them plotted

#plot summary data
table.sum.subset <- table.sum %>% filter(species %in% c("human", "mouse"))
ggplot(table.sum.subset, aes(x=sampleIP, y=x, fill=species, label=table.sum.subset$x)) +
  geom_col() +
  scale_fill_manual(values=c("deepskyblue","orchid3")) +
  geom_text(size = 6, position = position_stack(vjust = 0.5)) +
  theme_minimal()


#scatter plots of HA, input, and V5 uniquely mapping reads as a percentage of the total 
inputcontrol.data.frame$cell.species <- factor(inputcontrol.data.frame$cell.species, levels = c("mouse", "human")) #set factors and put them in the right order for plotting
inputcontrol.data.frame.exo <- inputcontrol.data.frame[49:50,]
inputcontrol.data.frame <- inputcontrol.data.frame[1:47,] # remove Exo genes from table
a <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = Input.Control, y = HA.IP1_human, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'human'),
             aes(x = Input.Control, y = HA.IP1_human, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

b <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = Input.Control, y = HA.IP2_human, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'human'),
             aes(x = Input.Control, y = HA.IP2_human, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

c <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = Input.Control, y = V5.IP1_mouse, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'mouse'),
             aes(x = Input.Control, y = V5.IP1_mouse, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

d <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = Input.Control, y = V5.IP2_mouse, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'mouse'),
             aes(x = Input.Control, y = V5.IP2_mouse, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

e <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = HA.IP1_human, y = V5.IP2_mouse, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'mouse'),
             aes(x = HA.IP1_human, y = V5.IP2_mouse, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

f <- ggplot(inputcontrol.data.frame) +
  geom_point(aes(x = HA.IP2_human, y = V5.IP1_mouse, color = cell.species)) +
  geom_point(data = subset(inputcontrol.data.frame, cell.species == 'mouse'),
             aes(x = HA.IP2_human, y = V5.IP1_mouse, color = cell.species)) +
  scale_colour_manual(values=c("orchid3", "deepskyblue"))

ggarrange(a, b, c, d, e, f, ncol = 2, nrow = 3, labels = c("A", "B", "C", "D", "E", "F"))

#make plots looking at enrichment using %uniquely mapped across chromosomes for mouse and human separately
mouse <- inputcontrol.data.frame %>% filter(cell.species == "mouse") 
mouse.normalized <- mouse[,c(9:15)]*100 # keep only percent uniquely mapping read value
# calculate the difference in % uniquely mapping reads between the IP and the input control
for (x in 1:ncol(mouse.normalized)) {
  mouse.normalized[,paste(colnames(mouse.normalized[x]), ".normalized", sep = "")] <- (mouse.normalized[x]-mouse.normalized[1])
}
mouse.normalized$chromosome <- mouse$chromosome

# same for human
human <- inputcontrol.data.frame %>% filter(cell.species == "human") 
human.normalized <- human[,c(9:15)]*100
for (x in 1:ncol(human.normalized)) {
  human.normalized[,paste(colnames(human.normalized[x]), ".normalized", sep = "")] <- (human.normalized[x]-human.normalized[1])
}
human.normalized$chromosome <- human$chromosome

# now for GFP & mcherry
exo.normalized <- inputcontrol.data.frame.exo[,c(9:15)]*100
for (x in 1:ncol(exo.normalized)) {
  exo.normalized[,paste(colnames(exo.normalized[x]), ".normalized", sep = "")] <- (exo.normalized[x] - exo.normalized[1])
}
exo.normalized$chromosome <- inputcontrol.data.frame.exo$chromosome

## merge the human and mouse dataframes
mouse.normalized$species <- "mouse"
human.normalized$species <- "human"
total.normalized <- rbind(mouse.normalized, human.normalized)

# plot HA samples only
HA.melt <- melt(total.normalized, id.vars = c("Input.Control", "species", "chromosome"), measure.vars = c("HA.IP1_human.normalized", "HA.IP2_human.normalized"))
HA.melt$label <- "a" # generically label all of them
HA.melt$label[c(47,94)] <- "MT" #specifically label the mito chromosome

b <- ggplot(HA.melt, aes(x=Input.Control, y=value, color=factor(species), shape=factor(variable), order=factor(species, levels = c("human", "mouse")), size=3)) + 
  geom_point(alpha=0.6) +
  geom_point(data = subset(HA.melt, label == "MT"), alpha = 0.6, size=2, color='black') +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38')) +
  scale_shape_manual(name = "IP order",
                     values = c("HA.IP1_human.normalized" = 1,
                                "HA.IP2_human.normalized" = 16),
                     labels = c('1st', '2nd')) +
  scale_x_continuous(limits=c(0, 11)) +
  scale_y_continuous(limits=c(-10, 10)) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

# plot V5 samples only
V5.melt <- melt(total.normalized, id.vars = c("Input.Control", "species", "chromosome"), measure.vars = c("V5.IP1_mouse.normalized", "V5.IP2_mouse.normalized"))
V5.melt$label <- "a" # generically label all of them
V5.melt$label[c(22,69)] <- "MT" #specifically label the mito chromosome


a <- ggplot(V5.melt, aes(x=Input.Control, y=value, color=factor(species), shape=factor(variable), order=factor(species, levels = c("human", "mouse")), size=3)) + 
  geom_point(alpha=0.6) +
  geom_point(data = subset(V5.melt, label == "MT"), alpha = 0.6, color='black', size=2) +
  scale_colour_manual(name = 'Genome', 
                      values =c('mouse'='darkorchid1',
                                'human'='deepskyblue'), 
                      labels = c('mouse'='mm20',
                                 'human'='hg38')) +
  scale_shape_manual(name = "IP order",
                     values = c("V5.IP1_mouse.normalized" = 1,
                                "V5.IP2_mouse.normalized" = 16),
                     labels = c('1st', '2nd')) +
  scale_x_continuous(limits=c(0, 11)) +
  scale_y_continuous(limits=c(-10, 10)) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

ggarrange(a, b)

# calculate the slopes of each line; off.target slopes are the depletion efficiency values
a <- lmList(V5.IP1_mouse.normalized ~ Input.Control|species, data=total.normalized) # splits data points by species
b <- lmList(V5.IP2_mouse.normalized ~ Input.Control|species, data=total.normalized)
c <- lmList(HA.IP1_human.normalized ~ Input.Control|species, data=total.normalized)
d <- lmList(HA.IP2_human.normalized ~ Input.Control|species, data=total.normalized)
summary(a)
summary(b)
summary(c)
summary(d)

# calculate on-target slopes after removing mitochondrial genes
# do this again without the mitochondrial genes
mouse.normalized.no.mito <- mouse.normalized %>% filter(chromosome %notin% "Mm_ChrM")
mouse.lm.V5IP1.ontarget.no.mito <- lm(V5.IP1_mouse.normalized ~ Input.Control, data=mouse.normalized.no.mito) 
mouse.lm.V5IP2.ontarget.no.mito <- lm(V5.IP2_mouse.normalized ~ Input.Control, data=mouse.normalized.no.mito)
summary(mouse.lm.V5IP1.ontarget.no.mito)
summary(mouse.lm.V5IP2.ontarget.no.mito)
human.normalized.no.mito <- human.normalized %>% filter(chromosome %notin% "Hs_ChrM")
human.lm.HAIP1.ontarget.no.mito <- lm(HA.IP1_human.normalized ~ Input.Control, data=human.normalized.no.mito)
human.lm.HAIP2.ontarget.no.mito <- lm(HA.IP2_human.normalized ~ Input.Control, data=human.normalized.no.mito)
summary(human.lm.HAIP1.ontarget.no.mito)
summary(human.lm.HAIP2.ontarget.no.mito)

# calculate predicted % of total library based on regressions from off-target genes
V5.human.depletion <- c("V5-1"=a$human$coefficients[2], "V5-2"=b$human$coefficients[2])
HA.mouse.depletion <- c("HA-1"=c$mouse$coefficients[2], "HA-2"=d$mouse$coefficients[2])

melt.predictor <- melt(exo.normalized, id.vars = c("chromosome", "Input.Control"), 
                       measure.vars = c("V5.IP1_mouse","V5.IP2_mouse", "HA.IP1_human", "HA.IP2_human"))
melt.predictor.off.target <- melt.predictor[c(2,4,5,7),]
melt.predictor.off.target$slope <- c(a$human$coefficients[2], 
                                     b$human$coefficients[2], 
                                     c$mouse$coefficients[2], 
                                     d$mouse$coefficients[2])
melt.predictor.off.target$predicted <- melt.predictor.off.target$Input.Control*(1+melt.predictor.off.target$slope) # calculate predicted value
melt.again <- melt(melt.predictor.off.target, id.vars = c("chromosome", "variable"), measure.vars = c("value", "predicted"))
colnames(melt.again) <- c("chromosome", "IP", "variable", "value")
melt.again$variable <- gsub("value", "measured", melt.again$variable)
ggplot(melt.again, aes(x=variable, y=value, color=chromosome)) +
  geom_point(size = 3) + 
  scale_y_continuous(limits=c(0, .01)) +
  labs(y = expression("% uniquely mapping reads")) +
  scale_colour_manual(name = "Sample",
                      values =c('mCherry'='darkorchid1',
                                'GFP'='deepskyblue')) +
  geom_line(aes(group = IP), color="black") +
  theme_minimal()
  

# look at gene level counts data from featureCounts
counts <- read.table(file = "/gpfs/commons/home/jgregory/Barnyard/barnyard.count.matrix.tab", header = T, sep = "\t")
isexpr <- rowSums(cpm(counts[,c(3:9)])>1) >= 3 # this keeps all genes with at least 0.5 cpm in at least 3 of the samples.
sum(isexpr, na.rm = TRUE) 
counts <- counts[isexpr,] # filter low genes

# add biomart data
human.geneInfo.biomaRt <- read.table(file = file.path("/gpfs/commons/home/jgregory/MN_astro_coCulture/20190918", "human.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F) # load biomaRt info from saved table instead of biomaRt; faster
mouse.geneInfo.biomaRt <- read.table(file = file.path("/gpfs/commons/home/jgregory/MN_astro_coCulture/20190918", "mouse.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F) # load biomaRt info from saved table instead of biomaRt; faster
human.geneInfo.biomaRt$chrom <- paste("Hs_chr", human.geneInfo.biomaRt$chromosome_name, sep = "")
mouse.geneInfo.biomaRt$chrom <- paste("Mm_chr", mouse.geneInfo.biomaRt$chromosome_name, sep = "")
#harmonize colnames for the mouse and human biomart dataframes
cols.to.use <- intersect(colnames(human.geneInfo.biomaRt), colnames(mouse.geneInfo.biomaRt)) # which columns overlap
combine.biomart.mouse <- mouse.geneInfo.biomaRt %>% dplyr::select(cols.to.use)
combine.biomart.human <- human.geneInfo.biomaRt %>% dplyr::select(cols.to.use)
total.biomart <- rbind(combine.biomart.human, combine.biomart.mouse) # combine using overlapping columns
counts.biomart <- counts # make a new dataframe for adding biomart info 
counts.biomart$ensembl_gene_id <- gsub("\\..*","", counts.biomart$Geneid) #remove version number
counts.biomart$species <- grepl("ENSG", counts.biomart$ensembl_gene_id)  # set species as TRUE for human and FALSE for mouse
genes.with.annotations <- intersect(counts.biomart$ensembl_gene_id, total.biomart$ensembl_gene_id) 
genes.without.annotations <- setdiff(counts.biomart$ensembl_gene_id, total.biomart$ensembl_gene_id)
counts.biomart <- inner_join(total.biomart, counts.biomart) # 14809 genes left with annotations


meltraw.counts.biomart <- melt(counts.biomart, id.vars = c("gene_biotype", "species", "chrom"), measure.vars = c("R1_input_control",
                                                                                                                 "R8_IP_V5", 
                                                                                                                 "R12_IP2_V5", 
                                                                                                                 "R11_IP2_HA",
                                                                                                                 "R7_IP1_HA"))
                                                                                                                 
#plot protein coding vs all other biotypes
biotypes <- unique(meltraw.counts.biomart$gene_biotype)
for (x in 2:length(biotypes)) {
  meltraw.counts.biomart$gene_biotype <- gsub(biotypes[x], "other", meltraw.counts.biomart$gene_biotype, fixed = T)
}
biotypes <- unique(meltraw.counts.biomart$gene_biotype) # first round didn't change all of them
for (x in 3:length(biotypes)) {
  meltraw.counts.biomart$gene_biotype <- gsub(biotypes[x], "other", meltraw.counts.biomart$gene_biotype)
}
unique(meltraw.counts.biomart$gene_biotype)
ggplot(meltraw.counts.biomart, aes(x=variable, y=value, fill=gene_biotype)) + geom_col(position="fill") + facet_grid(~species)

#summarize data
table.sum.raw <- aggregate(meltraw.counts.biomart$value, by=list(biotype=meltraw.counts.biomart$gene_biotype,
                                                                 species=meltraw.counts.biomart$species,
                                                                 sampleIP=meltraw.counts.biomart$variable), 
                           FUN=sum)
######################
#calculate percentage of library for each gene in the human library
######################

human.counts.biomart <- counts.biomart %>% filter(species=="TRUE")
for (x in 9:15) {
  human.counts.biomart[,paste(colnames(human.counts.biomart[x]), ".normalized", sep = "")] <- human.counts.biomart[x]/sum(human.counts.biomart[x])
}

#calculate change (IP-input control)
for (x in 18:23) {
  human.counts.biomart[,paste(colnames(human.counts.biomart[x]), ".difference", sep = "")] <- (human.counts.biomart[x]-human.counts.biomart[17])
}


# plot protein coding vs everything else for HA on human
melt.counts.biomart <- melt(human.counts.biomart, id.vars = c("R1_input_control.normalized", 
                                                              "gene_biotype", 
                                                              "species", 
                                                              "chrom", 
                                                              "ensembl_gene_id"), 
                            measure.vars = c("R11_IP2_HA.normalized.difference",
                                             "R7_IP1_HA.normalized.difference"))
#melt.counts.biomart <- melt.counts.biomart %>% filter(chrom %notin% c("Hs_chrMT")) #remove MT genes
table.sum <- aggregate(melt.counts.biomart$value, by=list(species=melt.counts.biomart$species, 
                                                          biotype=melt.counts.biomart$gene_biotype, 
                                                          sampleIP=melt.counts.biomart$variable), 
                       FUN=sum)
protein.coding.melt.counts.biomart <- melt.counts.biomart %>% filter(gene_biotype %in% c("protein_coding")) # extract protein coding
other.melt.counts.biomart <- melt.counts.biomart %>% filter(gene_biotype %notin% c("protein_coding")) # extract everything else
other.melt.counts.biomart$gene_biotype <- "other" #change biotype to other for everythign except protein coding
melt.counts.biomart <- rbind(protein.coding.melt.counts.biomart, other.melt.counts.biomart) # rejoin with changed biotype
melt.counts.biomart <- melt.counts.biomart %>% filter(species=="TRUE") # get only human
melt.counts.biomart$value <- melt.counts.biomart$value*100 # change ot percentage
melt.counts.biomart$chrom_label <- grepl("Hs_chrMT", melt.counts.biomart$chrom)

HA.plot <- ggplot(melt.counts.biomart, aes(x=R1_input_control.normalized, y=value, color=factor(gene_biotype), shape=factor(variable), size=3)) + 
  geom_point(alpha=0.4) +
  geom_point(data = subset(melt.counts.biomart, chrom_label == "TRUE"), alpha = 0.6, color="black", size = 2) +
  scale_colour_manual(name = 'gene biotype', 
                      values =c('other'='deeppink1',
                                'protein_coding'='cyan3',
                                'black'='black'),
                      labels = c('other', 'protein_coding', "mitochondrial")) +
  scale_shape_manual(name = "IP order",
                     values = c("R7_IP1_HA.normalized.difference" = 16,
                                "R11_IP2_HA.normalized.difference" = 17),
                     labels = c("1st", "2nd")) +
  #scale_x_continuous(limits=c(0, .02)) +
  scale_y_continuous(breaks=c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1)) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

#enriched genes
human.counts.biomart$avg <- (human.counts.biomart$R7_IP1_HA.normalized.difference + human.counts.biomart$R11_IP2_HA.normalized.difference)/2 # avg difference across two replicates
human.enriched.genes <- human.counts.biomart %>% filter(avg > 0)
NROW(human.enriched.genes)
table(human.enriched.genes$gene_biotype)
human.depleted.genes <- human.counts.biomart %>% filter(avg < 0)
NROW(human.depleted.genes)
table(human.depleted.genes$gene_biotype)

######################
#calculate percentage of library for each gene in the mouse library for V5 samples
######################

mouse.counts.biomart <- counts.biomart %>% filter(species=="FALSE")
for (x in 9:15) {
  mouse.counts.biomart[,paste(colnames(mouse.counts.biomart[x]), ".normalized", sep = "")] <- mouse.counts.biomart[x]/sum(mouse.counts.biomart[x])
}

#calculate change (IP-input control)
for (x in 18:23) {
  mouse.counts.biomart[,paste(colnames(mouse.counts.biomart[x]), ".difference", sep = "")] <- (mouse.counts.biomart[x]-mouse.counts.biomart[17])
}


# plot protein coding vs everything else for HA on human
melt.counts.biomart <- melt(mouse.counts.biomart, id.vars = c("R1_input_control.normalized", 
                                                              "gene_biotype", 
                                                              "species", 
                                                              "chrom", 
                                                              "ensembl_gene_id"), 
                            measure.vars = c("R8_IP_V5.normalized.difference",
                                             "R12_IP2_V5.normalized.difference"))
#melt.counts.biomart <- melt.counts.biomart %>% filter(chrom %notin% c("Hs_chrMT")) #remove MT genes
table.sum <- aggregate(melt.counts.biomart$value, by=list(species=melt.counts.biomart$species, 
                                                          biotype=melt.counts.biomart$gene_biotype, 
                                                          sampleIP=melt.counts.biomart$variable), 
                       FUN=sum)

protein.coding.melt.counts.biomart <- melt.counts.biomart %>% filter(gene_biotype %in% c("protein_coding")) # extract protein coding
other.melt.counts.biomart <- melt.counts.biomart %>% filter(gene_biotype %notin% c("protein_coding")) # extract everything else
other.melt.counts.biomart$gene_biotype <- "other" #change biotype to other for everythign except protein coding
melt.counts.biomart <- rbind(protein.coding.melt.counts.biomart, other.melt.counts.biomart) # rejoin with changed biotype
melt.counts.biomart <- melt.counts.biomart %>% filter(species=="FALSE") # get only mouse
melt.counts.biomart$value <- melt.counts.biomart$value*100 # change ot percentage
melt.counts.biomart$chrom_label <- grepl("Mm_chrMT", melt.counts.biomart$chrom)

V5.plot <- ggplot(melt.counts.biomart, aes(x=R1_input_control.normalized, y=value, color=factor(gene_biotype), shape=factor(variable), size = 3)) + 
  geom_point(alpha=0.4) +
  geom_point(data = subset(melt.counts.biomart, chrom_label == "TRUE"), alpha = 0.6, color="black", size = 2) +
  scale_colour_manual(name = 'gene biotype', 
                      values =c('other'='deeppink1',
                                'protein_coding'='cyan3',
                                'black'='black'),
                      labels = c('other', 'protein_coding', "mitochondrial")) +
  scale_shape_manual(name = "IP order",
                     values = c("R8_IP_V5.normalized.difference" = 16,
                                "R12_IP2_V5.normalized.difference" = 17),
                     labels = c("1st", "2nd")) +
  #scale_x_continuous(limits=c(0, .03)) +
  scale_y_continuous(breaks=c(-3,-2.5,-2,-1.5,-1,-0.5,0,0.5,1), limits = c(-2.5,0.5)) +
  theme_minimal() +
  labs(x = "% uniquely mapped reads: Input", 
       y = expression(Delta* "% uniquely mapping reads (IP - Input)")) +
  geom_hline(yintercept = 0, linetype = 'dotted', col = 'red')

#enriched genes
mouse.counts.biomart$avg <- (mouse.counts.biomart$R8_IP_V5.normalized.difference + mouse.counts.biomart$R12_IP2_V5.normalized.difference)/2 # avg difference across two replicates
mouse.enriched.genes <- mouse.counts.biomart %>% filter(avg > 0)
NROW(mouse.enriched.genes)
table(mouse.enriched.genes$gene_biotype)
mouse.depleted.genes <- mouse.counts.biomart %>% filter(avg < 0)
NROW(mouse.depleted.genes)
table(mouse.depleted.genes$gene_biotype)

ggarrange(HA.plot, V5.plot)
