#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# calculate species composition of coculture inputs and matched IPs
# from STAR (gene level unique reads), samtools (all unique reads), and protein coding (ensembl gene annotations
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
samtools.count.matrix <- samtools.count.matrix %>% select(all_of(matched.cocultures))
count.matrix <- count.matrix %>% select(all_of(c(matched.cocultures, "ensembl_gene_id", "species")))

#########################
## calculate the species composition of each 
## library using the counts from samtools and STAR
#############################
samtools.count.matrix <- samtools.count.matrix[,order(names(samtools.count.matrix))] #reorder the column names by library name
count.matrix <- count.matrix[,order(names(count.matrix))] #reorder the column names by Run 
colnames(samtools.count.matrix) == meta.data.subset$library_name #check that order matches
colnames(count.matrix[2:16]) == meta.data.subset$library_name # check order

# add column for species
samtools.count.matrix$cell.species <- "mouse" # set mouse for all of them
samtools.count.matrix$cell.species[23:nrow(samtools.count.matrix)] <- "human" # replace appropriate rows with human


# calculate percent human for each sample using STAR or samtools counts
human.reads <- samtools.count.matrix %>% filter(cell.species == "human")
mouse.reads <- samtools.count.matrix %>% filter(cell.species == "mouse")
star.human.reads <- count.matrix %>% filter(species == "TRUE")
star.mouse.reads <- count.matrix %>% filter(species == "FALSE")
library.sizes <- data.frame("samtools" = colSums(samtools.count.matrix[,1:15]), 
                            "STAR" = colSums(count.matrix[2:16]),
                            "samtools.percent.human" = (colSums(human.reads[,1:15])/colSums(samtools.count.matrix[1:47,1:15])),
                            "samtools.percent.mouse" = (colSums(mouse.reads[,1:15])/colSums(samtools.count.matrix[1:47,1:15])),
                            "library_name" = colnames(samtools.count.matrix)[1:15],
                            "star.percent.human" = (colSums(star.human.reads[,2:16])/colSums(count.matrix[1:113473,2:16])),
                            "star.percent.mouse" = (colSums(star.mouse.reads[,2:16])/colSums(count.matrix[1:113473,2:16]))) 

melt.library.sizes <- melt(library.sizes, id.vars = c("library_name"), 
                           measure.vars = c("samtools.percent.human", "samtools.percent.mouse")) 
a <- ggplot(melt.library.sizes, aes(x=factor(library_name, levels = matched.cocultures), y=value, fill = variable)) + 
  geom_col() +
  scale_fill_manual(name="species",
                    values=c("deepskyblue","orchid3"),
                    labels=c("human", "mouse")) +
  theme_minimal() +
  scale_x_discrete(labels = c("co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", 
                              "co-culture", "MN IP",
                              "co-culture", "MN IP")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("% uniquely mapped reads") +
  labs(title="Samtools")

melt.library.sizes <- melt(library.sizes, id.vars = c("library_name"), 
                           measure.vars = c("star.percent.human", "star.percent.mouse")) 

b <- ggplot(melt.library.sizes, aes(x=factor(library_name, levels = matched.cocultures), y=value, fill = variable)) + 
  geom_col() +
  scale_fill_manual(name="species",
                    values=c("deepskyblue","orchid3"),
                    labels=c("human", "mouse")) +
  theme_minimal() +
  scale_x_discrete(labels = c("co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", 
                              "co-culture", "MN IP",
                              "co-culture", "MN IP")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("% uniquely mapped reads") +
  labs(title="STAR")

ggarrange(a,b)

##################################
## make this plot again looking at 
## only protein coding genes not on the MT genome
#################################

# merge with annotation file
mouse.geneInfo.biomaRt <- read.table(file = file.path(project, "mouse.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F)
mouse.geneInfo.biomaRt <- mouse.geneInfo.biomaRt %>% select(-ensembl_gene_id_version) # get rid of this column
human.geneInfo.biomaRt <- read.table(file = file.path("/gpfs/commons/home/jgregory/MN_astro_coCulture/20190918", "human.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F) # load biomaRt info from saved table instead of biomaRt; faster

joined <- left_join(star.human.reads, human.geneInfo.biomaRt) # adding biomart info where available
dropped.genes <- setdiff(star.human.reads$ensembl_gene_id, joined$ensembl_gene_id) # get a list of the genes that are dropped, if any.

unique(joined$gene_biotype) # how many biotypes are there
genes.to.keep <- c("protein_coding")
human.filtered <- subset(joined, (gene_biotype %in% genes.to.keep)) # keep protein coding genes
human.filtered <- subset(human.filtered, !(chromosome_name %in% "MT")) # remove by location
human.filtered <- human.filtered %>% select(colnames(human.filtered)[1:17], external_gene_name, description, chromosome_name, gene_biotype)

# join mouse reads with biomaRt
mouse.joined <- left_join(star.mouse.reads, mouse.geneInfo.biomaRt)
mouse.dropped.genes <- setdiff(star.mouse.reads$ensembl_gene_id, mouse.joined$ensembl_gene_id)
mouse.filtered <- subset(mouse.joined, (gene_biotype %in% genes.to.keep))
mouse.filtered <- subset(mouse.filtered, !(chromosome_name %in% "MT"))

# make a master list will all protein coding genes
total.protein.coding <- rbind(human.filtered, mouse.filtered)

# look for duplicates and get rid of them
total.protein.coding$duplicate <- duplicated(total.protein.coding$ensembl_gene_id)
total.duplicated <- total.protein.coding %>% filter(duplicate == "TRUE")
test <- subset(total.protein.coding, (ensembl_gene_id %in% total.duplicated$ensembl_gene_id)) 
total.protein.coding <- total.protein.coding %>% filter(duplicate == "FALSE") # all duplicates are human and actual duplicates, can get rid of them

# compute %human reads using the protein coding genes and add it to the meta data
human.total <- colSums(human.filtered[2:16])
mouse.total <- colSums(mouse.filtered[2:16])
total.reads <- human.total + mouse.total
percent.human <- (human.total/total.reads)*100
percent.mouse <- (mouse.total/total.reads)*100

meta.data.subset$library_name == colnames(human.filtered)[2:16] # check that order is correct before adding to meta data
meta.data.subset$percent.human <- percent.human
meta.data.subset$percent.mouse <- percent.mouse
meta.data.subset$mouse.reads <- mouse.total
meta.data.subset$human.reads <- human.total

# plot 
melt.meta.percent.human.mouse <- melt(meta.data.subset, 
                                      id.vars = c("library_name"), 
                                      measure.vars = c("percent.human", "percent.mouse"))
c <- ggplot(melt.meta.percent.human.mouse, 
       aes(x=factor(library_name, levels = matched.cocultures), 
           y=value, fill=variable)) +
  geom_col() +
  scale_fill_manual(name="species",
                    values=c("deepskyblue","orchid3"),
                    labels=c("human", "mouse")) +
  theme_minimal() +
  scale_x_discrete(labels = c("co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", "mAstro IP",
                              "co-culture", "MN IP", 
                              "co-culture", "MN IP",
                              "co-culture", "MN IP")) +
  theme(axis.text.x = element_text(angle = 90)) +
  ylab("% uniquely mapped reads") +
  labs(title = "protein coding")

ggarrange(a,b,c)
