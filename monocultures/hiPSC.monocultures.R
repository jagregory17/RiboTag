#!/usr/bin/env Rscript

###############################################
#
# James Gregory (jgregory@nygenome.org), 
# analyze input and IP samples from monocultures
###############################################

library('limma')
library('edgeR')
library('biomaRt')
library('dplyr')
library('DESeq2')
library("ggplot2")
library("pheatmap")
library("tidyr")
library("reshape2")
library("tibble")


#clear environment
rm(list=ls())

# read in count data; monoculture reads aligned to hg38
counts <- read.table(file = "path.to.file/20200421.monoculture.counts.tab", stringsAsFactors = F, sep = "\t")
counts <- counts[,-1] # remove gene length
meta.data <- read.table(file = "path.to.file/20200422_MNastro_meta.data.tab", sep = "\t", stringsAsFactors = F)
colnames(meta.data) <- meta.data[1,] # make column names
meta.data <- meta.data[-1,1:24] # select desired rows/columns

# select columns corresponding to samples used for this analysis; sample names are found in counts file
meta.data <- meta.data %>% filter(library_name %in% colnames(counts)) # select needed meta data
rownames(meta.data) <- meta.data$library_name # add back rownames
counts <- counts %>% dplyr::select(meta.data$library_name) # reorder

# check ordering of counts and meta data
rownames(meta.data) == colnames(counts) 
counts.only.ordered <- counts

#create DESeq object to compare 
dds <- DESeqDataSetFromMatrix(countData=counts.only.ordered, colData=meta.data, design=~RNA_bin + Sample_Type) 
keep <- rowSums(counts(dds)) >= 100 #filtering criteria, same as recommended in DESeq vignette; creates a logical
sum(keep, na.rm = T)
dds <- dds[keep,] #get rid of genes with very low counts
dds <- estimateSizeFactors(dds) # use defualt method

#make a heatmap of the top 5000 genes; this plot looks junky since it hasn't been normlized
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:1000]
df <- as.data.frame(colData(dds)[,c("Sample_Type", "Cell.Type", "Antibody")])

#calculate a variance stabilized normalization
vst <- vst(dds, blind=TRUE) # faster version of rld and almost equivalent
vst.assay <- as.data.frame(assay(vst)) #extract the vst normalized counts
pheatmap(assay(vst)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)
plotPCA(vst, intgroup=c("Cell.Type"))
plotPCA(vst, intgroup=c("Sample_Type"))
#########################
## Differential expression
##########################

dds.multicell <- DESeq(dds)

# make contrast groups
resultsNames(dds.multicell)
TSPvsIP <- results(dds.multicell, contrast=c("Sample_Type","IP", "Input")) # set contrast
summary(TSPvsIP)
TSPvsIP.results <- as.data.frame(TSPvsIP) # make a data frame from the DESeq object specifying the contrast
TSPvsIP.results$ensembl_gene_id <- gsub("\\..*","", rownames(TSPvsIP.results)) #get rid of the version number
rownames(TSPvsIP.results) <- TSPvsIP.results$ensembl_gene_id # make rownames without version number
human.geneInfo.biomaRt <- read.table(file = "path.to.file/human.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F) # load biomaRt annotations
annotated.TSPvsIP <- left_join(TSPvsIP.results, human.geneInfo.biomaRt) # merge with annotations by ensembl gene id

###########################
# GSEA analysis looking at gene biotype enrichment
############################
library(fgsea)
TSPvsIP.results.filter <- annotated.TSPvsIP %>% filter(padj < 0.1) # p-value filter

# make a list of genes for the chromosomes and use this instead of pathways
genes.by.biotype <- list()
unique.biotypes <- unique(human.geneInfo.biomaRt$gene_biotype) # get list of chromosome names
for (i in 1:length(unique.biotypes)) {
  x <- human.geneInfo.biomaRt %>% filter(gene_biotype == unique.biotypes[i]) %>% dplyr::select(c("ensembl_gene_id"))
  x <- as.character(x[,c(1)])
  genes.by.biotype[[unique.biotypes[i]]] <- x
}

# rank genes
res2 <- TSPvsIP.results.filter %>% 
  dplyr::select(ensembl_gene_id, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>%
  summarize(stat=mean(stat))
res2

ranks <- deframe(res2) # convert to named vector
barplot(sort(ranks, decreasing = T)) # make a plot of the ranks
head(ranks, 20)

# run pre-ranked gene set enrichment analysis
fgseaRes <- fgsea(pathways=genes.by.biotype, stats=ranks, nperm=10000, minSize = 10)

# organize
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()

# plot
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  scale_fill_manual(values=c("deepskyblue","orchid3")) +
  coord_flip() +
  labs(x="gene biotype", y="Normalized Enrichment Score",
       title="Chromosome NES from GSEA") + 
  theme_minimal()

#########
## GSEA of chromosome enrichment using only protein coding genes
####################
TSPvsIP.results.filter <- annotated.TSPvsIP %>% filter(padj < 0.1) # p-value filter
# make a list of genes for the chromosomes and use this instead of pathways
genes.by.chromosome <- list()
unique.chromosomes <- unique(human.geneInfo.biomaRt$chromosome_name) # get list of chromosome names
protein.coding.human.geneInfo.biomaRt <- human.geneInfo.biomaRt %>% filter(gene_biotype == "protein_coding")
for (i in 1:length(unique.chromosomes)) {
  x <- protein.coding.human.geneInfo.biomaRt %>% filter(chromosome_name == unique.chromosomes[i]) %>% dplyr::select(c("ensembl_gene_id"))
  x <- as.character(x[,c(1)])
  genes.by.chromosome[[unique.chromosomes[i]]] <- x
}
# select only protein coding genes
TSPvsIP.results.filter <- TSPvsIP.results.filter %>% filter(ensembl_gene_id %in% protein.coding.human.geneInfo.biomaRt$ensembl_gene_id) # use only protein coding genes
# rank genes
res2 <- TSPvsIP.results.filter %>% 
  dplyr::select(ensembl_gene_id, stat) %>% 
  na.omit() %>% 
  distinct() %>% 
  group_by(ensembl_gene_id) %>% 
  summarize(stat=mean(stat))
res2

ranks <- deframe(res2)
barplot(sort(ranks, decreasing = T)) # make a plot of the ranks
fgseaRes <- fgsea(pathways=genes.by.chromosome, stats=ranks, nperm=10000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy %>% 
  dplyr::select(-leadingEdge, -ES, -nMoreExtreme) %>% 
  arrange(padj) %>% 
  DT::datatable()
#plot
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.01)) +
  scale_fill_manual(values=c("deepskyblue","orchid3")) +
  coord_flip() +
  labs(x="Chromosome", y="Normalized Enrichment Score",
       title="Chromosome NES from GSEA") + 
  theme_minimal()

############################
## Plotting specific marker genes
############################
human.geneInfo.biomaRt <- read.table(file = "path.to.file/human.geneInfo.biomaRt.tab"), header = T, sep = "\t", stringsAsFactors = F) # load biomaRt info from saved table instead of biomaRt; faster
vst <- vst(dds, blind=TRUE) # faster version of rld and almost equivalent
vst.assay <- as.data.frame(assay(vst)) #extract the vst normalized counts
vst.assay$ensembl_gene_id <- gsub("\\..*","", rownames(vst.assay)) #get rid of the version number
count.table <- left_join(vst.assay, human.geneInfo.biomaRt)
count.table$duplicate <- duplicated(count.table$ensembl_gene_id)
duplicate.genes <- count.table %>% filter(duplicate == "TRUE") # figure out which rows were duplicated
duplicated.rows <- count.table %>% filter(ensembl_gene_id %in% duplicate.genes$ensembl_gene_id) # see if they are actually duplicates; they are
count.table <- count.table %>% filter(duplicate == "FALSE") # remove duplicates

### plot all genes together
gaba.genes <- c("ENSG00000078018", "ENSG00000128683", "ENSG00000136750", "ENSG00000132688", "ENSG00000181449", "ENSG00000007372", "ENSG00000070748")
gaba.counts <- vst.assay %>% filter(ensembl_gene_id %in% gaba.genes) 
transpose.gaba.counts <- t(gaba.counts[,c(1:16)])
colnames(transpose.gaba.counts) <- gaba.counts$ensembl_gene_id
meta.gaba.counts <- cbind(transpose.gaba.counts, meta.data)

melt.counts <- melt(meta.gaba.counts, id.vars = c("Cell.Type", "Sample_Type", "paired.samples"), measure.vars = colnames(meta.gaba.counts[1:NCOL(transpose.gaba.counts)]))
melt.counts$sample.combo <- paste(melt.counts$Cell.Type, melt.counts$Sample_Type, sep = "-")
melt.counts$sample.combo <- factor(melt.counts$sample.combo, levels = c("MN-TSP", "MN-IP", "GABA-TSP", "GABA-IP", "NPC-TSP", "NPC-IP"))
melt.counts$variable <- factor(melt.counts$variable, levels = c("ENSG00000078018", "ENSG00000128683", "ENSG00000136750", "ENSG00000132688", "ENSG00000181449", "ENSG00000007372", "ENSG00000070748"))

ggplot(melt.counts, aes(x=variable, y=value)) + 
  #geom_line(aes(group = pairs)) + #connects paired samples
  geom_point(aes(color = Cell.Type, shape=as.factor(Sample_Type)), size = 3, position = position_dodge2(width = 0.45)) +
  theme_minimal() +
  geom_boxplot(aes(color=Cell.Type)) +
  theme(axis.title.x=element_text(hjust = 0.5),
        axis.title.y=element_text(hjust = 0.5),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="right",
        plot.title = element_text(hjust = 0.5)) +
  scale_shape_manual(name = "Sample Type",
                     values = c('IP'=1,
                                'Input'=16),
                     labels = c('IP', 'Input')) 

