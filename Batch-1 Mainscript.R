# Set the prefix for each output file name
outputPrefix <- "DEX_H9-neuron_DESeq2-batch1"

sampleFiles <- c("C1-13-output_basename.counts",
                 "C2-14-output_basename.counts",
                 "D1-16-output_basename.counts",
                 "D2-17-output_basename.counts",
                 "D3-18-output_basename.counts")

sampleNames <- c("C1-13","C2-14",
                 "D1-16","D2-17","D3-18")

sampleCondition <- c("6hr_DMSO","6hr_DMSO",
                     "6hr_DEX","6hr_DEX","6hr_DEX")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, 
                          condition = sampleCondition)

treatments = c("6hr_DMSO","6hr_DEX")

library("DESeq2")

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = ".",
                                       design = ~condition)
colData(ddsHTSeq)$condition <- factor(sampleCondition,
                                      levels = treatments)

#guts
dds <- DESeq(ddsHTSeq)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
res <- results(dds)

#convert ensemble_id_version to ensemble_id
library(stringr)

rownames(dds) <- str_replace(rownames(dds),
                             pattern = ".[0-9]+$",
                             replacement = "")

# results table will be generated using results() which will include:
# log2 fold changes, p values and adjusted p values
res <- results(dds, alpha=0.05)

#summarize result of differential gene expression experiment
summary(res)

#returns the names of the coefficient of interest
resultsNames(dds)

library("apeglm")

#shrinks down effect size of differentially expressed genes
resLFC <- lfcShrink(dds, coef="condition_6hr_DEX_vs_6hr_DMSO", type="apeglm")
plotMA(resLFC, ylim=c(-2,2))

#generate a plot of normalized counts for a sample gene. to be modified using ggplot2
topGene <- plotCounts(dds, gene=which.min(res$padj), returnData=TRUE)

# transform raw counts into normalized values
# DESeq2 has two options:  1) rlog transformed and 2) variance stabilization
# variance stabilization is very good for heatmaps, etc.
rld <- rlogTransformation(dds, blind=T)
vsd <- varianceStabilizingTransformation(dds, blind=T)

#this adds 1 to log2 values
ntd <- normTransform(dds)

library("vsn")
library("hexbin")

#generate a plot of log2+1 values
meanSdPlot(assay(ntd))

#both lines generate the plot of the transformed log values
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#independent filter rejection rate plot
plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)

#pvalue histogram plot
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`Null`="khaki", `Alternative`="powderblue")

#barplot of pvalue histogram
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="Frequency", xlab="P Value")
axis(side=1, at=c(0,50))
text(x = c(0, length(h1$counts)), y=0, labels="",
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

#Principal components plot shows additional but rough clustering of samples
library("genefilter")
library("ggplot2")
library("grDevices")

#Script to test plotCounts on a sample gene. I plan to write a multipurpose script
#that I will use to generate multiple plotCounts for a diverse genes of interest
topGene <- plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)

(ggplot(topGene, aes(x=condition, y=count, color=sampleCondition)) 
  + geom_point(position=position_jitter(width=0.12, height=0), size=3))


rv <- rowVars(assay(rld))
select <- order(rv, decreasing=T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))

# set condition
condition <- treatments
scores <- data.frame(pc$x, sampleCondition)

pcaData <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar")) 

pca1 <- (ggplot(pcaData, aes(x = PC1, y = PC2, label = rownames(pcaData), col = condition))
         + geom_point(size =5)
         + geom_text(vjust = 0, nudge_y = .3, col = "Black")
         + theme_bw()
         + xlab(paste0("PC1: ", percentVar[1], "% variance")) 
         + ylab(paste0("PC2: ", percentVar[2], "% variance")) 
         + ggtitle("PCA plot H9 Libraries - C3-15 removed"))

print(pca1)

ggsave("H9acuteDex-PCAplot-no-C3-15.svg", plot=pca1, width =6.5, height=5, unit='in')

#boxplot of cook's distance for all samples 
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)

#heatmaps, clusterring maps e.t.c
sampleDists <- dist(t(assay(rld)))

#heatmap of sample distance to guage how close the groups are to each other
library("RColorBrewer")
library("pheatmap")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition)
colnames(sampleDistMatrix) <- sampleNames
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, 
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#heatmap of normalized count data for all libraries
select <- order(rowMeans(counts(dds, normalized=TRUE)),
                decreasing=TRUE)[1:50]

df <- as.data.frame(colData(dds)["condition"])

pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col = df)

###################################################################
#scripts to extract result tables from the res DESeqResults object#
###################################################################

# plot the shrunken log2 fold changes against the raw changes:
plot(
  x=res$log2FoldChange,
  y=resLFC$log2FoldChange,pch=20,
  cex=.2,
  col=1+(res$padj < 0.05),
  xlab="raw log2 fold change",
  ylab="shrunken log2 fold change"
)
abline(0,1)

#Creates a copy of res which we call Top_Genes_1
Top_Genes_1 = res
View(as.data.frame(Top_Genes_1))

#Creates a copy of the shrunken res DESeqResults named Top_Genes_1LFC
Top_Genes_1LFC = resLFC
View(as.data.frame(Top_Genes_1LFC))

#Renames columns 3 and 4 of Top_Genes_1LFC as specified below
names(Top_Genes_1LFC)[2] <- 'ShrunkenLog2FoldChange'
names(Top_Genes_1LFC)[3] <- 'ShrunkenlfcSE'
names(Top_Genes_1LFC)[5] <- 'mod_padj'
View(as.data.frame(Top_Genes_1LFC))

resdata_1 <- merge(as.data.frame(Top_Genes_1), 
                   as.data.frame(Top_Genes_1LFC), 
                   by = 'row.names', sort = FALSE) 

View(resdata_1)

resdata_1$mod_padj <- NULL
resdata_1$pvalue.y <- NULL
resdata_1$baseMean.y <- NULL

names(resdata_1)[1] <- 'GeneID'
names(resdata_1)[2] <- 'baseMean'
names(resdata_1)[6] <- 'pvalue'

View(resdata_1)

rescount_1 = as.data.frame(counts(dds, normalized=TRUE))

library("dplyr")
rescount_1 <- tibble::rownames_to_column(rescount_1, var = "GeneID")
View(rescount_1)

resdata_1 <- merge(resdata_1, rescount_1, by = 'GeneID', sort = FALSE)
View(resdata_1)

#subset the genes that had padj less than 0.05
resdata_1_padj <- subset(resdata_1, padj<0.05)
resdata_1_padj <- resdata_1_padj[order(resdata_1_padj$ShrunkenLog2FoldChange, decreasing = TRUE),]
View(resdata_1_padj)

#subset the genes that had pvalues less than 0.05
resdata_1_pval <- subset(resdata_1, pvalue<0.05)
resdata_1_pval <- resdata_1_pval[!(is.na(resdata_1_pval$padj) | resdata_1_pval$padj=="NA"), ]
resdata_1_pval <- resdata_1_pval[order(resdata_1_pval$ShrunkenLog2FoldChange, decreasing = TRUE),]
View(resdata_1_pval)

#Extracting hgnc gene-names from GeneIDs, using biomaRt
library("biomaRt")

# select dataset to use
ensembl <- useEnsembl(biomart = "ensembl",
                      dataset = "hsapiens_gene_ensembl",
                      mirror = "useast")

# Extract hgnc gene names from ensembl ids for differentially expressed genes. 
b <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
           filters = "ensembl_gene_id",
           values = resdata_1_padj$GeneID,
           mart = ensembl)

names(b)[1] <- 'GeneID'

#Merge gene names with resdata_1_padj
resdata_1_padj <- merge(b, resdata_1_padj, all = TRUE)
View(resdata_1_padj)

# Extract hgnc gene names for pvalue cutoff genes
a <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
           filters = "ensembl_gene_id",
           values = resdata_1_pval,
           mart = ensembl)

names(a)[1] <- 'GeneID'

#Merge gene names(a) with resdata_1_pval
resdata_1_pval <- merge(a, resdata_1_pval)

# Extract hgnc gene names from ensembl ids for all expressed genes
d <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
           filters = "ensembl_gene_id",
           values = resdata_1$GeneID,
           mart = ensembl)

# Extract hgnc gene names from ensembl ids for all expressed genes
f <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
           filters = "ensembl_gene_id",
           values = resdata_1$GeneID,
           mart = ensembl)

names(d)[1] <- 'GeneID'
#remove duplicate row with erroneous hgnc name "LINC02210-CRHR1" tagged to "ENSG00000278232"
#ENSG00000278232 should only be associated with "CRHR1", not "LINC02210-CRHR1"
#This is a manual removal of the inappropriate row. this solution might not be effective-
#with future updates. removal done on 4/5/2021
d = d[-20026,]

#Merge gene names with resdata_1 to annotate all expressed genes 
resdata_1_annotated <- merge(d, resdata_1, all = TRUE)
resdata_1_annotated <- resdata_1_annotated[order(resdata_1_annotated$ShrunkenLog2FoldChange, decreasing = TRUE),]
View(resdata_1_annotated)

#get gene names and transcript lengths when they exist
c <- getBM(filter = "ensembl_gene_id",
           value = rownames(res),
           attributes = c("ensembl_gene_id","description","transcript_length","hgnc_symbol"),
           mart=ensembl)

#pick only the longest transcript for each gene ID
c <- group_by(c, ensembl_gene_id, hgnc_symbol) %>% 
  summarize(.,description = unique(description),
            transcript_length = max(transcript_length)) 

#remove duplicate row with erroneous hgnc name "LINC02210-CRHR1" tagged to "ENSG00000278232"
#ENSG00000278232 should only be associated with "CRHR1", not "LINC02210-CRHR1"
#This is a manual removal of the inappropriate row. this solution might not be effective-
#with future updates. removal done on 4/5/2021
c = c[-20026,]

#fetch go_terms of interest for all expressed genes. 
go_c <- getBM(filter = "ensembl_gene_id",
              value = rownames(res),
              attributes = c("ensembl_gene_id","description","go_id","name_1006","namespace_1003"),
              mart = ensembl)

resdata_1_go <- resdata_1[resdata_1$GeneID %in% c$ensembl_gene_id,]

# Extract hgnc gene names from ensembl ids for differentially expressed genes. 
e <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol"),
           filters = "ensembl_gene_id",
           values = resdata_1_go$GeneID,
           mart = ensembl)

names(e)[1] <- 'GeneID'

#remove duplicate row with erroneous hgnc name "LINC02210-CRHR1" tagged to "ENSG00000278232"
#ENSG00000278232 should only be associated with "CRHR1", not "LINC02210-CRHR1"
#This is a manual removal of the inappropriate row. this solution might not be effective-
#with future updates. removal done on 4/5/2021
e = e[-20026,]

#Merge gene names with resdata_1_go to annotate all expressed genes 
resdata_1_go_annotated <- merge(e, resdata_1_go, all = TRUE)

#Gene onthology enrichment analysis with goseq
library("goseq")

#0/1 vector for DE or not DE genes
de <- as.numeric(resdata_1_go_annotated$padj < 0.05)
names(de) <- resdata_1_go_annotated$GeneID
de[is.na(de)] = 0

# length of each gene (extracted earlier using biomart)
len <- c[[4]]

# first try to account for transcript length bias by calculating the
# probability of being DE based purely on gene length
pwf <- nullp(DEgenes=de, "hg19", "ensGene", bias.data=len)

# use the Wallenius approximation to calculate enrichment p-values
GO.wall <- goseq(pwf = pwf, gene2cat = go_c[,c(1,3)])

# do FDR correction on p-values using Benjamini-Hochberg, add to output object
GO.wall <- cbind(
  GO.wall,
  padj_overrepresented=p.adjust(GO.wall$over_represented_pvalue, method="BH"),
  padj_underrepresented=p.adjust(GO.wall$under_represented_pvalue, method="BH")
)

# identify ensembl gene IDs annotated with to 2nd from top enriched GO term
g <- go_c$go_id==GO.wall[2,1]
gids <- go_c[g,1]

# inspect DE results for those genes
resdata_1_go[gids,]

#A function to get the genelists of "numDFinCat" in GO.wall report
getGeneLists <- function(pwf, goterms, genome, ids){
  gene2cat <- getgo(rownames(pwf), genome, ids)
  cat2gene <- split(rep(names(gene2cat), sapply(gene2cat, length)),
                    unlist(gene2cat, use.names = FALSE))
  out <- list()
  for(term in goterms){
    tmp <- pwf[cat2gene[[term]],]
    tmp <- rownames(tmp[tmp$DEgenes > 0, ])
    out[[term]] <- tmp
  }
  out
}

#Generates a of genes of interest associated with each goterm
library("org.Hs.eg.db")
goList <- getGeneLists(pwf, GO.wall$category, "hg19", "ensGene")

#Add the list of GeneIDs to the GO.Wall table
GO.wall$EnsemblID <- sapply(GO.wall$category, function(x) paste0(goList[[x]], collapse = ","))

#Plot top 50 terms overrepresented in the data. 
ggData <- GO.wall[1:50,]
ggData$gene_ratio <- ggData$numDEInCat/ggData$numInCat

#sort the plot according to gene ratio
ggData$term <- factor(ggData$term, levels = ggData$term[order(ggData$gene_ratio)])

#plot top 50 overrepresented GO terms from the GOseq analysis -unfiltered data
gg1 <- ggplot(ggData, aes(x=term, y=gene_ratio, size=numDEInCat, fill=over_represented_pvalue))+
  geom_point(shape=21)+
  scale_size(range = c(2,10))+
  scale_fill_continuous(low = 'red', high = 'blue') +
  facet_grid(. ~ ontology)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))+
  coord_flip()

print(gg1)
dev.off()

#Plot top 20 terms overrepresented in the data.
ggData2 <- ggData[1:20,]

#sort the plot according to gene ratio
ggData2$term <- factor(ggData2$term, levels = ggData2$term[order(ggData2$gene_ratio)])

#plot top 20 overrepresented GO terms from the GOseq analysis -unfiltered data
gg2 <- ggplot(ggData2, aes(x=term, y=gene_ratio, size=numDEInCat, 
                           fill=over_represented_pvalue, color = over_represented_pvalue))+
  geom_point(shape=21)+
  scale_color_continuous(low = 'red', high = 'blue')+
  scale_fill_continuous(low = 'red', high = 'blue')+
  scale_size(range = c(2,10))+
  guides(color=FALSE)+
  labs(title = "H9 Batch1",
       x = "GO:term", 
       y = "Gene Ratio", 
       size = "#DE genes in GO:term", 
       fill="p value")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))+
  coord_flip()

print(gg2)

library("svglite")

ggsave("Top20Enriched_Go-terms_H9_batch1.svg", plot=gg2, width = 7.5, height = 7, unit = 'in')

#output ggData2 file as a table
write.table(ggData2, 
            file = paste0(outputPrefix, "_GOseq-Top20_result.txt"), 
            row.names = FALSE,
            sep = '\t')


#set of script to filter the results of GOseq
ggData_filtered <- subset(GO.wall, ggData$numInCat<=600)
ggData_filtered <- subset(ggData_filtered, ggData_filtered$numInCat>=20)
ggData_filtered <- subset(ggData_filtered, ggData_filtered$numDEInCat>1)
ggData_filtered <- subset(ggData_filtered, ggData_filtered$numInCat<=600)
ggData_filtered <- subset(ggData_filtered, ggData_filtered$over_represented_pvalue<=0.05)

#Plot top 50 terms overrepresented in the filtered results
ggData_filtered_top50 <- ggData_filtered[1:50,]
ggData_filtered_top50$gene_ratio <- ggData_filtered_top50$numDEInCat/ggData_filtered_top50$numInCat

#sort the plot according to gene ratio
ggData_filtered_top50$term <- factor(ggData_filtered_top50$term, levels = ggData_filtered_top50$term[order(ggData_filtered_top50$gene_ratio)])

#plot 50 overrepresented GO terms in the filtered result 
ggfiltered1 <- ggplot(ggData_filtered_top50, aes(x=term, y=gene_ratio, size=numDEInCat, fill=over_represented_pvalue))+
  geom_point(shape=21)+
  scale_size(range = c(2,10))+
  scale_fill_continuous(low = 'red', high = 'blue') +
  facet_grid(. ~ ontology)+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))+
  coord_flip()

print(ggfiltered1)




#output result file as a table
write.table(resdata_1_padj, 
            file = paste0(outputPrefix, "_padj_result1.txt"), 
            row.names = FALSE,
            sep = '\t')

#output pvalue result file as a table 
write.table(resdata_1_pval,
            file = paste0(outputPrefix, "_pvalue_result1.txt"),
            row.names = FALSE,
            sep = '\t')

#output file to be used for creating ranked genelist file for GSEA analysis 
write.table(resdata_1_annotated, 
            file = paste0(outputPrefix, "full_ranked_genelist.txt"),
            row.names = FALSE,
            sep = '\t')

#create a datafile of genes ranked by shrunkenlog2foldchange according to GSEA-prerank file creation specification with the addition of a 3rd column
gsas_prerank_file <- data.frame(resdata_1_annotated$GeneID, resdata_1_annotated$ShrunkenLog2FoldChange, resdata_1_annotated$hgnc_symbol)
names(gsas_prerank_file) <- c("GeneID","ShrunkenLog2FoldChange","hgnc_symbol")

write.table(gsas_prerank_file,
            file = paste0(outputPrefix, "preGSAS_rank_file.rnk.txt"),
            row.names = FALSE,
            sep = '\t')

#create a datafile of genes ranked by shrunkenlog2foldchange according to GSEA-prerank file creation specification
gsas_prerank_file2 <- data.frame(resdata_1_annotated$GeneID, resdata_1_annotated$ShrunkenLog2FoldChange)
names(gsas_prerank_file2) <- c("GeneID","ShrunkenLog2FoldChange")

write.table(gsas_prerank_file2,
            file = paste0(outputPrefix, "preGSAS_rank_file2.rnk.txt"),
            row.names = FALSE,
            sep = '\t')

#create a datafile of genes ranked by shrunkenlog2foldchange according to GSEA-prerank file creation specification
gsas_prerank_file3 <- data.frame(resdata_1_annotated$hgnc_symbol, resdata_1_annotated$ShrunkenLog2FoldChange)
names(gsas_prerank_file3) <- c("hgnc_symbol","ShrunkenLog2FoldChange")

write.table(gsas_prerank_file3,
            file = paste0(outputPrefix, "preGSAS_rank_file3.rnk.txt"),
            row.names = FALSE,
            sep = '\t')

#Merge gene names with resdata_1 to annotate all HGNC genes 
resdata_1_filtered_annotated <- merge(d, resdata_1)
resdata_1_filtered_annotated <- resdata_1_filtered_annotated[!(is.na(resdata_1_filtered_annotated$hgnc_symbol) | resdata_1_filtered_annotated$hgnc_symbol==""), ]
resdata_1_filtered_annotated <- resdata_1_filtered_annotated[order(resdata_1_filtered_annotated$ShrunkenLog2FoldChange, decreasing = TRUE),]

#create a datafile of genes ranked by shrunkenlog2foldchange according to GSEA-prerank file creation specification
gsas_prerank_file4 <- data.frame(resdata_1_filtered_annotated$hgnc_symbol, resdata_1_filtered_annotated$ShrunkenLog2FoldChange)
names(gsas_prerank_file4) <- c("hgnc_symbol","ShrunkenLog2FoldChange")

write.table(gsas_prerank_file4,
            file = paste0(outputPrefix, "preGSAS_rank_file4.rnk.txt"),
            row.names = FALSE,
            sep = '\t')

#subsets a list of genes of interest
select_genes_1 <- resdata_1_annotated[resdata_1_annotated$GeneID %in% c("ENSG00000004478","ENSG00000096060",
                                                                        "ENSG00000113580","ENSG00000078018",
                                                                        "ENSG00000160307","ENSG00000096384",
                                                                        "ENSG00000109971","ENSG00000111640"), ]

#specifies the order of gene ids in the new table
select_genes_1$GeneID  <- factor(select_genes_1$GeneID , levels = c("ENSG00000004478","ENSG00000096060",
                                                                    "ENSG00000113580","ENSG00000078018",
                                                                    "ENSG00000160307","ENSG00000096384",
                                                                    "ENSG00000109971","ENSG00000111640"))

#orders the list of genes according to specification 
select_genes_1 <- select_genes_1[order(select_genes_1$GeneID), ]

select_genes_1$log2FoldChange         <- NULL
select_genes_1$lfcSE                  <- NULL
select_genes_1$stat                   <- NULL
select_genes_1$stat                   <- NULL
select_genes_1$pvalue                 <- NULL
select_genes_1$padj                   <- NULL
select_genes_1$ShrunkenLog2FoldChange <- NULL
select_genes_1$ShrunkenlfcSE          <- NULL

write.table(select_genes_1,
            file = paste0(outputPrefix, "select_gene_panel_diagnosis.txt"),
            row.names = FALSE,
            sep = '\t')

############################################################
##Analysis of genes common to both H9 batch and H1 batch 1##
############################################################

#read the file witn GeneID of DEgenes shared in common with H1 Batch1
commonGenes <- read.delim("common_genes.txt")

#subset the data of common genes from the annotated result table
resdata_commonGenes <- resdata_1_annotated[resdata_1_annotated$GeneID %in% commonGenes$GeneID, ]

#add a new column that indicates the direction of lfc
resdata_commonGenes$Direction <- ifelse(resdata_commonGenes$ShrunkenLog2FoldChange>0, "up", "down")

#rearrange result columns to specification
resdata_commonGenes <- resdata_commonGenes[, colnames(resdata_commonGenes)[c(1:10,16,11:15)]]

#write out result table
write.table(resdata_commonGenes, 
            file = paste0(outputPrefix, "DEgenes_shared_with_H1_Batch1.txt"),
            row.names = FALSE,
            sep = '\t')

##########################################################################################################

#reorder the table using shrunkenlogfold2 as factor
resdata_commonGenes$hgnc_symbol <- factor(resdata_commonGenes$hgnc_symbol, levels = resdata_commonGenes$hgnc_symbol[order(resdata_commonGenes$ShrunkenLog2FoldChange)])

mycols_common2 <- rep("gray", nrow(resdata_commonGenes))
mycols_common2[ resdata_commonGenes$Direction == "up" ] <- "red"
mycols_common2[ resdata_commonGenes$Direction == "down" ]  <- "blue" 

ggCommon_H1_H9 <- ggplot(resdata_commonGenes, aes(x=hgnc_symbol, y=ShrunkenLog2FoldChange))+
  geom_bar(stat='identity', fill=mycols_common2)+
  labs(title = "H9 Batch1",
       x = "Genes", 
       y = "Shrunken Log2 Fold-Change")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))+
  coord_flip()

print(ggCommon_H1_H9)

ggsave("H9_batch1_commonGenes2.svg", plot = ggCommon_H1_H9, width = 1.5, height = 5, unit = 'in')

##############################################################################
##Analysis of genes common to all 3 batches, H9_batch1, H1_batch1, H1_batch2##
##############################################################################

#read the file witn GeneID of DEgenes shared in common with H1 Batch1
commonGenes_all <- read.delim("commonDEgenes-H9_batch1-H1_batch1-H1_batch2.txt")

#subset the data of common genes from the annotated result table
resdata_commonGenes_all <- resdata_1_annotated[resdata_1_annotated$GeneID %in% commonGenes_all$GeneID, ]

#add a new column that indicates the direction of lfc
resdata_commonGenes_all$Direction <- ifelse(resdata_commonGenes_all$ShrunkenLog2FoldChange>0, "up", "down")

#rearrange result columns to specification
resdata_commonGenes_all <- resdata_commonGenes_all[, colnames(resdata_commonGenes_all)[c(1:10,16,11:15)]]

#write out result table
write.table(resdata_commonGenes_all, 
            file = paste0(outputPrefix, "DEgenes_shared_among_all_batches.txt"),
            row.names = FALSE,
            sep = '\t')

############################################################################################################
##scripts to make a volcanoe plot of results with common genes highlighted.################################# 
############################################################################################################

#Setup our custom point color vector 
mycols <- rep("gray", nrow(resdata_1_annotated))
mycols[ resdata_1_annotated$pvalue < 0.05 ] <- "green"
mycols[ resdata_1_annotated$padj < 0.05 ]  <- "red" 
mycols[ resdata_1_annotated$GeneID %in% resdata_commonGenes_all$GeneID ] <- "blue"


#make a volcanoe plot
plot( resdata_1_annotated$ShrunkenLog2FoldChange,  -log(resdata_1_annotated$padj), 
      col=mycols,
      xlim = c(-1, 3),
      xlab="Log2(FoldChange)",
      ylab="-Log(Padj)")
text( resdata_commonGenes_all$ShrunkenLog2FoldChange, -log(resdata_commonGenes_all$padj), resdata_commonGenes_all$hgnc_symbol, 
      cex=0.6,
      pos=4,
      col="black")


###############################################################################################################################################################################################

#sort the plot according to gene ratio
resdata_commonGenes_all$hgnc_symbol <- factor(resdata_commonGenes_all$hgnc_symbol, levels = resdata_commonGenes_all$hgnc_symbol[order(resdata_commonGenes_all$ShrunkenLog2FoldChange)])

mycols_common <- rep("gray", nrow(resdata_commonGenes_all))
mycols_common[ resdata_commonGenes_all$Direction == "up" ] <- "red"
mycols_common[ resdata_commonGenes_all$Direction == "down" ]  <- "blue" 

#plot fold change and direction of genes common to all 3 batches 
ggCommon_all <- ggplot(resdata_commonGenes_all, aes(x=hgnc_symbol, y=ShrunkenLog2FoldChange))+
  geom_bar(stat='identity', fill=mycols_common)+
  labs(title = "H9 Batch1",
       x = "Genes", 
       y = "Shrunken Log2 Fold-Change")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=90))+
  coord_flip()

print(ggCommon_all)

ggsave("H9_batch1_commonGenes.svg", plot = ggCommon_all, width = 3, height = 11, unit = 'in')

########################################################################################
#Script to make Venn Diagram of all gene lists from H9 batch1, H1 batch1, and H1 batch2#
########################################################################################

#read in the DE genelists for H1 batch1, and H1 batch2
H1_batch1 <- read.delim("Dex_H1-neuron_DESeq2_batch1_padj_result1.txt")
H1_batch2 <- read.delim("Dex_H1-neuron_DESeq2_batch2_padj_result1.txt")

#create a new table of DE gene liist named "H9_batch1" for Venn Diagram construction. 
H9_batch1 <- resdata_1_padj

library(ggVennDiagram)

#create a list of all GeneIDs from each gene list
res_Venn <- list(H9_Batch1=sample(H9_batch1$GeneID), H1_Batch1=sample(H1_batch1$GeneID), H1_Batch2=sample(H1_batch2$GeneID))

#draw a Venn diagram showing common genes between all 3 genelists
resVenn_plot <- ggVennDiagram(res_Venn)

ggsave("resVenn_plot.svg", plot = resVenn_plot, width = 7, height = 7, unit = 'in')

library(ggvenn)

#draw the ggvenn version of the Venn diagram above
ggvenn(res_Venn)



library(fgsea)
library(data.table)
library(ggplot2)
library(tidyverse)

pathways_hallmark <- gmtPathways("msigdb.v7.4.symbols.gmt")

#turns gene rank file into a named vector
ranks <- deframe(gsas_prerank_file4)

#does the fgsea analysis 
fgseaRes <- fgsea(pathways=pathways_hallmark, stats=ranks, eps=0)

#Tidy up and sort the table using NES column
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))

topPathwaysUp <- fgseaRes[ES > 0][head(order(padj), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(padj), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

plotGseaTable(pathways_hallmark[topPathways], ranks, fgseaRes, 
              gseaParam=0.5)

collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      pathways_hallmark, ranks)

mainPathways <- fgseaRes[fgseaRes$pathway %in% collapsedPathways$mainPathways]

mainPathways <- mainPathways[order(mainPathways$NES, decreasing = TRUE),] 


#script to construct a table of top most enriched pathwya for upregulated or downregulated transcript 
topMainPathwaysUp <- mainPathways[ES > 0][tail(order(NES), n=10), pathway]
topMainPathwaysDown <- mainPathways[ES < 0][head(order(NES), n=10), pathway]
topMainPathways <- c(rev(topMainPathwaysUp), rev(topMainPathwaysDown))

#script to plot our constructed dataframe of top enriched pathways 
plotGseaTable(pathways_hallmark[topMainPathways], ranks, fgseaRes, 
              gseaParam=0.5)

#write out the table of collapsed main enriched pathways
fwrite(mainPathways, 
       file = paste0(outputPrefix, "_GSEA_mainPathways_Enrichment_result.txt"), 
       sep="\t", 
       sep2=c("", " ", ""))

#write out the table of full list of enriched pathways
fwrite(fgseaRes, 
       file = paste0(outputPrefix, "_GSEA_Complete_Enrichment_result.txt"), 
       sep="\t", 
       sep2=c("", " ", ""))

#######################################################################################################
##code to subset curated genesets from fgsea results###################################################
#######################################################################################################

#put into Hallmark_sets all rows where pathway starts with 'HALLMARK'
Hallmark_sets <- subset(fgseaRes, grepl('HALLMARK', pathway))

#sort Hallmark_sets by the normalized enrichment score
Hallmark_sets <- Hallmark_sets[order(Hallmark_sets$NES, decreasing = TRUE), ]

ggplot(Hallmark_sets, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

######################################################################################################

#put into 'cortico' filter Gene_sets all rows where pathway starts with 'cortico'
Cortic_filter_sets <- subset(fgseaRes, grepl('CORTICO', pathway))

#bind the dataframe Cortic_filter_sets to 'Dex' filter result
Cortic_filter_sets <- rbind(Cortic_filter_sets, subset(fgseaRes, grepl('DEX', pathway)))


#sort 'cortico' by the normalized enrichment score
Cortic_filter_sets <- Cortic_filter_sets[order(Cortic_filter_sets$NES, decreasing = TRUE), ]

ggplot(Cortic_filter_sets, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Filter of genesets with 'cortico' in their name-pathways NES from GSEA") + 
  theme_minimal()

######################################################################################################

#put into 'STEROID' filter Gene_sets all rows where pathway starts with 'STEROID'
Steroid_filter_sets <- subset(fgseaRes, grepl('STEROID', pathway))

#sort 'cortico' by the normalized enrichment score
Steroid_filter_sets <- Steroid_filter_sets[order(Steroid_filter_sets$NES, decreasing = TRUE), ]

ggplot(Steroid_filter_sets, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Filter of genesets with 'Steroid' in their name-pathways NES from GSEA") + 
  theme_minimal()

 #####################################################################################################

#script to plot the GOBP Cellular response to corticosteroid stimulus
plotEnrichment(pathways_hallmark[["GOBP_CELLULAR_RESPONSE_TO_CORTICOSTEROID_STIMULUS"]],
               ranks) + labs(title="GOBP CELLULAR RESPONSE TO CORTICOSTEROID STIMULUS")

#script to plot the fan embryonic CTX ex4 excitatory neuron
plotEnrichment(pathways_hallmark[["FAN_EMBRYONIC_CTX_EX_4_EXCITATORY_NEURON"]],
               ranks) + labs(title="FAN EMBRYONIC CTX EX4 EXCITATORY NEURON")

######################################################################################################

#script to plot the GOCC set for synapes
plotEnrichment(pathways_hallmark[["GOCC_SYNAPSE"]],
               ranks) + labs(title="GOCC SYNAPSE")

#script to plot the GOCC set for anchoring junction
plotEnrichment(pathways_hallmark[["GOCC_ANCHORING_JUNCTION"]],
               ranks) + labs(title="GOCC ANCHORING JUNCTION")

#script to plot the GOCC set for endoplasmic reticulum
plotEnrichment(pathways_hallmark[["GOCC_ENDOPLASMIC_RETICULUM"]],
               ranks) + labs(title="GOCC ENDOPLASMIC RETICULUM")

######################################################################################################

#create a ranked file similar to rank4 file, but wth ensembl gene id, instead of hgnc symbols
gsas_prerank_file5 <- data.frame(resdata_1_filtered_annotated$GeneID, resdata_1_filtered_annotated$ShrunkenLog2FoldChange)
names(gsas_prerank_file5) <- c("ensembl_gene_id","ShrunkenLog2FoldChange")

# Extract hgnc gene names from ensembl ids for all expressed genes
d2 <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),
            filters = "ensembl_gene_id",
            values = gsas_prerank_file5$ensembl_gene_id,
            mart = ensembl)

gsas_prerank_file5 <- merge(d2, gsas_prerank_file5)
gsas_prerank_file5 <- gsas_prerank_file5[!(is.na(gsas_prerank_file5$entrezgene_id) | gsas_prerank_file5$entrezgene_id==""), ]
gsas_prerank_file5 <- gsas_prerank_file5[order(gsas_prerank_file5$ShrunkenLog2FoldChange, decreasing = TRUE),]


library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

geneList <- gsas_prerank_file5
geneList$ensembl_gene_id <- NULL
geneList$hgnc_symbol <- NULL

data(geneList)

gene <- names(geneList)[abs(geneList) > 2]

c5 <- read.gmt("msigdb.v7.4.entrez.gmt")

em2 <- GSEA(geneList, TERM2GENE = c5, eps = 0, verbose=FALSE)
