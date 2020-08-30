# This script was used to import quantification files into R envrinment via tximport and
# to do differential gene expression analysis (DEG) on them via deseq2.
# This script was used for the C57Bl/6J wildtype vs SwissOF1 wildtype comparison
# Author(s): Dorottya Ralbovszki
# Created:2020.01.20.

rm(list = ls()) #Clear workspace

# importing sample info of samples
sampleinfo <- read.delim("sample_metadata_WT_new.txt")
View(sampleinfo)
sampleinfo

# importing quantification files
dir <- list.files("salmon_data/WT/")
quant_files <- list.files("salmon_data/WT/", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
names(quant_files) <- dir
quant_files

# checking the imported files
library(readr)
quants <- read_tsv(quant_files[1])
head(quants)


# creating transcript database from the gencode mouse genome
library(GenomicFeatures)
txdb <- makeTxDbFromGFF(file="gencode.vM24.annotation.gtf",
                        organism="Mus musculus")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")
head(tx2gene)


# importing and summarizing quantification files into matrix using tximport
library(tximport)
tx2gene <- tx_map
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

names(txi)
head(txi$counts)

# DEG analysis using deseq2
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design <- ~ Strain)

# exporting TPM
tpm <- txi$abudance
write.csv(tpm, file = "tmp_values_WT_new.csv", quote = FALSE)

# cheking in how many samples genes are expressed
is_expressed <- assay(dds) >= 5
head(is_expressed)
sum(is_expressed[1,])
sum(is_expressed[2,])
hist(rowSums(is_expressed),main="Number of samples a gene is expressed in",xlab="Sample Count")

# filtering out genes that had a lower read number than 5 when the read number of all samples (6) was summarized
keep <- rowSums(counts(dds)) >= 5
dds <- dds[keep,]

# visualising count distributions
boxplot(assay(dds))
boxplot(log10(assay(dds)))

# setting the right comparison conditions
dds$Strain = relevel(dds$Strain, "SWISSOF1")
# DEG with the new comparison settings
dds <- DESeq(dds)
# extracting result table from DEG analysis
res_wt <- results(dds, contrast=c("Strain", "C57", "SWISSOF1"))
resultsNames(dds)

# plot counts of smallest p-value
plotCounts(dds, gene=which.min(res_wt$padj), intgroup="Strain")

# remove string after period to get actual ENSEMBL ID
tmp = gsub("\\..*","",row.names(res_wt)) 
row.names(res_wt) = tmp
head(row.names(res_wt))

# order the results by p-values
res_wtOrdered <- res_wt[order(res_wt$pvalue),]
# save separetaly the down-, and upregulated genes
resup_wt <- subset(res_wt, log2FoldChange>0)
resdown_wt <- subset(res_wt, log2FoldChange<0)

# getting the number of significant genes
sum(res_wt$padj < 0.05, na.rm=TRUE)

# summary of analysis
summary(res_wt)
res05 <- results(dds, alpha=0.05)
summary(res05)

# filtering out too high log2fold changes which means that a gene was only expressed in 1 sample/strain
keep_logfold_p <- res_wtOrdered$log2FoldChange <= 10
res_wtOrdered <- res_wtOrdered[keep_logfold_p,]
keep_logfold_n <- res_wtOrdered$log2FoldChange >= -10
res_wtOrdered <- res_wtOrdered[keep_logfold_n,]

# checking filtered result table
res_wtOrdered

# export results into a CSV file
write.csv( as.data.frame(res_wtOrdered), file="results_genotype_effect_new_0406.csv" )

# getting the number of significant genes after filtering
sum(res_wtOrdered$padj < 0.05, na.rm=TRUE)

# annotating result table with gene sybol, entrez ID and gene name
library("AnnotationDbi")
library("org.Mm.eg.db")

annots_symbol <- select(org.Mm.eg.db, keys = rownames(res_wtOrdered), column = "SYMBOL", keytype = "ENSEMBL")
annots_entrez <- select(org.Mm.eg.db, keys = rownames(res_wtOrdered), column = "ENTREZID", keytype = "ENSEMBL")
annots_name <- select(org.Mm.eg.db, keys = rownames(res_wtOrdered), column = "GENENAME", keytype = "ENSEMBL")

# exporting annotated results into a csv file
write.csv( as.data.frame(annots_name), file="annots_name_wt.csv" )
write.csv( as.data.frame(annots_entrez), file="annots_entrez_wt.csv" )
write.csv( as.data.frame(annots_symbol), file="annots_symbol_wt0407.csv" )


# log fold change shrinkage for visualization and ranking 
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Strain_C57Bl.6_vs_SWISSOF1", type="apeglm")
resLFC



# MA-plot showing the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples 
# coloured points showing p-value less than 0.1
# comparint the raw and the log fold change shrinkage
plotMA(res_wt, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# variance stabilizing transformations (VST)
vsd <- vst(dds, blind=FALSE)

# plot PCA of the transformed data
plotPCA(vsd, intgroup=c("Strain"))

# creating sample-to-sample distance plot
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Strain","Genotype")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Strain, vsd$Genotype, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# creating heatmap of the  20 genes with the highest variance across samples
topVarGenes <- head( order( rowVars( assay(vsd) ), 
                            decreasing=TRUE ), 20 )
pheatmap( assay(vsd)[ topVarGenes, ], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#p-value histogram (diagnostic plot for multiple testing)
use <- res_wt$baseMean > metadata(res_wt)$filterThreshold
h <- hist(res_wt$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c("powderblue")
barplot(height = rbind( h$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")

hist(res_wt$pvalue[res_wt$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# Check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd)), col="blue")





# creating significant genes heatmap
# only worked if sample metadata was created this way
sampleTable <- data.frame(Strain = factor(rep(c("C57", "SWISSOF1"), 
                                                each = 3)))

rownames(sampleTable) <- colnames(txi$counts)
sampleTable

# extracting significant genes
expmatrix_DESeq <- DESeq2::rlog(dds, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)

library(gplots)
genes_significant <- expmatrix[c("ENSMUSG00000000308.14",
                                 "ENSMUSG00000000384.15",
                                 "ENSMUSG00000003476.16",
                                 "ENSMUSG00000003929.11",
                                 "ENSMUSG00000004895.9",
                                 "ENSMUSG00000005973.6",
                                 "ENSMUSG00000009734.18",
                                 "ENSMUSG00000012483.4",
                                 "ENSMUSG00000012519.14",
                                 "ENSMUSG00000017929.13",
                                 "ENSMUSG00000017978.18",
                                 "ENSMUSG00000018470.8",
                                 "ENSMUSG00000018698.15",
                                 "ENSMUSG00000019828.13",
                                 "ENSMUSG00000020142.12",
                                 "ENSMUSG00000021070.6",
                                 "ENSMUSG00000021319.7",
                                 "ENSMUSG00000021337.8",
                                 "ENSMUSG00000021880.7",
                                 "ENSMUSG00000022103.10",
                                 "ENSMUSG00000024014.8",
                                 "ENSMUSG00000024140.10",
                                 "ENSMUSG00000024942.17",
                                 "ENSMUSG00000025789.9",
                                 "ENSMUSG00000025870.10",
                                 "ENSMUSG00000026018.12",
                                 "ENSMUSG00000026098.13",
                                 "ENSMUSG00000026113.17",
                                 "ENSMUSG00000026185.8",
                                 "ENSMUSG00000026237.5",
                                 "ENSMUSG00000026516.8",
                                 "ENSMUSG00000026688.5",
                                 "ENSMUSG00000026765.12",
                                 "ENSMUSG00000027224.14",
                                 "ENSMUSG00000027274.16",
                                 "ENSMUSG00000027400.11",
                                 "ENSMUSG00000027792.11",
                                 "ENSMUSG00000028487.18",
                                 "ENSMUSG00000028602.12",
                                 "ENSMUSG00000028635.7",
                                 "ENSMUSG00000028656.14",
                                 "ENSMUSG00000028757.4",
                                 "ENSMUSG00000028901.13",
                                 "ENSMUSG00000029288.11",
                                 "ENSMUSG00000029754.13",
                                 "ENSMUSG00000030532.6",
                                 "ENSMUSG00000031212.3",
                                 "ENSMUSG00000031297.14",
                                 "ENSMUSG00000031391.18",
                                 "ENSMUSG00000032259.8",
                                 "ENSMUSG00000032271.13",
                                 "ENSMUSG00000032643.12",
                                 "ENSMUSG00000032679.12",
                                 "ENSMUSG00000033808.16",
                                 "ENSMUSG00000034652.12",
                                 "ENSMUSG00000034796.14",
                                 "ENSMUSG00000035277.15",
                                 "ENSMUSG00000035726.8",
                                 "ENSMUSG00000035929.11",
                                 "ENSMUSG00000036131.12",
                                 "ENSMUSG00000037025.11",
                                 "ENSMUSG00000037962.7",
                                 "ENSMUSG00000037990.18",
                                 "ENSMUSG00000038257.9",
                                 "ENSMUSG00000039231.18",
                                 "ENSMUSG00000039488.15",
                                 "ENSMUSG00000039579.15",
                                 "ENSMUSG00000039977.16",
                                 "ENSMUSG00000040998.18",
                                 "ENSMUSG00000041911.3",
                                 "ENSMUSG00000042369.8",
                                 "ENSMUSG00000042501.12",
                                 "ENSMUSG00000042770.8",
                                 "ENSMUSG00000042772.15",
                                 "ENSMUSG00000043671.14",
                                 "ENSMUSG00000044068.7",
                                 "ENSMUSG00000044566.15",
                                 "ENSMUSG00000044708.5",
                                 "ENSMUSG00000046500.8",
                                 "ENSMUSG00000047766.15",
                                 "ENSMUSG00000050558.13",
                                 "ENSMUSG00000050711.7",
                                 "ENSMUSG00000051246.3",
                                 "ENSMUSG00000051397.5",
                                 "ENSMUSG00000051747.15",
                                 "ENSMUSG00000052926.16",
                                 "ENSMUSG00000053310.11",
                                 "ENSMUSG00000055202.11",
                                 "ENSMUSG00000055301.8",
                                 "ENSMUSG00000055675.6",
                                 "ENSMUSG00000056596.8",
                                 "ENSMUSG00000058400.13",
                                 "ENSMUSG00000058897.18",
                                 "ENSMUSG00000059040.5",
                                 "ENSMUSG00000061414.8",
                                 "ENSMUSG00000063698.9",
                                 "ENSMUSG00000064329.13",
                                 "ENSMUSG00000064330.9",
                                 "ENSMUSG00000066361.3",
                                 "ENSMUSG00000066438.6",
                                 "ENSMUSG00000068396.9",
                                 "ENSMUSG00000070056.6",
                                 "ENSMUSG00000070880.10",
                                 "ENSMUSG00000071369.11",
                                 "ENSMUSG00000071470.4",
                                 "ENSMUSG00000072437.4",
                                 "ENSMUSG00000072812.4",
                                 "ENSMUSG00000073876.3",
                                 "ENSMUSG00000074269.10",
                                 "ENSMUSG00000074731.3",
                                 "ENSMUSG00000074735.2",
                                 "ENSMUSG00000075705.12",
                                 "ENSMUSG00000078503.9",
                                 "ENSMUSG00000078735.4",
                                 "ENSMUSG00000078954.9",
                                 "ENSMUSG00000079588.3",
                                 "ENSMUSG00000079685.10",
                                 "ENSMUSG00000086600.8",
                                 "ENSMUSG00000087369.1",
                                 "ENSMUSG00000092116.1",
                                 "ENSMUSG00000094686.1",
                                 "ENSMUSG00000095595.2",
                                 "ENSMUSG00000096449.2",
                                 "ENSMUSG00000096995.2",
                                 "ENSMUSG00000097462.7",
                                 "ENSMUSG00000101969.2",
                                 "ENSMUSG00000104178.1",
                                 "ENSMUSG00000105987.4",
                                 "ENSMUSG00000111785.1",
                                 "ENSMUSG00000116819.1")
                               ,]
# ordering genes
genes_significant <- genes_significant[order(rowMeans(genes_significant),
                                             decreasing = TRUE),]
# plotting heatmap
pheatmap::pheatmap(genes_significant, 
                   cluster_rows=FALSE, 
                   show_rownames=FALSE, 
                   show_colnames = TRUE,
                   cluster_cols=TRUE,
                   annotation_col = sampleTable,
                   border_color = NA)

