# This script was used to import quantification files into R envrinment via tximport and
# to do differential gene expression analysis (DEG) on them via deseq2.
# This script was used for the C57Bl/6J En1+/- vs SwissOF1 En1+/- comparison
# Author(s): Dorottya Ralbovszki
# Created:2020.01.20.

rm(list = ls()) #Clear workspace

# importing sample info of samples
sampleinfo <- read.delim("sample_metadata_mut_new.txt")
View(sampleinfo)
sampleinfo

# importing quantification files
dir <- list.files("salmon_data/MUT/")
quant_files <- list.files("salmon_data/MUT/", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
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
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")

names(txi)
head(txi$counts)

# DEG analysis using deseq2
library(DESeq2)
dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design <- ~ Strain)

# exporting TPM
tpm <- txi$abudance
write.csv(tpm, file = "tmp_values_MUT_new.csv", quote = FALSE)

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
res_mut <- results(dds, contrast=c("Strain", "C57", "SWISSOF1"))
resultsNames(dds)

# plot counts of smallest p value
plotCounts(dds, gene=which.min(res_mut$padj), intgroup="Strain")

# remove string after period to get actual ENSEMBL ID
tmp = gsub("\\..*","",row.names(res_mut)) 
row.names(res_mut) = tmp
head(row.names(res_mut))


# order the results by p values
res_mutOrdered <- res_mut[order(res_mut$pvalue),]
# save separetaly the down-, and upregulated genes
resup_mut <- subset(res_mut, log2FoldChange>0)
resdown_mut <- subset(res_mut, log2FoldChange<0)

# getting the number of significant genes
sum(res_mut$padj < 0.05, na.rm=TRUE)

# summary of analysis
summary(res_mut)
res05 <- results(dds, alpha=0.05)
summary(res05)

# filtering out too high log2fold changes which means that a gene was only expressed in 1 sample/strain
keep_logfold_p <- res_mutOrdered$log2FoldChange <= 10
res_mutOrdered <- res_mutOrdered[keep_logfold_p,]
keep_logfold_n <- res_mutOrdered$log2FoldChange >= -10
res_mutOrdered <- res_mutOrdered[keep_logfold_n,]

# checking filtered result table
res_mutOrdered

# export results into a CSV file
write.csv( as.data.frame(res_mutOrdered), file="results_genotype_effect_mut_new_0406.csv" )

# getting the number of significant genes after filtering
sum(res_mutOrdered$padj < 0.05, na.rm=TRUE)

# annotating result table with gene sybol, entrez ID and gene name
library("AnnotationDbi")
library("org.Mm.eg.db")

annots_symbol <- select(org.Mm.eg.db, keys = rownames(res_mutOrdered), column = "SYMBOL", keytype = "ENSEMBL")
annots_entrez <- select(org.Mm.eg.db, keys = rownames(res_mutOrdered), column = "ENTREZID", keytype = "ENSEMBL")
annots_name <- select(org.Mm.eg.db, keys = rownames(res_mutOrdered), column = "GENENAME", keytype = "ENSEMBL")

# exporting annotated results into a csv file
write.csv( as.data.frame(annots_name), file="annots_name_mut0407.csv" )
write.csv( as.data.frame(annots_entrez), file="annots_entrez_mut0407.csv" )
write.csv( as.data.frame(annots_symbol), file="annots_symbol_mut0407.csv" )


# log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Strain_C57_vs_SWISSOF1", type="apeglm")



# MA-plot showing the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples 
# coloured points showing p-value less than 0.1
# comparint the raw and the log fold change shrinkage
plotMA(res_mut, ylim=c(-2,2))
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

# creating heatmap of the  20 genes with the highest variance across samples
topVarGenes <- head( order( rowVars( assay(vsd) ), 
                            decreasing=TRUE ), 20 )
pheatmap( assay(vsd)[ topVarGenes, ], scale="row",
          trace="none", dendrogram="column",
          col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255))

#p-value histogram (diagnostic plot for multiple testing)
use <- res_mut$baseMean > metadata(res_mut)$filterThreshold
h <- hist(res_mut$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c("powderblue")
barplot(height = rbind( h$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")

hist(res_mut$pvalue[res_mut$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# Let's add a blue horizontal line that corresponds to the median logCPM
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
genes_significant <- expmatrix[c("ENSMUSG00000000805.18",
                                 "ENSMUSG00000001300.16",
                                 "ENSMUSG00000005447.12",
                                 "ENSMUSG00000013523.13",
                                 "ENSMUSG00000019865.9",
                                 "ENSMUSG00000019865.9",
                                 "ENSMUSG00000020173.17",
                                 "ENSMUSG00000020253.15",
                                 "ENSMUSG00000020524.16",
                                 "ENSMUSG00000020723.3",
                                 "ENSMUSG00000021193.10",
                                 "ENSMUSG00000021567.15",
                                 "ENSMUSG00000021867.16",
                                 "ENSMUSG00000021969.8",
                                 "ENSMUSG00000022048.8",
                                 "ENSMUSG00000023439.11",
                                 "ENSMUSG00000024565.10",
                                 "ENSMUSG00000024713.16",
                                 "ENSMUSG00000024734.8",
                                 "ENSMUSG00000025804.5",
                                 "ENSMUSG00000026638.15",
                                 "ENSMUSG00000027217.13",
                                 "ENSMUSG00000027400.11",
                                 "ENSMUSG00000027570.15",
                                 "ENSMUSG00000028023.16",
                                 "ENSMUSG00000028370.7",
                                 "ENSMUSG00000028558.14",
                                 "ENSMUSG00000028841.14",
                                 "ENSMUSG00000028883.17",
                                 "ENSMUSG00000029005.4",
                                 "ENSMUSG00000029193.7",
                                 "ENSMUSG00000029428.13",
                                 "ENSMUSG00000030123.15",
                                 "ENSMUSG00000030235.17",
                                 "ENSMUSG00000030761.16",
                                 "ENSMUSG00000030792.8",
                                 "ENSMUSG00000031558.15",
                                 "ENSMUSG00000032076.20",
                                 "ENSMUSG00000032854.12",
                                 "ENSMUSG00000033597.9",
                                 "ENSMUSG00000035513.19",
                                 "ENSMUSG00000037362.8",
                                 "ENSMUSG00000038007.14",
                                 "ENSMUSG00000038173.15",
                                 "ENSMUSG00000038738.15",
                                 "ENSMUSG00000039106.6",
                                 "ENSMUSG00000039126.10",
                                 "ENSMUSG00000039735.16",
                                 "ENSMUSG00000040543.16",
                                 "ENSMUSG00000041607.17",
                                 "ENSMUSG00000041975.17",
                                 "ENSMUSG00000042379.8",
                                 "ENSMUSG00000045573.9",
                                 "ENSMUSG00000046610.15",
                                 "ENSMUSG00000049336.16",
                                 "ENSMUSG00000049744.15",
                                 "ENSMUSG00000050148.9",
                                 "ENSMUSG00000054457.5",
                                 "ENSMUSG00000055301.8",
                                 "ENSMUSG00000056158.14",
                                 "ENSMUSG00000056608.14",
                                 "ENSMUSG00000057123.14",
                                 "ENSMUSG00000058174.7",
                                 "ENSMUSG00000059325.14",
                                 "ENSMUSG00000059839.9",
                                 "ENSMUSG00000060550.16",
                                 "ENSMUSG00000060550.16",
                                 "ENSMUSG00000060860.8",
                                 "ENSMUSG00000061762.12",
                                 "ENSMUSG00000063171.4",
                                 "ENSMUSG00000066705.7",
                                 "ENSMUSG00000069072.9",
                                 "ENSMUSG00000070583.1",
                                 "ENSMUSG00000070605.4",
                                 "ENSMUSG00000075296.5",
                                 "ENSMUSG00000076498.2",
                                 "ENSMUSG00000078591.1",
                                 "ENSMUSG00000082361.6",
                                 "ENSMUSG00000091705.8",
                                 "ENSMUSG00000093483.1",
                                 "ENSMUSG00000097039.8",
                                 "ENSMUSG00000097431.2",
                                 "ENSMUSG00000115529.1")
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
