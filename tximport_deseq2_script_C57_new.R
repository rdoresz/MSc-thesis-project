# This script was used to import quantification files into R envrinment via tximport and
# to do differential gene expression analysis (DEG) on them via deseq2.
# This script was used for the C57Bl/6J wildtype vs C57Bl/6J En1+/- comparison
# Author(s): Dorottya Ralbovszki
# Created:2020.01.20.

rm(list = ls()) #Clear workspace

# importing sample info of samples
sampleinfo <- read.delim("sample_metadata_C57_new_new.txt")
View(sampleinfo)
sampleinfo

# importing quantification files
dir <- list.files("salmon_data/C57/")
quant_files <- list.files("salmon_data/C57", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
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
dds <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design <- ~ Genotype)

# exporting TPM
tpm <- txi$abudance
# write.csv(tpm, file = "tmp_values_C57_new.csv", quote = FALSE)

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
dds$Genotype = relevel(dds$Genotype,"wildtype")
# DEG with the new comparison settings
dds <- DESeq(dds)
# extracting result table from DEG analysis
res_c57 <- results(dds, contrast=c("Genotype", "mutant", "wildtype"))
resultsNames(dds)

# plot counts of smallest p value
plotCounts(dds, gene=which.min(res_c57$padj), intgroup="Genotype")

# remove string after period to get actual ENSEMBL ID
tmp = gsub("\\..*","",row.names(res_c57)) 
row.names(res_c57) = tmp
head(row.names(res_c57))

# order the results by p-values
res_c57Ordered <- res_c57[order(res_c57$pvalue),]
# save separetaly the down-, and upregulated genes
resup_c57 <- subset(res_c57, log2FoldChange>0)
resdown_c57 <- subset(res_c57, log2FoldChange<0)

# getting the number of significant genes
sum(res_c57$padj < 0.05, na.rm=TRUE)

# summary of analysis
summary(res_c57)
res05 <- results(dds, alpha=0.05)
summary(res05)

# filtering out too high log2fold changes which means that a gene was only expressed in 1 sample/strain
keep_logfold_p <- res_c57Ordered$log2FoldChange <= 10
res_c57Ordered <- res_c57Ordered[keep_logfold_p,]
keep_logfold_n <- res_c57Ordered$log2FoldChange >= -10
res_c57Ordered <- res_c57Ordered[keep_logfold_n,]

# export results into a CSV file
write.csv( as.data.frame(res_c57Ordered), file="results_allele_effect_c57_new_0406.csv" )

# getting the number of significant genes after filtering
sum(res_c57Ordered$padj < 0.05, na.rm=TRUE)

# annotating result table with gene sybol, entrez ID and gene name
library("AnnotationDbi")
library("org.Mm.eg.db")

annots_symbol <- select(org.Mm.eg.db, keys = rownames(res_c57Ordered), column = "SYMBOL", keytype = "ENSEMBL")
annots_entrez <- select(org.Mm.eg.db, keys = rownames(res_c57Ordered), column = "ENTREZID", keytype = "ENSEMBL")
annots_name <- select(org.Mm.eg.db, keys = rownames(res_c57Ordered), column = "GENENAME", keytype = "ENSEMBL")

# exporting annotated results into a csv file
write.csv( as.data.frame(annots_name), file="annots_name_c57.csv" )
write.csv( as.data.frame(annots_entrez), file="annots_entrez_c57.csv" )
write.csv( as.data.frame(annots_symbol), file="annots_symbol_c57_0407.csv" )


# log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Genotype_mutant_vs_wildtype", type="apeglm")
resLFC



# MA-plot showing the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples 
# coloured points showing p-value less than 0.1
# comparint the raw and the log fold change shrinkage
plotMA(res_c57, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# variance stabilizing transformations (VST)
vsd <- vst(dds, blind=FALSE)

# plot PCA of the transformed data
plotPCA(vsd, intgroup=c("Genotype"))

# creating sample-to-sample distance plot
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Genotype")])

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=sampleTable_c57)

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
use <- res_c57$baseMean > metadata(res_c57)$filterThreshold
h <- hist(res_c57$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c("powderblue")
barplot(height = rbind( h$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")

hist(res_c57$pvalue[res_c57$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(assay(vsd)), col="blue")




# creating significant genes heatmap
# only worked if sample metadata was created this way
sampleTable <- data.frame(Genotype = factor(rep(c("WT", "MUT"), 
                                                each = 3)))

rownames(sampleTable) <- colnames(txi$counts)
sampleTable

# extracting significant genes
expmatrix_DESeq <- DESeq2::rlog(dds, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)

library(gplots)

genes_significant <- expmatrix[c("ENSMUSG00000001240.13",
                                 "ENSMUSG00000001930.17",
                                 "ENSMUSG00000001946.14",
                                 "ENSMUSG00000004221.16",
                                 "ENSMUSG00000004936.8",
                                 "ENSMUSG00000005299.6",
                                 "ENSMUSG00000016427.7",
                                 "ENSMUSG00000019874.11",
                                 "ENSMUSG00000020218.11",
                                 "ENSMUSG00000020950.10",
                                 "ENSMUSG00000021587.5",
                                 "ENSMUSG00000021647.7",
                                 "ENSMUSG00000021872.8",
                                 "ENSMUSG00000022090.10",
                                 "ENSMUSG00000022667.18",
                                 "ENSMUSG00000024681.11",
                                 "ENSMUSG00000025362.6",
                                 "ENSMUSG00000026048.16",
                                 "ENSMUSG00000026413.12",
                                 "ENSMUSG00000027678.17",
                                 "ENSMUSG00000028234.6",
                                 "ENSMUSG00000029917.15",
                                 "ENSMUSG00000030279.15",
                                 "ENSMUSG00000030413.7",
                                 "ENSMUSG00000031997.9",
                                 "ENSMUSG00000032368.14",
                                 "ENSMUSG00000033960.6",
                                 "ENSMUSG00000034055.16",
                                 "ENSMUSG00000034243.17",
                                 "ENSMUSG00000034723.11",
                                 "ENSMUSG00000034892.8",
                                 "ENSMUSG00000035329.7",
                                 "ENSMUSG00000037025.11",
                                 "ENSMUSG00000037400.17",
                                 "ENSMUSG00000037526.7",
                                 "ENSMUSG00000037784.14",
                                 "ENSMUSG00000038291.16",
                                 "ENSMUSG00000038600.12",
                                 "ENSMUSG00000038738.15",
                                 "ENSMUSG00000039672.12",
                                 "ENSMUSG00000039706.11",
                                 "ENSMUSG00000039714.9",
                                 "ENSMUSG00000041559.7",
                                 "ENSMUSG00000041773.8",
                                 "ENSMUSG00000046480.6",
                                 "ENSMUSG00000049928.15",
                                 "ENSMUSG00000052305.6",
                                 "ENSMUSG00000052305.6",
                                 "ENSMUSG00000053930.13",
                                 "ENSMUSG00000054409.5",
                                 "ENSMUSG00000056306.5",
                                 "ENSMUSG00000058443.5",
                                 "ENSMUSG00000060063.9",
                                 "ENSMUSG00000063260.2",
                                 "ENSMUSG00000067288.13",
                                 "ENSMUSG00000067870.5",
                                 "ENSMUSG00000068117.10",
                                 "ENSMUSG00000068523.12",
                                 "ENSMUSG00000069917.7",
                                 "ENSMUSG00000069919.7",
                                 "ENSMUSG00000073940.3",
                                 "ENSMUSG00000073982.11",
                                 "ENSMUSG00000074695.3",
                                 "ENSMUSG00000075330.4",
                                 "ENSMUSG00000078735.4",
                                 "ENSMUSG00000079018.10",
                                 "ENSMUSG00000079641.3",
                                 "ENSMUSG00000089679.1",
                                 "ENSMUSG00000097622.2",
                                 "ENSMUSG00000100627.6",
                                 "ENSMUSG00000102422.1",
                                 "ENSMUSG00000102644.5",
                                 "ENSMUSG00000113771.1")
                            ,]
rownames(genes_significant)<- c("Ramp2",
                                "Vwf",
                                "Esam",
                                "Ikbkg",
                                "Map2k1",
                                "Letm1",
                                "Ndufa1",
                                "Fabp7",
                                "Wif1",
                                "Foxg1",
                                "Pcsk1",
                                "Cartpt",
                                "Rnase10",
                                "Pdlim2",
                                "Cd200r1",
                                "Ms4a3",
                                "Rps26",
                                "Ercc5",
                                "Pkp1",
                                "Ncoa3",
                                "Rps20",
                                "Qrfprl",
                                "C2cd5",
                                "Pglyrp1",
                                "Trpc6",
                                "Zic1",
                                "Jcad",
                                "Phka1",
                                "Golgb1",
                                "Tmx4",
                                "Rps29",
                                "Fbxo33",
                                "Foxa2",
                                "Atp11b",
                                "Atg14",
                                "Dzip1l",
                                "Snx25",
                                "Atp6v0a4",
                                "Shank1",
                                "Kcne2",
                                "Ldb2",
                                "Cplx3",
                                "Fmod",
                                "Enc1",
                                "Scn4b",
                                "Glp2r",
                                "Hbb-b1",
                                "Hbb-bs",
                                "Shisa6",
                                "Tmem74",
                                "Sertm1",
                                "Rpl10-ps3",
                                "Alox5ap",
                                "Syt10",
                                "Rps28",
                                "Rpl31-ps8",
                                "Mei1",
                                "Gng5",
                                "Hba-a2",
                                "Hba-a1",
                                "Hbb-bt",
                                "Rhog",
                                "Il22",
                                "A930003A15Rik",
                                "Il11ra2",
                                "Ly6c1",
                                "Rpl39",
                                "Gm16299",
                                "A330033J07Rik",
                                "A830008E24Rik",
                                "Iqschfp",
                                "Thap6",
                                "Gm36529")
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

# exporting dds table into csv file
write.csv( as.data.frame(counts(dds)), file="dds_c570507.csv" )