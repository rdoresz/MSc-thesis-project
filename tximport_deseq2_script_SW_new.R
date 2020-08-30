# This script was used to import quantification files into R envrinment via tximport and
# to do differential gene expression analysis (DEG) on them via deseq2.
# This script was used for the C57Bl/6J wildtype vs C57Bl/6J En1+/- comparison
# Author(s): Dorottya Ralbovszki
# Created:2020.01.20.

rm(list = ls()) #Clear workspace

# importing sample info of samples
sampleinfo <- read.delim("sample_metadata_SW_new.txt")
View(sampleinfo)
sampleinfo

# importing quantification files
dir <- list.files("salmon_data/SW/")
quant_files <- list.files("salmon_data/SW/", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
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
write.csv(tpm, file = "tmp_values_SW_new.csv", quote = FALSE)

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
res_sw <- results(dds, contrast=c("Genotype", "mutant", "wildtype"))
resultsNames(dds)

# plot counts of smallest p value
plotCounts(dds, gene=which.min(res_sw$padj), intgroup="Genotype")

# remove string after period to get actual ENSEMBL ID
tmp = gsub("\\..*","",row.names(res_sw)) 
row.names(res_sw) = tmp
head(row.names(res_sw))


# order the results by p values
res_swOrdered <- res_sw[order(res_sw$pvalue),]
# save separetaly the down-, and upregulated genes
resup_sw <- subset(res_sw, log2FoldChange>0)
resdown_sw <- subset(res_sw, log2FoldChange<0)

# getting the number of significant genes
sum(res_sw$padj < 0.05, na.rm=TRUE)

# summary of analysis
summary(res_sw)
res05 <- results(dds, alpha=0.05)
summary(res05)

# filtering out too high log2fold changes which means that a gene was only expressed in 1 sample/strain
keep_logfold_p <- res_swOrdered$log2FoldChange <= 10
res_swOrdered <- res_swOrdered[keep_logfold_p,]
keep_logfold_n <- res_swOrdered$log2FoldChange >= -10
res_swOrdered <- res_swOrdered[keep_logfold_n,]

# checking filtered result table
res_swOrdered


# export results into a CSV file
write.csv( as.data.frame(res_swOrdered), file="results_allele_effect_sw_new_0406.csv" )

# getting the number of significant genes after filtering
sum(res_swOrdered$padj < 0.05, na.rm=TRUE)

# annotating result table with gene sybol, entrez ID and gene name
library("AnnotationDbi")
library("org.Mm.eg.db")

annots_symbol <- select(org.Mm.eg.db, keys = rownames(res_swOrdered), column = "SYMBOL", keytype = "ENSEMBL")
annots_entrez <- select(org.Mm.eg.db, keys = rownames(res_swOrdered), column = "ENTREZID", keytype = "ENSEMBL")
annots_name <- select(org.Mm.eg.db, keys = rownames(res_swOrdered), column = "GENENAME", keytype = "ENSEMBL")

# exporting annotated results into a csv file
write.csv( as.data.frame(annots_name), file="annots_name_sw_0407.csv" )
write.csv( as.data.frame(annots_entrez), file="annots_entrez_sw.csv" )
write.csv( as.data.frame(annots_symbol), file="annots_symbol_sw_0407.csv" )


# log fold change shrinkage for visualization and ranking 
resultsNames(dds)

resLFC



# MA-plot showing the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples 
# coloured points showing p-value less than 0.1
# comparint the raw and the log fold change shrinkage
plotMA(res_sw, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# variance stabilizing transformations (VST)
vsd <- vst(dds, blind=FALSE)

# plot PCA of the transformed data
plotPCA(vsd, intgroup=c("Genotype"))

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
use <- res_sw$baseMean > metadata(res_sw)$filterThreshold
h <- hist(res_sw$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c("powderblue")
barplot(height = rbind( h$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")

hist(res_sw$pvalue[res_sw$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")


# Check distributions of samples using boxplots
boxplot(assay(vsd), xlab="", ylab="Log2 counts per million",las=2,main="Normalised Distributions")
# Let's add a blue horizontal line that corresponds to the median logCPM
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
genes_significant <- expmatrix[c("ENSMUSG00000028023.16",
                                 "ENSMUSG00000028298.10",
                                 "ENSMUSG00000025329.3",
                                 "ENSMUSG00000068697.7",
                                 "ENSMUSG00000074577.9",
                                 "ENSMUSG00000026969.3",
                                 "ENSMUSG00000026638.15",
                                 "ENSMUSG00000100241.1",
                                 "ENSMUSG00000105960.1",
                                 "ENSMUSG00000101028.1",
                                 "ENSMUSG00000035513.19",
                                 "ENSMUSG00000037362.8",
                                 "ENSMUSG00000055775.16",
                                 "ENSMUSG00000042379.8",
                                 "ENSMUSG00000031738.14",
                                 "ENSMUSG00000007034.15",
                                 "ENSMUSG00000047182.6",
                                 "ENSMUSG00000028883.17",
                                 "ENSMUSG00000019865.9",
                                 "ENSMUSG00000001027.7",
                                 "ENSMUSG00000000120.6",
                                 "ENSMUSG00000030123.15",
                                 "ENSMUSG00000025013.15",
                                 "ENSMUSG00000024734.8",
                                 "ENSMUSG00000086017.1",
                                 "ENSMUSG00000020908.14",
                                 "ENSMUSG00000072812.4",
                                 "ENSMUSG00000005892.4",
                                 "ENSMUSG00000027270.14",
                                 "ENSMUSG00000049630.6",
                                 "ENSMUSG00000004885.5",
                                 "ENSMUSG00000047586.4",
                                 "ENSMUSG00000019890.4",
                                 "ENSMUSG00000060550.16",
                                 "ENSMUSG00000060550.16",
                                 "ENSMUSG00000037600.16",
                                 "ENSMUSG00000071489.1",
                                 "ENSMUSG00000024526.9",
                                 "ENSMUSG00000061762.12",
                                 "ENSMUSG00000046410.10",
                                 "ENSMUSG00000024713.16",
                                 "ENSMUSG00000048251.15",
                                 "ENSMUSG00000047810.9",
                                 "ENSMUSG00000079499.9",
                                 "ENSMUSG00000078952.9",
                                 "ENSMUSG00000031558.15",
                                 "ENSMUSG00000021732.14",
                                 "ENSMUSG00000005672.12",
                                 "ENSMUSG00000056418.3",
                                 "ENSMUSG00000040118.15",
                                 "ENSMUSG00000031767.13",
                                 "ENSMUSG00000026730.12",
                                 "ENSMUSG00000069132.3",
                                 "ENSMUSG00000041959.14",
                                 "ENSMUSG00000046999.2",
                                 "ENSMUSG00000079017.3",
                                 "ENSMUSG00000029361.18",
                                 "ENSMUSG00000031548.7",
                                 "ENSMUSG00000099061.1",
                                 "ENSMUSG00000064215.13",
                                 "ENSMUSG00000030677.8",
                                 "ENSMUSG00000030792.8",
                                 "ENSMUSG00000028862.6",
                                 "ENSMUSG00000034758.12",
                                 "ENSMUSG00000021070.6",
                                 "ENSMUSG00000040312.14",
                                 "ENSMUSG00000041380.13",
                                 "ENSMUSG00000050505.7",
                                 "ENSMUSG00000036902.11",
                                 "ENSMUSG00000097785.2",
                                 "ENSMUSG00000036545.9",
                                 "ENSMUSG00000059327.9",
                                 "ENSMUSG00000020524.16",
                                 "ENSMUSG00000048027.9",
                                 "ENSMUSG00000069072.9",
                                 "ENSMUSG00000039252.11",
                                 "ENSMUSG00000037143.17",
                                 "ENSMUSG00000029608.10",
                                 "ENSMUSG00000057729.12",
                                 "ENSMUSG00000041736.7",
                                 "ENSMUSG00000071379.2",
                                 "ENSMUSG00000000058.6",
                                 "ENSMUSG00000029552.19",
                                 "ENSMUSG00000028194.15",
                                 "ENSMUSG00000022342.6",
                                 "ENSMUSG00000039474.13",
                                 "ENSMUSG00000002980.14",
                                 "ENSMUSG00000020836.15",
                                 "ENSMUSG00000021319.7",
                                 "ENSMUSG00000023571.4",
                                 "ENSMUSG00000029765.12",
                                 "ENSMUSG00000056158.14",
                                 "ENSMUSG00000024044.19",
                                 "ENSMUSG00000004698.11",
                                 "ENSMUSG00000063887.13",
                                 "ENSMUSG00000019828.13",
                                 "ENSMUSG00000025790.14",
                                 "ENSMUSG00000042589.18",
                                 "ENSMUSG00000031772.17",
                                 "ENSMUSG00000026156.8",
                                 "ENSMUSG00000044816.10",
                                 "ENSMUSG00000026185.8",
                                 "ENSMUSG00000020728.17",
                                 "ENSMUSG00000059203.10",
                                 "ENSMUSG00000047904.6",
                                 "ENSMUSG00000020953.17",
                                 "ENSMUSG00000040998.18",
                                 "ENSMUSG00000061859.17",
                                 "ENSMUSG00000030551.14",
                                 "ENSMUSG00000027996.13",
                                 "ENSMUSG00000027577.14",
                                 "ENSMUSG00000046922.7",
                                 "ENSMUSG00000020000.7",
                                 "ENSMUSG00000090223.2",
                                 "ENSMUSG00000058897.18",
                                 "ENSMUSG00000034796.14",
                                 "ENSMUSG00000026147.16",
                                 "ENSMUSG00000026765.12",
                                 "ENSMUSG00000050447.15",
                                 "ENSMUSG00000028003.6",
                                 "ENSMUSG00000028584.3",
                                 "ENSMUSG00000027859.10",
                                 "ENSMUSG00000037771.11",
                                 "ENSMUSG00000027985.14",
                                 "ENSMUSG00000117655.1",
                                 "ENSMUSG00000022103.10",
                                 "ENSMUSG00000041449.16",
                                 "ENSMUSG00000056306.5",
                                 "ENSMUSG00000026787.3",
                                 "ENSMUSG00000046500.8",
                                 "ENSMUSG00000051397.5",
                                 "ENSMUSG00000030270.11",
                                 "ENSMUSG00000070880.10",
                                 "ENSMUSG00000036526.8",
                                 "ENSMUSG00000057818.8",
                                 "ENSMUSG00000018698.15",
                                 "ENSMUSG00000027168.21",
                                 "ENSMUSG00000042514.11",
                                 "ENSMUSG00000024985.20",
                                 "ENSMUSG00000043091.9",
                                 "ENSMUSG00000051747.15",
                                 "ENSMUSG00000086365.2",
                                 "ENSMUSG00000066363.12",
                                 "ENSMUSG00000029754.13",
                                 "ENSMUSG00000047746.14",
                                 "ENSMUSG00000058400.13")
                               ,]
rownames(genes_significant)<- c("Pitx2",
                                "Cga",
                                "Padi1",
                                "Myoz1",
                                "Ripor3",
                                "Fam166a",
                                "Irf6",
                                "Slc18a3",
                                "6330410L21Rik",
                                "Gm28723",
                                "Ntng2",
                                "Ccn3",
                                "Myh8",
                                "Esm1",
                                "Irx6",
                                "Slc44a4",
                                "Irs3",
                                "Sema3a",
                                "Nmbr",
                                "Scn4a",
                                "Ngfr",
                                "Plxnd1",
                                "Tll2",
                                "Zp1",
                                "Gm11872",
                                "Myh3",
                                "Ahnak2",
                                "Trh",
                                "Lamp5",
                                "C1ql3",
                                "Crabp2",
                                "Nccrp1",
                                "Nts",
                                "H2-Q7",
                                "H2-Q9",
                                "Kdf1",
                                "Ptgdr",
                                "Cidea",
                                "Tac1",
                                "Kcnk6",
                                "Pcsk5",
                                "Bcl11b",
                                "Ccdc88b",
                                "6530402F18Rik",
                                "Lncenc1",
                                "Slit2",
                                "Fgf10",
                                "Kit",
                                "BC043934",
                                "Cacna2d1",
                                "Nudt7",
                                "Pter",
                                "Nxph2",
                                "S100a10",
                                "1110032F04Rik",
                                "Ifi27l2a",
                                "Nos1",
                                "Sfrp1",
                                "Gm28050",
                                "Ifi27",
                                "Kif22",
                                "Dkkl1",
                                "Map3k6",
                                "Tle6",
                                "Bdkrb2",
                                "Cchcr1",
                                "Htr2c",
                                "Pcdh20",
                                "Neto2",
                                "B230217O12Rik",
                                "Adamts2",
                                "Eda",
                                "Gria1",
                                "Rgmb",
                                "Slc7a14",
                                "Lgi2",
                                "Cfap61",
                                "Rph3a",
                                "Prtn3",
                                "Tspo",
                                "Hpcal1",
                                "Cav2",
                                "Tes",
                                "Ddah1",
                                "Kcnv1",
                                "Wfs1",
                                "Bcam",
                                "Coro6",
                                "Sfrp4",
                                "C1qtnf12",
                                "Plxna4",
                                "Car10",
                                "Epb41l3",
                                "Hdac9",
                                "Nlgn1",
                                "Grm1",
                                "Slco3a1",
                                "Cux2",
                                "Cntnap4",
                                "B3gat2",
                                "D630023F18Rik",
                                "Igfbp5",
                                "Cep112",
                                "Il1rapl2",
                                "Sstr2",
                                "Coch",
                                "Npnt",
                                "Patj",
                                "Nr2f2",
                                "Sfrp2",
                                "Chrna4",
                                "Gpr6",
                                "Moxd1",
                                "Pcp4",
                                "Col25a1",
                                "Cpne7",
                                "Col9a1",
                                "Lypd6b",
                                "Lypd6",
                                "Lrat",
                                "Lrrc38",
                                "Ngf",
                                "Slc32a1",
                                "Lef1",
                                "6720468P15Rik",
                                "Gfra2",
                                "Serpina3h",
                                "Sertm1",
                                "Gad2",
                                "Tafa4",
                                "Tacstd2",
                                "Cpne9",
                                "Gad1",
                                "Card11",
                                "1600029O15Rik",
                                "Lhx1",
                                "Pax6",
                                "Klhl14",
                                "Tcf7l2",
                                "Tuba1c",
                                "Ttn",
                                "4930593C16Rik",
                                "Serpina3f",
                                "Dlx6",
                                "Fbxo40",
                                "Qrfpr")

# ordering genes
genes_significant <- genes_significant[order(rowMeans(genes_significant),
                                             decreasing = TRUE),]
# plotting heatmap
pheatmap::pheatmap(genes_significant, 
                   cluster_rows=FALSE, 
                   show_rownames=FALSE, 
                   show_colnames = TRUE,
                   cluster_cols=TRUE,
                   annotation_col = sampleTable)
