# This script was used to import quantification files into R envrinment via tximport and
# to do differential gene expression analysis (DEG) on them via deseq2.
# This script was used for analysing all 12 samples.
# Author(s): Dorottya Ralbovszki
# Created:2020.01.20.

rm(list = ls()) #Clear workspace

# importing sample info of samples
sampleinfo <- read.delim("metadataall.txt")
View(sampleinfo)
sampleinfo

# importing quantification files
dir <- list.files("salmon_data/all/")
quant_files <- list.files("salmon_data/all/", pattern = "quant.sf", recursive = TRUE, full.names = TRUE)
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
write.csv(tx2gene, file = "tx2gene_multi.csv", row.names = FALSE, quote = FALSE)
txi <- tximport(quant_files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = FALSE, ignoreAfterBar = TRUE, countsFromAbundance = "lengthScaledTPM")
table(tx_map$TXNAME %in% quants$Name)

names(txi)
head(txi$counts)

# DEG analysis using deseq2
library(DESeq2)
dds1 <- DESeqDataSetFromTximport(txi, colData = sampleinfo, design <- ~ Genotype + Strain + Genotype:Strain)

# exporting TPM
tpm <- txi$abudance
write.csv(tpm, file = "tmp_values_multi.csv", quote = FALSE)

# cheking in how many samples genes are expressed
is_expressed <- assay(dds) >= 5
head(is_expressed)
sum(is_expressed[1,])
sum(is_expressed[2,])
hist(rowSums(is_expressed),main="Number of samples a gene is expressed in",xlab="Sample Count")

# filtering out genes that had a lower read number than 5 when the read number of all samples (6) was summarized
keep <- rowSums(counts(dds1)) >= 5
dds1 <- dds1[keep,]

# visualising count distributions
boxplot(assay(dds))
boxplot(log10(assay(dds)))

# setting the right comparison conditions
dds1$Genotype = relevel(dds1$Genotype,"wildtype")
dds1$Strain = relevel(dds1$Strain, "SWISSOF1")
# DEG with the new comparison settings
dds2 <- DESeq(dds1)
# extracting result table from DEG analysis
res_multi <- results(dds2, contrast=list(c("Genotype_mutant_vs_wildtype", 
                                      "Genotypemutant.StrainC57")))

# checking if comparison condition was right
resultsNames(dds2)


# plot counts of smallest p value
plotCounts(dds, gene=which.min(res_multi$padj), intgroup=c("Strain", "Genotype"))

# remove string after period to get actual ENSEMBL ID
tmp = gsub("\\..*","",row.names(res_multi)) 
row.names(res_multi) = tmp
head(row.names(res_multi))

# order the results by p values
res_multiOrdered <- res_multi[order(res_multi$pvalue),]
# save separetaly the down-, and upregulated genes
resup_multi <- subset(res_multi, log2FoldChange>0)
resdown_multi <- subset(res_multi, log2FoldChange<0)

# getting the number of significant genes
sum(res_multi$padj < 0.05, na.rm=TRUE)

# summary of analysis
summary(res_multi)
res05 <- results(dds, alpha=0.05)
summary(res05)

# filtering out too high log2fold changes which means that a gene was only expressed in 1 sample/strain
keep_logfold_p <- res_multiOrdered$log2FoldChange <= 10
res_multiOrdered <- res_multiOrdered[keep_logfold_p,]
keep_logfold_n <- res_multiOrdered$log2FoldChange >= -10
res_multiOrdered <- res_multiOrdered[keep_logfold_n,]

# checking filtered result table
res_multiOrdered

# export results into a CSV file
write.csv( as.data.frame(res_multiOrdered), file="res_multi0508.csv" )

# getting the number of significant genes after filtering
sum(res_multiOrdered$padj < 0.05, na.rm=TRUE)

# annotating result table with gene sybol, entrez ID and gene name
library("AnnotationDbi")
library("org.Mm.eg.db")

annots_symbol <- select(org.Mm.eg.db, keys = rownames(res_multiOrdered), column = "SYMBOL", keytype = "ENSEMBL")
annots_entrez <- select(org.Mm.eg.db, keys = rownames(res_multiOrdered), column = "ENTREZID", keytype = "ENSEMBL")
annots_name <- select(org.Mm.eg.db, keys = rownames(res_multirdered), column = "GENENAME", keytype = "ENSEMBL")

# exporting annotated results into a csv file
write.csv( as.data.frame(annots_name), file="annots_name_c57.csv" )
write.csv( as.data.frame(annots_entrez), file="annots_entrez_c57.csv" )
write.csv( as.data.frame(annots_symbol), file="annots_symbol_multi0507.csv" )



# log fold change shrinkage for visualization and ranking
resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Strain_C57_vs_SWISSOF1", type="apeglm")
resLFC



# MA-plot showing the log2 fold changes attributable to a given variable over the mean of normalized counts for all the samples 
# coloured points showing p-value less than 0.1
# comparint the raw and the log fold change shrinkage
plotMA(res_multi, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# variance stabilizing transformations (VST)
vsd <- vst(dds2, blind=FALSE)

# plot PCA of the transformed data
plotPCA(vsd, intgroup=c("Strain", "Genotype"))


# creating significant genes heatmap
# only worked if sample metadata was created this way
sampleTable <- data.frame(Genotype = factor(rep(c("WT", "MUT", "WT", "MUT"), 
                                                     each = 3)), Strain = factor(rep(c("SWISSOF1", "C57"), 
                                                                                     each = 6)))
# testing heatmap plotting
rownames(sampleTable) <- colnames(txi$counts)
sampleTable

# extracting test genes
expmatrix_DESeq <- DESeq2::rlog(dds2, fitType="local")
expmatrix <- SummarizedExperiment::assay(expmatrix_DESeq)
select <- order(rowMeans(expmatrix),decreasing=TRUE)[1:30]

library(gplots)
heatmap(expmatrix[select,])
genes_interest <- expmatrix[c("ENSMUSG00000069917.7",
                              "ENSMUSG00000073940.3",
                              "ENSMUSG00000069919.7",
                              "ENSMUSG00000052305.6",
                              "ENSMUSG00000038600.12",
                              "ENSMUSG00000016427.7",
                              "ENSMUSG00000040280.10",
                              "ENSMUSG00000053930.13",
                              "ENSMUSG00000020950.10",
                              "ENSMUSG00000038738.15")
                            ,]
rownames(genes_interest)<- c("Hba-a2","Hbb-bt","Hba-a1","Hbb-bs",
                             "Atp6v0a4","Ndufa1","Ndufa4l2","Shisa6",
                             "Foxg1","Shank1")

# ordering test genes
genes_interest <- genes_interest[order(rowMeans(genes_interest),
                                       decreasing = TRUE),]

# plotting test heatmap
pheatmap::pheatmap(genes_interest, 
                   cluster_rows=FALSE, 
                   show_rownames=TRUE, 
                   show_colnames = TRUE,
                   cluster_cols=TRUE,
                   annotation_col = sampleTable,
                   clustering_method = "median")



# my heatmap containing significant genes from the 4 comparisons
genes_interest <- expmatrix[c("ENSMUSG00000000058.6",
                              "ENSMUSG00000000120.6",
                              "ENSMUSG00000000308.14",
                              "ENSMUSG00000000384.15",
                              "ENSMUSG00000000805.18",
                              "ENSMUSG00000001027.7",
                              "ENSMUSG00000001240.13",
                              "ENSMUSG00000001300.16",
                              "ENSMUSG00000001930.17",
                              "ENSMUSG00000001946.14",
                              "ENSMUSG00000002980.14",
                              "ENSMUSG00000003476.16",
                              "ENSMUSG00000003929.11",
                              "ENSMUSG00000004221.16",
                              "ENSMUSG00000004698.11",
                              "ENSMUSG00000004885.5",
                              "ENSMUSG00000004895.9",
                              "ENSMUSG00000004936.8",
                              "ENSMUSG00000005299.6",
                              "ENSMUSG00000005447.12",
                              "ENSMUSG00000005672.12",
                              "ENSMUSG00000005892.4",
                              "ENSMUSG00000005973.6",
                              "ENSMUSG00000007034.15",
                              "ENSMUSG00000009734.18",
                              "ENSMUSG00000012483.4",
                              "ENSMUSG00000012519.14",
                              "ENSMUSG00000013523.13",
                              "ENSMUSG00000016427.7",
                              "ENSMUSG00000017929.13",
                              "ENSMUSG00000017978.18",
                              "ENSMUSG00000018470.8",
                              "ENSMUSG00000018698.15",
                              "ENSMUSG00000019828.13",
                              "ENSMUSG00000019865.9",
                              "ENSMUSG00000019874.11",
                              "ENSMUSG00000019890.4",
                              "ENSMUSG00000020000.7",
                              "ENSMUSG00000020142.12",
                              "ENSMUSG00000020173.17",
                              "ENSMUSG00000020218.11",
                              "ENSMUSG00000020253.15",
                              "ENSMUSG00000020524.16",
                              "ENSMUSG00000020723.3",
                              "ENSMUSG00000020728.17",
                              "ENSMUSG00000020836.15",
                              "ENSMUSG00000020908.14",
                              "ENSMUSG00000020950.10",
                              "ENSMUSG00000020953.17",
                              "ENSMUSG00000021070.6",
                              "ENSMUSG00000021193.10",
                              "ENSMUSG00000021319.7",
                              "ENSMUSG00000021337.8",
                              "ENSMUSG00000021567.15",
                              "ENSMUSG00000021587.5",
                              "ENSMUSG00000021647.7",
                              "ENSMUSG00000021732.14",
                              "ENSMUSG00000021867.16",
                              "ENSMUSG00000021872.8",
                              "ENSMUSG00000021880.7",
                              "ENSMUSG00000021969.8",
                              "ENSMUSG00000022048.8",
                              "ENSMUSG00000022090.10",
                              "ENSMUSG00000022103.10",
                              "ENSMUSG00000022342.6",
                              "ENSMUSG00000022667.18",
                              "ENSMUSG00000023439.11",
                              "ENSMUSG00000023571.4",
                              "ENSMUSG00000024014.8",
                              "ENSMUSG00000024044.19",
                              "ENSMUSG00000024140.10",
                              "ENSMUSG00000024526.9",
                              "ENSMUSG00000024565.10",
                              "ENSMUSG00000024681.11",
                              "ENSMUSG00000024713.16",
                              "ENSMUSG00000024734.8",
                              "ENSMUSG00000024942.17",
                              "ENSMUSG00000024985.20",
                              "ENSMUSG00000025013.15",
                              "ENSMUSG00000025329.3",
                              "ENSMUSG00000025362.6",
                              "ENSMUSG00000025789.9",
                              "ENSMUSG00000025790.14",
                              "ENSMUSG00000025804.5",
                              "ENSMUSG00000025870.10",
                              "ENSMUSG00000026018.12",
                              "ENSMUSG00000026048.16",
                              "ENSMUSG00000026098.13",
                              "ENSMUSG00000026113.17",
                              "ENSMUSG00000026147.16",
                              "ENSMUSG00000026156.8",
                              "ENSMUSG00000026185.8",
                              "ENSMUSG00000026237.5",
                              "ENSMUSG00000026413.12",
                              "ENSMUSG00000026516.8",
                              "ENSMUSG00000026638.15",
                              "ENSMUSG00000026688.5",
                              "ENSMUSG00000026730.12",
                              "ENSMUSG00000026765.12",
                              "ENSMUSG00000026787.3",
                              "ENSMUSG00000026969.3",
                              "ENSMUSG00000027168.21",
                              "ENSMUSG00000027217.13",
                              "ENSMUSG00000027224.14",
                              "ENSMUSG00000027270.14",
                              "ENSMUSG00000027274.16",
                              "ENSMUSG00000027400.11",
                              "ENSMUSG00000027570.15",
                              "ENSMUSG00000027577.14",
                              "ENSMUSG00000027678.17",
                              "ENSMUSG00000027792.11",
                              "ENSMUSG00000027859.10",
                              "ENSMUSG00000027985.14",
                              "ENSMUSG00000027996.13",
                              "ENSMUSG00000028003.6",
                              "ENSMUSG00000028023.16",
                              "ENSMUSG00000028194.15",
                              "ENSMUSG00000028234.6",
                              "ENSMUSG00000028298.10",
                              "ENSMUSG00000028370.7",
                              "ENSMUSG00000028487.18",
                              "ENSMUSG00000028558.14",
                              "ENSMUSG00000028584.3",
                              "ENSMUSG00000028602.12",
                              "ENSMUSG00000028635.7",
                              "ENSMUSG00000028656.14",
                              "ENSMUSG00000028757.4",
                              "ENSMUSG00000028841.14",
                              "ENSMUSG00000028862.6",
                              "ENSMUSG00000028883.17",
                              "ENSMUSG00000028901.13",
                              "ENSMUSG00000029005.4",
                              "ENSMUSG00000029193.7",
                              "ENSMUSG00000029288.11",
                              "ENSMUSG00000029361.18",
                              "ENSMUSG00000029428.13",
                              "ENSMUSG00000029552.19",
                              "ENSMUSG00000029608.10",
                              "ENSMUSG00000029754.13",
                              "ENSMUSG00000029765.12",
                              "ENSMUSG00000029917.15",
                              "ENSMUSG00000030123.15",
                              "ENSMUSG00000030235.17",
                              "ENSMUSG00000030270.11",
                              "ENSMUSG00000030279.15",
                              "ENSMUSG00000030413.7",
                              "ENSMUSG00000030532.6",
                              "ENSMUSG00000030551.14",
                              "ENSMUSG00000030677.8",
                              "ENSMUSG00000030761.16",
                              "ENSMUSG00000030792.8",
                              "ENSMUSG00000031212.3",
                              "ENSMUSG00000031297.14",
                              "ENSMUSG00000031391.18",
                              "ENSMUSG00000031548.7",
                              "ENSMUSG00000031558.15",
                              "ENSMUSG00000031738.14",
                              "ENSMUSG00000031767.13",
                              "ENSMUSG00000031772.17",
                              "ENSMUSG00000031997.9",
                              "ENSMUSG00000032076.20",
                              "ENSMUSG00000032259.8",
                              "ENSMUSG00000032271.13",
                              "ENSMUSG00000032368.14",
                              "ENSMUSG00000032643.12",
                              "ENSMUSG00000032679.12",
                              "ENSMUSG00000032854.12",
                              "ENSMUSG00000033597.9",
                              "ENSMUSG00000033808.16",
                              "ENSMUSG00000033960.6",
                              "ENSMUSG00000034055.16",
                              "ENSMUSG00000034243.17",
                              "ENSMUSG00000034652.12",
                              "ENSMUSG00000034723.11",
                              "ENSMUSG00000034758.12",
                              "ENSMUSG00000034796.14",
                              "ENSMUSG00000034892.8",
                              "ENSMUSG00000035277.15",
                              "ENSMUSG00000035329.7",
                              "ENSMUSG00000035513.19",
                              "ENSMUSG00000035726.8",
                              "ENSMUSG00000035929.11",
                              "ENSMUSG00000036131.12",
                              "ENSMUSG00000036526.8",
                              "ENSMUSG00000036545.9",
                              "ENSMUSG00000036902.11",
                              "ENSMUSG00000037025.11",
                              "ENSMUSG00000037143.17",
                              "ENSMUSG00000037362.8",
                              "ENSMUSG00000037400.17",
                              "ENSMUSG00000037526.7",
                              "ENSMUSG00000037600.16",
                              "ENSMUSG00000037771.11",
                              "ENSMUSG00000037784.14",
                              "ENSMUSG00000037962.7",
                              "ENSMUSG00000037990.18",
                              "ENSMUSG00000038007.14",
                              "ENSMUSG00000038173.15",
                              "ENSMUSG00000038257.9",
                              "ENSMUSG00000038291.16",
                              "ENSMUSG00000038600.12",
                              "ENSMUSG00000038738.15",
                              "ENSMUSG00000039106.6",
                              "ENSMUSG00000039126.10",
                              "ENSMUSG00000039231.18",
                              "ENSMUSG00000039252.11",
                              "ENSMUSG00000039474.13",
                              "ENSMUSG00000039488.15",
                              "ENSMUSG00000039579.15",
                              "ENSMUSG00000039672.12",
                              "ENSMUSG00000039706.11",
                              "ENSMUSG00000039714.9",
                              "ENSMUSG00000039735.16",
                              "ENSMUSG00000039977.16",
                              "ENSMUSG00000040118.15",
                              "ENSMUSG00000040312.14",
                              "ENSMUSG00000040543.16",
                              "ENSMUSG00000040998.18",
                              "ENSMUSG00000041380.13",
                              "ENSMUSG00000041449.16",
                              "ENSMUSG00000041559.7",
                              "ENSMUSG00000041607.17",
                              "ENSMUSG00000041736.7",
                              "ENSMUSG00000041773.8",
                              "ENSMUSG00000041911.3",
                              "ENSMUSG00000041959.14",
                              "ENSMUSG00000041975.17",
                              "ENSMUSG00000042369.8",
                              "ENSMUSG00000042379.8",
                              "ENSMUSG00000042501.12",
                              "ENSMUSG00000042514.11",
                              "ENSMUSG00000042589.18",
                              "ENSMUSG00000042770.8",
                              "ENSMUSG00000042772.15",
                              "ENSMUSG00000043091.9",
                              "ENSMUSG00000043671.14",
                              "ENSMUSG00000044068.7",
                              "ENSMUSG00000044566.15",
                              "ENSMUSG00000044708.5",
                              "ENSMUSG00000044816.10",
                              "ENSMUSG00000045573.9",
                              "ENSMUSG00000046410.10",
                              "ENSMUSG00000046480.6",
                              "ENSMUSG00000046500.8",
                              "ENSMUSG00000046610.15",
                              "ENSMUSG00000046922.7",
                              "ENSMUSG00000046999.2",
                              "ENSMUSG00000047182.6",
                              "ENSMUSG00000047586.4",
                              "ENSMUSG00000047746.14",
                              "ENSMUSG00000047766.15",
                              "ENSMUSG00000047810.9",
                              "ENSMUSG00000047904.6",
                              "ENSMUSG00000048027.9",
                              "ENSMUSG00000048251.15",
                              "ENSMUSG00000049336.16",
                              "ENSMUSG00000049630.6",
                              "ENSMUSG00000049744.15",
                              "ENSMUSG00000049928.15",
                              "ENSMUSG00000050148.9",
                              "ENSMUSG00000050447.15",
                              "ENSMUSG00000050505.7",
                              "ENSMUSG00000050558.13",
                              "ENSMUSG00000050711.7",
                              "ENSMUSG00000051246.3",
                              "ENSMUSG00000051397.5",
                              "ENSMUSG00000051747.15",
                              "ENSMUSG00000052305.6",
                              "ENSMUSG00000052926.16",
                              "ENSMUSG00000053310.11",
                              "ENSMUSG00000053930.13",
                              "ENSMUSG00000054409.5",
                              "ENSMUSG00000054457.5",
                              "ENSMUSG00000055202.11",
                              "ENSMUSG00000055301.8",
                              "ENSMUSG00000055675.6",
                              "ENSMUSG00000055775.16",
                              "ENSMUSG00000056158.14",
                              "ENSMUSG00000056306.5",
                              "ENSMUSG00000056418.3",
                              "ENSMUSG00000056596.8",
                              "ENSMUSG00000056608.14",
                              "ENSMUSG00000057123.14",
                              "ENSMUSG00000057729.12",
                              "ENSMUSG00000057818.8",
                              "ENSMUSG00000058174.7",
                              "ENSMUSG00000058400.13",
                              "ENSMUSG00000058443.5",
                              "ENSMUSG00000058897.18",
                              "ENSMUSG00000059040.5",
                              "ENSMUSG00000059203.10",
                              "ENSMUSG00000059325.14",
                              "ENSMUSG00000059327.9",
                              "ENSMUSG00000059839.9",
                              "ENSMUSG00000060063.9",
                              "ENSMUSG00000060550.16",
                              "ENSMUSG00000060860.8",
                              "ENSMUSG00000061414.8",
                              "ENSMUSG00000061762.12",
                              "ENSMUSG00000061859.17",
                              "ENSMUSG00000063171.4",
                              "ENSMUSG00000063260.2",
                              "ENSMUSG00000063698.9",
                              "ENSMUSG00000063887.13",
                              "ENSMUSG00000064215.13",
                              "ENSMUSG00000064329.13",
                              "ENSMUSG00000064330.9",
                              "ENSMUSG00000066361.3",
                              "ENSMUSG00000066363.12",
                              "ENSMUSG00000066438.6",
                              "ENSMUSG00000066705.7",
                              "ENSMUSG00000067288.13",
                              "ENSMUSG00000067870.5",
                              "ENSMUSG00000068117.10",
                              "ENSMUSG00000068396.9",
                              "ENSMUSG00000068523.12",
                              "ENSMUSG00000068697.7",
                              "ENSMUSG00000069072.9",
                              "ENSMUSG00000069132.3",
                              "ENSMUSG00000069917.7",
                              "ENSMUSG00000069919.7",
                              "ENSMUSG00000070056.6",
                              "ENSMUSG00000070583.1",
                              "ENSMUSG00000070605.4",
                              "ENSMUSG00000070880.10",
                              "ENSMUSG00000071369.11",
                              "ENSMUSG00000071379.2",
                              "ENSMUSG00000071470.4",
                              "ENSMUSG00000071489.1",
                              "ENSMUSG00000072437.4",
                              "ENSMUSG00000072812.4",
                              "ENSMUSG00000073876.3",
                              "ENSMUSG00000073940.3",
                              "ENSMUSG00000073982.11",
                              "ENSMUSG00000074269.10",
                              "ENSMUSG00000074577.9",
                              "ENSMUSG00000074695.3",
                              "ENSMUSG00000074731.3",
                              "ENSMUSG00000074735.2",
                              "ENSMUSG00000075296.5",
                              "ENSMUSG00000075330.4",
                              "ENSMUSG00000075705.12",
                              "ENSMUSG00000076498.2",
                              "ENSMUSG00000078503.9",
                              "ENSMUSG00000078591.1",
                              "ENSMUSG00000078735.4",
                              "ENSMUSG00000078952.9",
                              "ENSMUSG00000078954.9",
                              "ENSMUSG00000079017.3",
                              "ENSMUSG00000079018.10",
                              "ENSMUSG00000079499.9",
                              "ENSMUSG00000079588.3",
                              "ENSMUSG00000079641.3",
                              "ENSMUSG00000079685.10",
                              "ENSMUSG00000082361.6",
                              "ENSMUSG00000086017.1",
                              "ENSMUSG00000086365.2",
                              "ENSMUSG00000086600.8",
                              "ENSMUSG00000087369.1",
                              "ENSMUSG00000089679.1",
                              "ENSMUSG00000090223.2",
                              "ENSMUSG00000091705.8",
                              "ENSMUSG00000092116.1",
                              "ENSMUSG00000093483.1",
                              "ENSMUSG00000094686.1",
                              "ENSMUSG00000095595.2",
                              "ENSMUSG00000096449.2",
                              "ENSMUSG00000096995.2",
                              "ENSMUSG00000097039.8",
                              "ENSMUSG00000097431.2",
                              "ENSMUSG00000097462.7",
                              "ENSMUSG00000097622.2",
                              "ENSMUSG00000097785.2",
                              "ENSMUSG00000099061.1",
                              "ENSMUSG00000100241.1",
                              "ENSMUSG00000100627.6",
                              "ENSMUSG00000101028.1",
                              "ENSMUSG00000101969.2",
                              "ENSMUSG00000102422.1",
                              "ENSMUSG00000102644.5",
                              "ENSMUSG00000104178.1",
                              "ENSMUSG00000105960.1",
                              "ENSMUSG00000105987.4",
                              "ENSMUSG00000111785.1",
                              "ENSMUSG00000113771.1",
                              "ENSMUSG00000115529.1",
                              "ENSMUSG00000116819.1",
                              "ENSMUSG00000117655.1")
                            ,]
genes_interest <- genes_interest[order(rowMeans(genes_interest),
                                       decreasing = TRUE),]
# plotting heatmap
pheatmap::pheatmap(genes_interest, 
                   cluster_rows=FALSE, 
                   show_rownames=FALSE, 
                   show_colnames = TRUE,
                   cluster_cols=TRUE,
                   annotation_col = sampleTable)

# exporting dds2 into csv file
write.csv( as.data.frame(counts(dds2)), file="dds_multi0507.csv" )