  #Set working directry.
setwd(dir = "Research Rheumtaoid Arhtritis")
getwd()
  #Install packages
install.packages('BiocManager')
BiocManager::install('Rsubread', force = TRUE)
library(Rsubread)
browseVignettes('Rsubread')

  #Making the index from the referencegene. Index is needed for the aligning.
buildindex(
  basename = 'GRCh38 data/GRCh38_index',
  reference = 'GRCh38 data/Homo_sapiens.GRCh38.dna.primary_assembly.fa',
  memory = 7000,
  indexSplit = TRUE)


  #Mapping of the sample against the human genome. 
align.norm1 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785819_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785819_2_subset40k.fastq", output_file = "norm1.BAM")
align.norm2 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785820_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785820_2_subset40k.fastq", output_file = "norm2.BAM")
align.norm3 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785828_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785828_2_subset40k.fastq", output_file = "norm3.BAM")
align.norm4 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785831_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785831_2_subset40k.fastq", output_file = "norm4.BAM")
align.RA1 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785979_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785979_2_subset40k.fastq", output_file = "RA1.BAM")
align.RA2 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785980_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785980_2_subset40k.fastq", output_file = "RA2.BAM")
align.RA3 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785986_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785986_2_subset40k.fastq", output_file = "RA3.BAM")
align.RA4 <- align(index = "GRCh38_index", readfile1 = "Data_RA_raw/SRR4785988_1_subset40k.fastq", readfile2 = "Data_RA_raw/SRR4785988_2_subset40k.fastq", output_file = "RA4.BAM")


  #Loading the Rsamtools for sorting and indexing.
library(Rsamtools)

  #Filenames of the samples
samples <- c('norm1', 'norm2', 'norm3', 'norm4', 'RA1', 'RA2', 'RA3', 'RA4')


  #Sorting and indexing of the BAM-file for each sample.
lapply(samples, function(s) {sortBam(file = paste0(s, '.BAM'), destination = paste0(s, '.sorted'))
})

  #Creating BAM.bai
lapply(samples, function(s) {
  bam_file <- paste0(s, ".BAM")
  sorted_file_prefix <- paste0(s, ".sorted")
  sortBam(file = bam_file, destination = sorted_file_prefix)
  indexBam(paste0(sorted_file_prefix, ".bam"))})

library(readr)
library(dplyr)
library(Rsamtools)
library(Rsubread)

  #Define a vector with the names of the BAM-files. Each BAM contains reads of an RNA-seq-experiment.
allsamples <- c('norm1.sorted.bam', 
                'norm2.sorted.bam', 
                'norm3.sorted.bam', 
                'norm4.sorted.bam', 
                'RA1.sorted.bam',
                'RA2.sorted.bam',
                'RA3.sorted.bam',
                'RA4.sorted.bam')

  #Creating a count matrix.
count_matrix <- featureCounts(
  files = allsamples,
  annot.ext = "Homo_sapiens.GRCh38.114.gtf",
  isPairedEnd = TRUE,
  isGTFAnnotationFile = TRUE,
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE)

count_matrix <- read.table("count_matrix_RA_casus.txt", header = TRUE, row.names = 1)

head(count_matrix)
dim(count_matrix)
class(count_matrix)

colnames(count_matrix) <- c('Norm1', 'Norm2', 'Norm3', 'Norm4', 'RA1', 'RA2', 'RA3', 'RA4')

write.csv(count_matrix, "Edited_countmatrix.casusRA.csv")

treatment <- c("Norm", "Norm", "Norm", "Norm", "RA", "RA", "RA", "RA")
treatment_table <- data.frame(treatment)
rownames(treatment_table) <- c('Norm1', 'Norm2', 'Norm3', 'Norm4', 'RA1', 'RA2', 'RA3', 'RA4')

BiocManager::install('DESeq2')
library(DESeq2)
BiocManager::install('KEGGREST')
library(KEGGREST)


  #Creating a DESeqDataSet.
dds <- DESeqDataSetFromMatrix(countData = round(count_matrix),
                              colData = treatment_table,
                              design = ~ treatment)

  #Analysis.
dds <- DESeq(dds)
results_RA <- results(dds)

  #Saving results in a file
write.table(results_RA, file = 'Results_case_RA.csv', row.names = TRUE, col.names = TRUE)

  #Filter up-and down-regulating processes.
sum(results_RA$padj < 0.05 & results_RA$log2FoldChange > 1, na.rm = TRUE)
sum(results_RA$padj < 0.05 & results_RA$log2FoldChange < -1, na.rm = TRUE)

highest_fold_change <- results_RA[order(results_RA$log2FoldChange, decreasing = TRUE), ]
lowest_fold_change <- results_RA[order(results_RA$log2FoldChange, decreasing = FALSE), ]
lowest_p_waarde <- results_RA[order(results_RA$padj, decreasing = FALSE), ]

head(lowest_fold_change)
head(lowest_p_waarde)
head(highest_fold_change)

#Creating a volcano plot.
if (!requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  BiocManager::install("EnhancedVolcano")
}
library(EnhancedVolcano)

EnhancedVolcano(results_RA,
                lab = rownames(results_RA),
                x = 'log2FoldChange',
                y = 'padj')

  #Alternative plot without p-value cutoff, all genes.
EnhancedVolcano(results_RA,
                lab = rownames(results_RA),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0)

dev.copy(png, 'Volcanoplot_casus_RA.png', 
         width = 8,
         height = 10,
         units = 'in',
         res = 500)
dev.off()

  #Creating a hsa pathway. 
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

results_pathview <- results_RA
results_pathview[1] <- NULL
results_pathview[2:5] <- NULL
results_pathview

pathview(
  gene.data = results_pathview,
  pathway.id = "hsa05323",  
  species = "hsa",          
  gene.idtype = "SYMBOL",     
  limit = list(gene = 2),
  low = 'red',
  high = 'green'
)

head(results_pathview)

keggFind("pathway", "rheumatoid arthritis")

keggLink("pathway", "hsa:283149")