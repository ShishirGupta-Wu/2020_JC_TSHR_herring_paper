#########################################################################################################
#                                                        
# R script that creates boxplots corresponding to the relative expression (in counts-per-million) of 
# 15 genes of interest in brain (BR) and saccus vasculosus (BSH) tissues of 14 Atlantic herring samples
#
# The RNA-seq data was mapped to the genome Ch_v.2.0.2 using the aligner GSNAP in SNP-tolerant mode. 
# Read counts were assigned to genes using the program featureCounts. Raw read counts were normalized 
# using the TMM-method implemented in the R package edgeR. 
#
# Created by: Angela P. Fuentes-Pardo, PhD. E-mail: apfuentesp@gmail.com
# Postdoctoral researcher at Leif Andersson Lab, Uppsala University
# Created on: 2020-01-15, updated on: 2020-08-26
#                                                        
#########################################################################################################

# Clean environment space
rm(list = ls())

# --------------------------------------------------------------------------------
# Set general input files and parameters
# --------------------------------------------------------------------------------
workingDir <- "~/Dropbox/PostDoc_UU/Projects/Herring/RNA-seq/DGE_genes_interest/TSHR_Junfeng/geneExpression-TSHR-others/2020-08-26-BSH-BR-tissue"
# rc_Ch_v2.0.2_gsnap_BSHBR_fc
rc <- "~/Dropbox/PostDoc_UU/Projects/Herring/RNA-seq/DGE_genes_interest/TSHR_Junfeng/data/featureCounts/Angela/Her-16samples-RNAseq-BSH-BR.gsnap.s2.q20.gene_id.featureCounts_edited.txt"
# rc_Ch_v2.0.2_star_BSHBR_fc
#rc <- "~/Dropbox/PostDoc_UU/Projects/Herring/RNA-seq/DGE_genes_interest/TSHR_Junfeng/data/featureCounts/Angela/merged_gene_counts_2020-01-13-Ch_v2.0.2_star_BSHBR_edited.txt"
metaData <- "~/Dropbox/PostDoc_UU/Projects/Herring/RNA-seq/Spawning/analysis/04-DGE/metaData_Spawning_BR-BSH_RNAseq_Herring.csv"
excludeSamples <- c("Ber16.T17F.BSH", "HAS.5F.BR")  # Because of low TIN value, Nima
sampleOrder <- c("Ber16.T10F.2.BR","Ber16.T13F.3.BR","Ber16.T17F.BR","Ber16.T18F.BR","Ber16.T10F.BSH","Ber16.T13F.BSH","Ber16.T18F.BSH","HAS.8F.BR","HAS.9F.BR","HAS.14F.BR","HAS.5F.BSH","HAS.8F.BSH","HAS.9F.BSH","HAS.14F.BSH")
minCPMvalue <- 1
minCPMsampleNumber <- 1
genesOfInterest <- "~/Dropbox/PostDoc_UU/Projects/Herring/RNA-seq/DGE_genes_interest/TSHR_Junfeng/geneExpression-TSHR-others/2020-08-26-BSH-BR-tissue/list_genes_interest_2020-08-26.txt"
plotIntermediates <- FALSE
alignerName <- "GSNAP"  # STAR
genomeVersion <- "Ch_v2.0.2"

# --------------------------------------------------------------------------------
# Data format preprocessing
# --------------------------------------------------------------------------------

# Set working directory
setwd(workingDir)

### Load the count table ---------------- 
# The count data shows read counts across samples and genes. The columns denote samples and rows denote genes
cr <- read.table(rc, header = TRUE, row.names = 1, sep = "\t", stringsAsFactors = FALSE, quote = "", fill = FALSE)
head(cr)
str(cr)

# Save gene lengths as a separate dataframe and eliminate this column from the read count dataframe (tequired for rpkm calculation)
#crGeneLength <- data.frame(ensemble_gene_id = rownames(cr), Length = cr$Length, stringsAsFactors = FALSE)
#cr$Length <- NULL
#head(crGeneLength)
#colnames(cr)
#dim(cr)

### Load the metaData file (useful for normalization, DGE, but not required in this case) ----------------
# Each row corresponds to a sample, and each column to sample attributes such tissue type, sequencing batch, lane, etc.
mr <- read.csv2(metaData, header = TRUE, stringsAsFactors = FALSE)
rownames(mr) <- mr$Sample_ID  # Assign sample IDs as rownames
head(mr)
str(mr)

# Exclude samples if required (e.g. because of low TIN value, or poor RIN)
if(!is.null(excludeSamples)) {
  cr <- cr[ ,!(names(cr) %in% excludeSamples)]
  colnames(cr)
  mr <- mr[!(rownames(mr) %in% excludeSamples), ]
  rownames(mr)
}

# Reorder samples (columns) if required
if(!is.null(sampleOrder)) {
  cr <- cr[ ,sampleOrder]
  colnames(cr)
  mr <- mr[sampleOrder, ]
  rownames(mr)
}

# --------------------------------------------------------------------------------
# Raw read count filtering
# --------------------------------------------------------------------------------

# As dataset were subsetted, assign them to the correspondent df.
cr_subset <- cr
mr_subset <- mr

### Create a boxplot to visualise the distribution of raw counts
if(plotIntermediates == TRUE) {
  boxplot(log10(as.matrix(cr_subset) + 1), ylab = expression("Log"[10]~"Read counts"), las = 2,
          main = "Raw data")
}
# When the median values are zero across all samples, the data set would benefit from a low count filtering.

### Check if any samples need to be discarded based on the number of genes detected 
# Create a barplot of genes detected across samples.
if(plotIntermediates == TRUE) { 
  barplot(colSums(cr_subset > 5), ylab = "Number of detected genes", las = 2)
  abline(h = median(colSums(cr_subset > 5))) 
}

# Create a similar plot for detection rate across genes.
if(plotIntermediates == TRUE) { barplot(rowSums(cr_subset > 5), xlab = "Genes", ylab = "Number of samples", las = 2, names.arg = "")
  abline(h = median(rowSums(cr_subset > 5))) 
}
# Some of the genes are not detected across most samples. These genes can be discarded 

### Remove genes with low counts ----------------
# A minimum of 1 count-per-million (CPM) across 1 sample was required (for DGE, 2 samples is generally used since each of test groups consist of 3 samples).
library(edgeR)

keepgenes <- rowSums(edgeR::cpm(cr_subset) > minCPMvalue) >= minCPMsampleNumber
cf <- cr_subset[keepgenes, ]  # The missingness in the data set is reduced

if(plotIntermediates == TRUE) {
  if(logBase == 10) {
    boxplot(log10(as.matrix(cf) + 1), ylab = expression("Log"[10]~"Read counts"), las = 2,
            main = paste0("Raw data filtered - minimum ", minCPMvalue," CPM across ", minCPMsampleNumber ," samples"))
  }
  if(logBase == 2) {
    boxplot(log2(as.matrix(cf) + 1), ylab = expression("Log"[2]~"Read counts"), las = 2,
            main = paste0("Raw data filtered - minimum ", minCPMvalue," CPM across ", minCPMsampleNumber ," samples"))
  }
}

# Explore how many genes were removed.
length(row.names(cr_subset))  # Number of genes before
#[1] 27879
length(row.names(cf))  # Number of genes after filtering
#[1] 21171
length(row.names(cr_subset)) - length(row.names(cf))
#[1] 6708 genes were removed

# Check that the labels are the same order in counts and metaData
all.equal(colnames(cf), rownames(mr_subset))
mr_subset <- mr_subset[match(colnames(cf),mr_subset$Sample_ID), ]
all.equal(colnames(cf), rownames(mr_subset))
# Save the filtered data
#write.table(cf,"./counts_filtered.txt", col.names=T, quote=F, sep="\t", dec=".")


# --------------------------------------------------------------------------------
# Calculate CPM-TMM-normalized counts with edgeR
# --------------------------------------------------------------------------------

### Obtain normalised expression values using the TMM method implemented in edgeR ----------------
# The TMM (trimmed mean of M-values) method uses a weighted trimmed mean of the log expression ratios between samples

# Set the factor giving group membership for columns of y.
group <- factor(colnames(cf))  # in this case each sample is a group, as we are not calculating differential expression between experimental groups
#group <- factor(mr_subset$Tissue)  # set tissues as groups
# Create a DGEList object from the read count data
expr <- edgeR::DGEList(counts = cf, group = group)
# Obtain library size normalization factors using the TMM method
normfact <- edgeR::calcNormFactors(expr, method = "TMM")

### Estimate gene expression values as reads per Kilobase per Million (rpkm).
#crGeneLength <- crGeneLength[keepgenes, ]  # subset the df to the genes that passed the min. cpm and sample count
#all.equal(rownames(normfact$counts), crGeneLength$ensemble_gene_id)  # Confirm the gene IDs are in the same order in both df
#normfact$genes <- data.frame(Length = crGeneLength$Length)  # Assign gene lengths to the normfact object
#rpkm_tmm <- edgeR::rpkm(normfact, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)  # log=TRUE produces log2-transformed normalised counts per million (log2 CPM) 
#rpkm_tmm <- edgeR::rpkm(normfact, gene.length = crGeneLength$Length, normalized.lib.sizes = TRUE, log = TRUE, prior.count = 1)
#boxplot(as.matrix(rpkm_tmm), ylab = expression("Log"[2]~"Read counts"), las = 2, main = "RPKM-TMM")

### Transform the TMM-normalized gene expression values to counts-per-million (cpm) ----------------
# Note that cpm() gives expression values normalized by library size and TMM. https://rdrr.io/bioc/edgeR/src/R/cpm.R
# rpkm() gives expression values normalized by library size, TMM and gene length. https://support.bioconductor.org/p/77193/, https://www.biostars.org/p/317701/
ctmm <- edgeR::cpm(normfact, normalized.lib.sizes = TRUE, 
                   log = FALSE  # FALSE because TRUE produces log2-transformed normalised counts per million (log2 CPM)
                   #log = TRUE, prior.count = 1  # s prior.count = 1 or 2 is required when calculating log2 counts to avoid inf values
)   
#boxplot(as.matrix(ctmm), ylab = expression("Log"[2]~"Read counts"), las = 2, main = "CPM-TMM")
#boxplot(log2(as.matrix(ctmm) + 1), ylab = expression("Log"[2]~"Read counts"), las = 2, main = "CPM-TMM")

# Save CPM-TMM normalized counts in a text file
ctmm_rc <- edgeR::cpm(normfact, normalized.lib.sizes = TRUE, log = FALSE)  # CPM TMM-normalized, no log2 transformed
write.table(ctmm_rc, file = paste0("CPM-TMM-normalized-expression_Herring_RNAseq-BSHBR.",genomeVersion,".",alignerName,"_bulkGenes.txt"), quote = FALSE, row.names = TRUE)


# --------------------------------------------------------------------------------
# Generate a joint plot using ggplot2 graphics
# --------------------------------------------------------------------------------

# Convert the df with normalized read counts to long format for plotting ----------------
library(ggplot2)
library(reshape2)

logCount_tmm <- log2(as.matrix(ctmm) + 1)
logCount_tmm_melt <- reshape2::melt(logCount_tmm) 
names(logCount_tmm_melt)[1:2] <- c("gene_id", "sample")
logCount_tmm_melt$tissue <- rep(NA, nrow(logCount_tmm_melt))
logCount_tmm_melt[grepl(".BR", logCount_tmm_melt$sample),"tissue"] <- "BR"
logCount_tmm_melt[grepl(".BSH", logCount_tmm_melt$sample),"tissue"] <- "BSH"
head(logCount_tmm_melt)
dim(logCount_tmm_melt)
#[1] 296394      4

# Load information of genes of interest (gene symbol and ensembl gene id are required) ----------------
genesOfInterest_df <- read.table(genesOfInterest, header = TRUE, stringsAsFactors = FALSE)
genesOfInterest_df
dim(genesOfInterest_df)
gene_list <- genesOfInterest_df$ensembl_gene_id  # Save the list of gene ids as a separate vector

# Filter dataset to only retain the genes of interest
logCount_tmm_melt_subset <- logCount_tmm_melt[logCount_tmm_melt$gene_id %in% gene_list, ]
head(logCount_tmm_melt_subset)
dim(logCount_tmm_melt_subset)
#[1] 196   4
#[1] 182   4

# Assign the correspondent gene name (symbol) to each gene
logCount_tmm_melt_subset$gene_name <- rep(NA, nrow(logCount_tmm_melt_subset))  # Create a new column for gene names

for(i in seq(1:nrow(genesOfInterest_df))) {  # loop over the file listing the information of genes of interest
  #i=1
  logCount_tmm_melt_subset[grepl(genesOfInterest_df[i,"ensembl_gene_id"], logCount_tmm_melt_subset$gene_id), "gene_name"] <- genesOfInterest_df[i,"gene_name"]
}
head(logCount_tmm_melt_subset)

# Set the column gene_name as a factor and set the gene order wanted (levels).
logCount_tmm_melt_subset$gene_name <- factor(logCount_tmm_melt_subset$gene_name, 
                                             #levels = c("TSHR","DIO2","RHO","OPN3","OPN4","OPN4a","OPN4xa","OPN4xb","OPN6a","OPN6b","OPN7a","OPN7b","OPN7d","OPN8a","OPN8c")
                                             #levels = c("TSHR","RHO","OPN3","OPN4","OPN4a","OPN4xa","OPN4xb","OPN6a","OPN6b","OPN7a","OPN7b","OPN7d","OPN8a","OPN8c")
                                             levels = c("RHO","OPN3","OPN4","OPN4a","OPN4xa","OPN4xb","OPN6a","OPN6b","OPN7a","OPN7b","OPN7d","OPN8a","OPN8c")
                                             )

# Make the plot
library(ggplot2)

p <- ggplot(data = logCount_tmm_melt_subset, aes(x = tissue, y = value, fill = tissue)) +
  geom_boxplot() +
  geom_point(size = 2, pch = 1, stroke = 0.8) +
  facet_grid(. ~ gene_name) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), 
        #legend.position = "none", 
        text = element_text(size = 20),
        strip.text.x = element_text(size = 19, face = "bold.italic")
        ) + 
  scale_fill_manual(values = c("#636363", "#bdbdbd")) +  # c("#fdcdac", "#cbd5e8") "#f4cae4" light pink, #b3e2cd" light green, "#bababa" medium gray
  #ggtitle(paste0(geneName, " (", highlightGene, ")")) +
  labs(y = expression(log[2] ~ ("TMM-normalized counts per million")), 
       x = "", 
       fill = "Tissue type") + 
  scale_y_continuous(limits = c(0,7), breaks = seq(0,7, by = 1))

p

# Save the plot as a PDF file.
pdf(paste0("Comparison_CPM-TMM-normalized-expression_onlyOPSINgenes_Herring_RNAseq_BSHBR_",genomeVersion,"_",alignerName,".pdf"),
    #paste0("Comparison_CPM-TMM-normalized-expression_TSHR-OPSINgenes_Herring_RNAseq_BSHBR_",genomeVersion,"_",alignerName,".pdf"),
    #width = 11, height = 8.5, paper = "a4r"
    width = 20, height = 7
)
print(p)
dev.off()

# Save CPM-TMM normalized counts for only the genes of interest in a separate file
ctmm_subset <- ctmm[rownames(ctmm) %in% gene_list, ]
dim(ctmm_subset)
write.table(ctmm_subset, file = paste0("CPM-TMM-normalized-expression_Herring_RNAseq-BSHBR.",genomeVersion,".",alignerName,"_15genes.txt"), 
            quote = FALSE, row.names = TRUE, col.names = NA, sep = "\t")

sessionInfo()
