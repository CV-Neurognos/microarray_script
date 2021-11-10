#General Bioconductor packages
library(Biobase)
library(oligoClasses)

#Annotation and data import packages
library(ArrayExpress)
library(GEOquery)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)

#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)

#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

#Plotting and color options packages
library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)

#Formatting/documentation packages
#library(rmarkdown)
#library(BiocStyle)
library(dplyr)
library(tidyr)

#Helpers:
library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)
#library(devtools)

######################## DEFINE CONSTANTS AND DIRECTORIES #################

GSE_QUERY <- "GSE157239"
working_directory_path <- getwd()
OUTPUT_QC <- paste (working_directory_path,"QC", sep = "/")
OUTPUT_QC_2 <- paste (working_directory_path,"QC_2", sep = "/")
OUTPUT_RESULT <- paste (working_directory_path,"RESULT", sep = "/")
raw_data_dir <- tempdir()

if (!dir.exists(raw_data_dir)) {
  dir.create(raw_data_dir)
}

if (!dir.exists(OUTPUT_RESULT)) {
  dir.create(OUTPUT_RESULT)
}

################### DOWNLOAD FILES #########################################
#Array_express
#anno_AE <- getAE("E-MTAB-2967", path = raw_data_dir, type = "raw")

#Genome Expression Omnibus
### getGEO has a bug , corrected by readr::local_edition(1)
readr::local_edition(1)
gset <- getGEO(GSE_QUERY, GSEMatrix =TRUE, AnnotGPL=FALSE)
#if (length(gset) > 1) idx <- grep("GPL18058", attr(gset, "names")) else idx <- 1
gset <- gset[[1]]


# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))


###################   THIS PIPELINE IS FOR AFFYMETRIX DATASETS #############

# add Sample and Data Relationship Format (SDRF)

sdrf_location <- file.path(raw_data_dir, "E-MTAB-2967.sdrf.txt")
SDRF <- read.delim(sdrf_location)

rownames(SDRF) <- SDRF$Array.Data.File
SDRF <- AnnotatedDataFrame(SDRF)


raw_data <- oligo::read.celfiles(filenames = file.path(raw_data_dir,
                                                       SDRF$Array.Data.File),
                                 verbose = FALSE, phenoData = SDRF)
stopifnot(validObject(raw_data))


Biobase::pData(raw_data) <- Biobase::pData(raw_data)[, c("Source.Name",
                                                         "Characteristics.individual.",
                                                         "Factor.Value.disease.",
                                                         "Factor.Value.phenotype.")]

# Quality control of the raw data
#First: transform expression set into log2 
exp_raw <- log2(Biobase::exprs(raw_data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)
percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])
dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     Disease = pData(raw_data)$Factor.Value.disease.,
                     Phenotype = pData(raw_data)$Factor.Value.phenotype.,
                     Individual = pData(raw_data)$Characteristics.individual.)
#Second: graphic PCA

ggplot(dataGG, aes(PC1, PC2)) +
  geom_point(aes(shape = Disease, colour = Phenotype)) +
  ggtitle("PCA plot of the log-transformed raw expression data") +
  xlab(paste0("PC1, VarExp: ", percentVar[1], "%")) +
  ylab(paste0("PC2, VarExp: ", percentVar[2], "%")) +
  theme(plot.title = element_text(hjust = 0.5))+
  coord_fixed(ratio = sd_ratio) +
  scale_shape_manual(values = c(4,15)) +
  scale_color_manual(values = c("darkorange2", "dodgerblue4"))

#Alternative: BOXplot directly from expressionSet                  
oligo::boxplot(raw_data, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

# Another QC assessment: arrayQualityMetrics() -- EXCELLENT: WORKS on ExpressionSets

arrayQualityMetrics(expressionset = raw_data,
                    outdir = OUTPUT_QC,
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))



## Background adjustment, calibration, summarization and annotation
## Deconvolution for background correction, quantile normalization and summarization with RMA (robust multichip average)

palmieri_eset_norm <- oligo::rma(raw_data, target = "core")



## post-normalization QC
exp_palmieri <- Biobase::exprs(palmieri_eset_norm)
arrayQualityMetrics(expressionset = raw_data,
                    outdir = OUTPUT_QC_2,
                    force = TRUE, do.logtransform = TRUE,
                    intgroup = c("Factor.Value.disease.", "Factor.Value.phenotype."))


#Filtering low intensity wells // low hibridized probes

palmieri_medians <- rowMedians(Biobase::exprs(palmieri_eset_norm))

hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

man_threshold <- 4

hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

abline(v = man_threshold, col = "coral4", lwd = 2)


no_of_samples <-
  table(paste0(pData(palmieri_eset_norm)$Factor.Value.disease., "_",
               pData(palmieri_eset_norm)$Factor.Value.phenotype.))
idx_man_threshold <- apply(Biobase::exprs(palmieri_eset_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)
samples_cutoff <- min(no_of_samples)
palmieri_manfiltered <- subset(palmieri_eset_norm, idx_man_threshold)

#### Annotation of the transcript clusters

anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(palmieri_manfiltered)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

### Remove multiple mappings  

anno_grouped <- group_by(anno_palmieri, PROBEID)
anno_summarized <-
  dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))

anno_filtered <- filter(anno_summarized, no_of_matches > 1)

probe_stats <- anno_filtered

ids_to_exlude <- (featureNames(palmieri_manfiltered) %in% probe_stats$PROBEID)

palmieri_final <- subset(palmieri_manfiltered, !ids_to_exlude)

fData(palmieri_final)$PROBEID <- rownames(fData(palmieri_final))

fData(palmieri_final) <- left_join(fData(palmieri_final), anno_palmieri)

rownames(fData(palmieri_final)) <- fData(palmieri_final)$PROBEID


######################            A linear model for the data           ######################




###################   THIS PIPELINE IS FOR GENERAL USE #############
dim(gset)
head(pData(gset))

#hibridization distribution

exprs_row_meadians <- rowMedians(exprs(gset))

hist_res <- hist(exprs_row_meadians, 50, col = "cornsilk1", freq = TRUE,
                 main = "Histogram of the median intensities",
                 border = "antiquewhite4",
                 xlab = "Median intensities")

#pre-Normalization QC

#Alternative: BOXplot directly from expressionSet
exp_raw <- log2(exprs(gset))

oligo::boxplot(exp_raw, target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

#Normalization

exprs(gset) <- normalizeBetweenArrays(exp_raw)


#post-Normalization QC

oligo::boxplot(exprs(gset), target = "core",
               main = "Boxplot of log2-intensitites for the raw data")

# Change names to construct the model.matrix, For each microarray dataset tune this part.

pData(gset)[,"title"]<- gsub("control","HD", pData(gset)[,"title"])
pData(gset)[,"diagnosis:ch1"] <- c("HD", "HD" ,"HD" ,"AD" ,"HD", "AD" ,  "HD" , "AD" , "AD" ,"HD", "AD" ,"HD" , "AD" ,"AD","AD" , "HD")

# assign samples to groups and set up design matrix
gs <- factor(pData(gset)[,"diagnosis:ch1"] )
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

gset <- gset[complete.cases(exprs(gset)), ]

v <- vooma(gset, design, plot=T)

# OR weights by group
# v <- voomaByGroup(gset, group=groups, design, plot=T, cex=0.1, pch=".", col=1:nlevels(gs))

v$genes <- fData(gset) # attach gene annotations

# fit linear model
fit  <- lmFit(v)

# set up contrasts of interest and recalculate model coefficients
group_contrast = unique(pData(gset)[,"diagnosis:ch1"], incomparables = FALSE)
cts <- c(paste(group_contrast[1],"-",group_contrast[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","miRNA_ID","SPOT_ID"))


# Visualize and quality control test results.
# Build histogram of P-values for all genes. Normal test
# assumption is that most genes are not differentially expressed.

tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

output_dir_file = paste(paste(OUTPUT_RESULT,GSE_QUERY, sep="/"),".csv",sep="")

write.csv(tT2, output_dir_file , row.names = FALSE)

