library(limma)
library(GEOquery)
library(convert)
library(plyr)


#################################### GET GSE && CONVERT GSE Series Matrix files as an ExpressionSet #####


# load series and platform data from GEO
# setwd("C:/Users/Valls-testedDX/Desktop/micrarray_test")
wdpath <- "C:/Users/Valls-testedDX/Desktop/micrarray_test"
a = getGEOSuppFiles("GSE147232")
untar(tarfile = "C:/Users/Valls-testedDX/Desktop/micrarray_test/GSE147232/GSE147232_RAW.tar", exdir = "C:/Users/Valls-testedDX/Desktop/micrarray_test")

# unzip each gz file
unz_files <- list.files(wdpath ,pattern = ".gz" , full.names = TRUE)
ldply(.data= unz_files, .fun = unzip , exdir=wdpath)

###################################### .gpr files processing #############################################

genepixFiles <- list.files(pattern = ".gpr")

RG <- read.maimages(genepixFiles, source="genepix", other.columns=c("F635 SD","B635 SD",
                                                                    "F532 SD","B532 SD","B532 Mean","B635 Mean","F Pixels","B Pixels"))

#readGPRHeader(genepixFiles[1])

#Background Correction for GenePix 
RGmodel <- kooperberg(RG)
MA <- normalizeWithinArrays(RGmodel)

#convert MA to ExpressionSet

eset = as(MA, "ExpressionSet")
new_names <- c("CTR_1" , "CTR_2" , "CTR_3" , "CTR_4",	"CTR_5",	"MCI_1",	"MCI_2",	"MCI_3",	"MCI_4","MCI_5")
eset@phenoData[,1] <- new_names


