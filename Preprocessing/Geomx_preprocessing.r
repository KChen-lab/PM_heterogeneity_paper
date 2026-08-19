#load the libraries
library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)
library(standR)
library(SummarizedExperiment)
library(EDASeq)
library(RUVSeq)
library(DESeq2)
library(edgeR)
library(limma)

# ====Load the data====

#set file path to location of parent folder containing dcc, pkc, and annotation files
datadir <- file.path("/path/to/files")

#automatically list files in each directory for use
DCCFiles <- dir(file.path(datadir, "dcc"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "pkc"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)
sampleAnnotFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
                       full.names = TRUE, recursive = TRUE)

#load data to create data object
mesoData <- 
  readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = sampleAnnotFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))
						 
#remove metastatic sample from the data set so it doesn't skew data
mesoData <- mesoData[, pData(mesoData)$PrimaryTumor_Vs_Metastatic != "Metastatic"]

#extract necessary info to make into summarized experiment for standr visualizations
counts <- exprs(mesoData)

#calculate library size
lib_size <- colSums(counts)

#add library size to phenodata
pData(mesoData)$lib_size <- lib_size

#shift all counts to 1 to enable downstream transformations
mesoData <- shiftCountsOne(mesoData, useDALogic = TRUE)

#Segment QC
#default parameters found in (), study specific ones chosen after QC visualization
QC_params <- list(minSegmentReads = 1000,  #minimum number of reads (1000)
                  percentTrimmed = 80,     #minimum % of reads trimmed (80%)
                  percentStitched = 80,    #minimum % of reads stitched (80%)
                  percentAligned = 80,     #minimum % of reads aligned (80%)
                  percentSaturation = 50,  #minimum sequencing saturation (50%)
                  minNegativeCount = 1,    #minimum negative control counts (10)
                  maxNTCCount = 9000,      #maximum counts observed in NTC well (1000)
                  minNuclei = 100,         #minimum # of nuclei estimated (100),
                  minArea = 5000)          #minimum segment area (5000)

mesoData <- setSegmentQCFlags(mesoData, qcCutoffs = QC_params)

#ensure expected pkcs have been loaded
pkcs <- annotation(mesoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

#collate all QC results
QCResults <- protocolData(mesoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <- 
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

##remove flagged segments
mesoData <- mesoData[, QCResults$QCStatus == "PASS"]

##Probe QC
mesoData <- setBioProbeQCFlags(mesoData,
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20),
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(mesoData)[["QCFlags"]]

#define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))

#exclude outliers
#subset object to exclude all that did not pass ratio & global testing
#there will always be at least one probe for every gene target, this doesn't remove genes from your data
ProbeQCPassed <- subset(mesoData,
                        fData(mesoData)[["QCFlags"]][, c("LowProbeRatio")] == FALSE &
                          fData(mesoData)[["QCFlags"]][, c("GlobalGrubbsOutlier")] == FALSE)

#check dimensions of probeQC
dim(ProbeQCPassed)
mesoData <- ProbeQCPassed

#create gene-level count data
#check how many unique targets the object has
length(unique(featureData(mesoData) [["TargetName"]]))

#collapse to targets
target_mesoData <- aggregateCounts(mesoData)
dim(target_mesoData)
exprs(target_mesoData)[1:5, 1:2]

## limit of Quantification

#calculated based on distribution of negative control probes
#intended to approximate quantifiable limit of gene expression per segment
#threshold is typically 2 geometric std above the geometric mean

#define LOQ std threshold and minimum value
cutoff <- 2
minLOQ <- 2

#calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_mesoData))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if(all(vars[1:2] %in% colnames(pData(target_mesoData)))) {
    LOQ[, module] <- pmax(minLOQ,
                          pData(target_mesoData)[, vars[1]] *
                            pData(target_mesoData)[, vars[2]] ^ cutoff)
  }
}

pData(target_mesoData)$LOQ <- LOQ

##filtering

#filter out segments/genes with abnormally low signal
#important for focusing on true biologically meaningful data

#determine the number of genes detected in each segment
LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_mesoData)$Module == module
  Mat_i <- t(esApply(target_mesoData[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

#ensure ordering since this is stored outside of the geomxset
LOQ_Mat <- LOQ_Mat[fData(target_mesoData)$TargetName, ]

## segment gene detection

#save detection rate information to pheno data
pData(target_mesoData)$GenesDetected <-
  colSums(LOQ_Mat, na.rm = TRUE)
pData(target_mesoData)$GeneDetectionRate <- 
  pData(target_mesoData)$GenesDetected / nrow(target_mesoData)

#determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_mesoData)$DetectionThreshold <-
  cut(pData(target_mesoData)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

#remove samples with less than 10% of genes detected
target_mesoData <- target_mesoData[, pData(target_mesoData)$GeneDetectionRate >= .1]
dim(target_mesoData)

##Gene Detection Rates

#load libraries
library(scales)

#calculate detection rate
LOQ_Mat <- LOQ_Mat[, colnames(target_mesoData)]
fData(target_mesoData)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_mesoData)$DetectionRate <- 
  fData(target_mesoData)$DetectedSegments / nrow(pData(target_mesoData))

#define gene of interest detection table
goi <- c("PDCD1", "CD274", "CD8A", "CD68", "EPCAM", 
         "KRT18", "NPHS1", "NPHS2", "CALB1", "CLDN8")
goi_df <- data.frame(Gene = goi,
                     Number = fData(target_mesoData)[goi, "DetectedSegments"],
                     DetectionRate = percent(fData(target_mesoData)[goi, "DetectionRate"]))

#subset to target genes detected in at least 10% of samples
#manually include negative control probe for downstream usage
negativeProbefData <- subset(fData(target_mesoData), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_mesoData <-
  target_mesoData[fData(target_mesoData)$DetectionRate >= 0.1 | 
                    fData(target_mesoData)$TargetName %in% neg_probes, ]
dim(target_mesoData)

#retain only detected genes of interest
goi <- goi[goi %in% rownames(target_mesoData)]

#get housekeeping genes from the nCounter panel
housekeeping <- c("FCF1","DHX16","ABCF1","EDC3","HDAC3",
  "CNOT10","SAP130","AGK","AMMECR1L","SF3A3","COG7",
  "TMUB2","ZC3H14","DDX50","MRPS5","ZNF346","MTMR14","ERCC3",
  "EIF2B4","TLK2","NUBP1","USP39","NOL7","CNOT4","ZNF143","PRPF38A")

#set counts variable to the expression matrix from the target_mesoData
counts <- exprs(target_mesoData)

#set metadata to the pData from target_mesoData
metadata <- pData(target_mesoData)

#double check: all housekeeping genes are in rownames
housekeeping <- housekeeping[housekeeping %in% rownames(counts)]
stopifnot(length(housekeeping) > 0)

#create a new set
set <- newSeqExpressionSet(as.matrix(round(counts)),
                           phenoData = AnnotatedDataFrame(metadata))

#run RUVg set experiment
set_ruv <- RUVg(set, housekeeping, k = 1)

metadata <- pData(target_mesoData)
pheno <- pData(set_ruv[3])

#add RUVg covariates to limma design
W_1 <- pheno$W_1
segments <- metadata$segment
repeats <- metadata$MRN
histology <- metadata$Histology
sex <- metadata$Sex_F_M
age <- metadata$Age_med_67
smoking <- metadata$Smoking_Y_N
slide_id <- metadata$"slide name"
asbestos <- metadata$Exposure_to_asbestos_Y_N
hereditary <- metadata$Hereditary_mesothelioma_Y_N
stage <- metadata$Stage_AJCC8
recurrence <- metadata$Early_Recurrence_Y_N
neochem <- metadata$NeoChem_Y_N

#create differentially expressed genes list 
y <- DGEList(counts = counts, group = segments)

#include others in factor model
table(y$samples$group)
y$samples$segments <- segments
y$samples$W_1 <- W_1
y$samples$repeats <- repeats
y$samples$histology <- histology
y$samples$sex <- sex
y$samples$age <- age
y$samples$smoking <- smoking
y$samples$slide_id <- slide_id
y$samples$asbestos <- asbestos
y$samples$hereditary <- hereditary
y$samples$stage <- stage
y$samples$recurrence <- recurrence
y$samples$neochem <- neochem

#filter and convert to logCPM
keep <- filterByExpr(y)
summary(keep)
y <- y[keep, ,keep.lib.sizes = FALSE]

#apply trimmed mean of M-values (TMM) normalization
#calculate scaling factors to convert raw library sizes into effective library sizes
y <- calcNormFactors(y, method = "TMM")

#create design matrix
design <- model.matrix(~0 + W_1 + segments)
colnames(design)
fit <- voomLmFit(y, design, block = y$samples$repeats, sample.weights = TRUE, plot = TRUE)

#remove batch effects associated with multiple slides from the same patient
fit$EList$E_2 <- limma::removeBatchEffect(fit$EList$E, batch = metadata$"slide name")
normalized_data <- fit$EList$E_2
