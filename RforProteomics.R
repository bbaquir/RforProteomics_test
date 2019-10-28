#test out suggested package for proteomic analysis
#https://bioconductor.org/packages/3.9/data/experiment/html/RforProteomics.html

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

BiocManager::install("RforProteomics", dependencies = TRUE)

library("RforProteomics")

install.packages("mzR", "IPPD", "biomaRt", "MSGFplus")

## gets the vignette source
rfile <- system.file("doc/RforProteomics.R",
                     package = "RforProteomics")
rfile

#prepare the working environment
library("RColorBrewer") ## Color palettes
library("ggplot2")  ## Convenient and nice plotting
library("reshape2") ## Flexibly reshape data

## load the required packages
library("mzR") ## the software package
library("msdata") ## the data package
## below, we extract the releavant example file
## from the local 'msdata' installation
filepath <- system.file("microtofq", package = "msdata")
file <- list.files(filepath, pattern="MM14.mzML",
                   full.names=TRUE, recursive = TRUE)
## creates a commection to the mzML file
mz <- openMSfile(file)
## demonstraction of data access
basename(fileName(mz))

## [1] "MM14.mzML"
runInfo(mz)
## $scanCount
## [1] 112
## 
## $lowMz
## [1] 0
## 
## $highMz
## [1] 0
## 
## $dStartTime
## [1] 270.334
## 
## $dEndTime
## [1] 307.678
## 
## $msLevels
## [1] 1
## 
## $startTimeStamp
## [1] NA
instrumentInfo(mz)
## $manufacturer
## [1] "Unknown"
## 
## $model
## [1] "instrument model"
## 
## $ionisation
## [1] "electrospray ionization"
## 
## $analyzer
## [1] "mass analyzer type"
## 
## $detector
## [1] "detector type"
## 
## $software
## [1] "so_in_0 "
## 
## $sample
## [1] "MM14_20uMsa_0"
## 
## $source
## [1] ""
## once finished, it is good to explicitely
## close the connection
close(mz)

file <- system.file("mzid", "Tandem.mzid.gz", package="msdata")
mzid <- openIDfile(file)
mzid
## Identification file handle.
## Filename:  Tandem.mzid.gz 
## Number of psms:  171

softwareInfo(mzid)
## [1] "xtandem x! tandem CYCLONE (2010.06.01.5) "  
## [2] "ProteoWizard MzIdentML 3.0.501 ProteoWizard"
enzymes(mzid)
##      name nTermGain cTermGain minDistance missedCleavages
## 1 Trypsin         H        OH           0               1
names(psms(mzid))
##  [1] "spectrumID"               "chargeState"             
##  [3] "rank"                     "passThreshold"           
##  [5] "experimentalMassToCharge" "calculatedMassToCharge"  
##  [7] "sequence"                 "modNum"                  
##  [9] "isDecoy"                  "post"                    
## [11] "pre"                      "start"                   
## [13] "end"                      "DatabaseAccess"          
## [15] "DBseqLength"              "DatabaseSeq"             
## [17] "DatabaseDescription"      "spectrum.title"          
## [19] "acquisitionNum"
head(psms(mzid))[, 1:13]
##   spectrumID chargeState rank passThreshold experimentalMassToCharge
## 1   index=12           3    1         FALSE                 903.7209
## 2  index=285           3    1         FALSE                 792.3792
## 3   index=83           3    1         FALSE                 792.5295
## 4   index=21           3    1         FALSE                 850.0782
## 5  index=198           3    1         FALSE                 527.2592
## 6   index=13           2    1         FALSE                 724.8816
##   calculatedMassToCharge                 sequence modNum isDecoy post pre
## 1               903.4032  LCYIALDFDEEMKAAEDSSDIEK      2   FALSE    S   K
## 2               792.3899   KDLYGNVVLSGGTTMYEGIGER      1   FALSE    L   R
## 3               792.3899   KDLYGNVVLSGGTTMYEGIGER      1   FALSE    L   R
## 4               849.7635 VIDENFGLVEGLMTTVHAATGTQK      1   FALSE    V   K
## 5               527.2849          GVGGAIVLVLYDEMK      1   FALSE    R   R
## 6               724.3771            HAVGGRYSSLLCK      1    TRUE    D   K
##   start end
## 1   217 239
## 2   292 313
## 3   292 313
## 4   842 865
## 5   297 311
## 6   392 404

library("mzID")
mzids <- list.files(system.file('extdata', package = 'mzID'),
                    pattern = '*.mzid', full.names = TRUE)
mzids

id <- mzID(mzids[1])
id

ids <- mzID(mzids[1:2])
ids

fid <- flatten(id)
names(fid)
dim(fid)

library("MSnbase")
## uses a simple dummy test included in the package
mzXML <- dir(system.file(package="MSnbase",dir="extdata"),
             full.name=TRUE,
             pattern="mzXML$")
basename(mzXML)

## reads the raw data into and MSnExp instance
raw <- readMSData(mzXML, verbose = FALSE, centroided = TRUE)
raw

## Extract a single spectrum
raw[[3]]

plot(raw, full = TRUE)
plot(raw[[3]], full = TRUE, reporters = iTRAQ4)

###########Quantitiative Proteomics
## Experiment information
library("rpx")
px1 <- PXDataset("PXD000001")
px1

pxfiles(px1)
## Downloading the mzTab data
mztab <- pxget(px1, "PXD000001_mztab.txt")
mztab

## Load mzTab peptide data
qnt <- readMzTabData(mztab, what = "PEP", version = "0.9")

sampleNames(qnt) <- reporterNames(TMT6)
head(exprs(qnt))

## remove missing values
qnt <- filterNA(qnt)
processingData(qnt)

## combine into proteins
## - using the 'accession' feature meta data
## - sum the peptide intensities
protqnt <- combineFeatures(qnt,
                           groupBy = fData(qnt)$accession,
                           fun = sum)

cls <- brewer.pal(5, "Set1")
matplot(t(tail(exprs(protqnt), n = 5)), type = "b",
        lty = 1, col = cls,
        ylab = "Protein intensity (summed peptides)",
        xlab = "TMT reporters")
legend("topright", tail(featureNames(protqnt), n=5),
       lty = 1, bty = "n", cex = .8, col = cls)

qntS <- normalise(qnt, "sum")
qntV <- normalise(qntS, "vsn")
qntV2 <- normalise(qnt, "vsn")

acc <- c("P00489", "P00924",
         "P02769", "P62894",
         "ECA")

idx <- sapply(acc, grep, fData(qnt)$accession)
idx2 <- sapply(idx, head, 3)
small <- qntS[unlist(idx2), ]

idx3 <- sapply(idx, head, 10)
medium <- qntV[unlist(idx3), ]

m <- exprs(medium)
colnames(m) <- c("126", "127", "128",
                 "129", "130", "131")
rownames(m) <- fData(medium)$accession
rownames(m)[grep("CYC", rownames(m))] <- "CYT"
rownames(m)[grep("ENO", rownames(m))] <- "ENO"
rownames(m)[grep("ALB", rownames(m))] <- "BSA"
rownames(m)[grep("PYGM", rownames(m))] <- "PHO"
rownames(m)[grep("ECA", rownames(m))] <- "Background"

cls <- c(brewer.pal(length(unique(rownames(m)))-1, "Set1"),
         "grey")
names(cls) <- unique(rownames(m))
wbcol <- colorRampPalette(c("white", "darkblue"))(256)
heatmap(m, col = wbcol, RowSideColors=cls[rownames(m)])

dfr <- data.frame(exprs(small),
                  Protein = as.character(fData(small)$accession),
                  Feature = featureNames(small),
                  stringsAsFactors = FALSE)
colnames(dfr) <- c("126", "127", "128", "129", "130", "131",
                   "Protein", "Feature")
dfr$Protein[dfr$Protein == "sp|P00924|ENO1_YEAST"] <- "ENO"
dfr$Protein[dfr$Protein == "sp|P62894|CYC_BOVIN"]  <- "CYT"
dfr$Protein[dfr$Protein == "sp|P02769|ALBU_BOVIN"] <- "BSA"
dfr$Protein[dfr$Protein == "sp|P00489|PYGM_RABIT"] <- "PHO"
dfr$Protein[grep("ECA", dfr$Protein)] <- "Background"
dfr2 <- melt(dfr)
## Using Protein, Feature as id variables
ggplot(aes(x = variable, y = value, colour = Protein),
       data = dfr2) +
  geom_point() +
  geom_line(aes(group=as.factor(Feature)), alpha = 0.5) +
  facet_grid(. ~ Protein) + theme(legend.position="none") +
  labs(x = "Reporters", y = "Normalised intensity")

