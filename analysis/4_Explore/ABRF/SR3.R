library(reshape2)
library(ggplot2) 

raw <- read.csv(file="iPRG_ppm20c_rt6-15-nosingle.csv")

head(raw)

##############################
#### 1. Pre-processing before analysis in MSstats
##############################

##############################
## 1.1 use only modified peptide sequence
raw<-raw[, -2]
colnames(raw) <- c("ProteinName", "PeptideSequence", "PrecursorCharge", "PrecursorMz","FragmentIon", 
                 "ProductCharge","ProductMz", "IsotopeLabelType", "Condition", "BioReplicate", 
                 "Run", "Intensity", "StandardType", "Truncated", "QValue")

## remove decoy protein name
proname <- unique(raw$ProteinName)
length(proname)
decoy <- proname[grep('DECOY', proname)] # 70

raw <- raw[-which(raw$ProteinName %in% decoy), ]
length(unique(raw$ProteinName)) #3027

## check whether spike-in proteins are measured or not
proname[which(proname %in% c("sp|P44015|VAC2_YEAST", "sp|P55752|ISCB_YEAST", 
                                "sp|P44374|SFG2_YEAST", "sp|P44983|UTR6_YEAST", 
                                "sp|P44683|PGA4_YEAST", "sp|P55249|ZRT4_YEAST"))]


## class of intensity is factor, change it as numeric
class(raw$Intensity)
raw$Intensity <- as.numeric(as.character(raw$Intensity))

## zero value is available
range(raw[!is.na(raw$Intensity), "Intensity"])

## there are NAs and zeros in skyline output
sum(is.na(raw$Intensity)) # 231
sum(!is.na(raw$Intensity) & raw$Intensity == 0) # 37279


##############################
## 1.2 remove truncated peaks
sum(raw$Truncated == "True") ## 795 truncated peaks
raw[raw$Truncated == "True", "Intensity"] <- NA ## replace with NA
sum(is.na(raw$Intensity)) # now 1026 NAs


##############################
## 1.3 Sum for isotopic peaks per peptide and precursor charge

## add the column for unique peptide and precursor
raw$pepprecursor <- paste(raw$PeptideSequence, raw$PrecursorCharge, sep="_")

raw <- raw[!is.na(raw$Intensity), ]

## sum of mooisotopic peaks
data_w <- dcast( pepprecursor ~ Run, data=raw, value.var='Intensity', fun.aggregate=function(x) sum(x, na.rm=TRUE), fill=NA_real_) 
head(data_w)
dim(data_w)

## make long format
newdata <- melt(data_w, id.vars=c('pepprecursor'))
colnames(newdata)[colnames(newdata) %in% c("variable","value")] <- c('Run','Intensity')

## assignn protein name
uniinfo <- unique(raw[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "pepprecursor")])		
raw <- merge(newdata, uniinfo, by="pepprecursor")


##############################
## 1.4 assign the annotation

## read annotation information
annot <- data.frame(Run=levels(raw$Run), 
                    Condition = c('Condition1', 'Condition1', 'Condition1', 
                                  'Condition2', 'Condition2', 'Condition2', 
                                  'Condition3', 'Condition3', 'Condition3', 
                                  'Condition4', 'Condition4', 'Condition4'),
                    BioReplicate = c('1', '1', '1', '2', '2', '2', '3', '3', '3', '4', '4', '4'))

## merge it by Run
raw <- merge(raw, annot, by="Run")
head(raw)

## add other required information
raw$FragmentIon <- "sum"
raw$ProductCharge <- NA
raw$IsotopeLabelType <- "L"

unique(raw[, c("Run", "BioReplicate", "Condition")])

colnames(raw)

raw <- raw[, c(4,5,6,9,10,11,7,8,1,3)]

##############################
## preprocessing is done.


##############################
#### 2. Analysis in MSstats
##############################

library(MSstats)

##############################
## 2.1 Transformation = log2
##     Normalization = equalize median
##     Model-based run-level summarization after imputation

quant <- dataProcess(raw,
                     normalization='equalizeMedians',
                     summaryMethod="TMP",
                     cutoffCensored="minFeature",
                     censoredInt="0",
                     MBimpute=TRUE, 
                     skylineReport=TRUE)


##############################
## 2.2 Data visualization

dataProcessPlots(quant, type="QCplot", ylimDown=0, width=7, height=7, which.Protein = 1, address="Skyline_intensity_new_")

dataProcessPlots(quant, type="Profileplot", ylimDown=0, featureName="NA", width=7, height=7, address="Skyline_intensity_new_")

dataProcessPlots(quant, type="Conditionplot", ylimDown=0, featureName="NA", width=7, height=7, address="Skyline_intensity_new_")



##############################
## 2.3 Model-based comparison and adjustment for multiple testing

comparison1<-matrix(c(-1,1,0,0),nrow=1)
comparison2<-matrix(c(-1,0,1,0),nrow=1)
comparison3<-matrix(c(-1,0,0,1),nrow=1)
comparison4<-matrix(c(0,-1,1,0),nrow=1)
comparison5<-matrix(c(0,-1,0,1),nrow=1)
comparison6<-matrix(c(0,0,-1,1),nrow=1)
comparison<-rbind(comparison1, comparison2, comparison3, comparison4, comparison5, comparison6)
row.names(comparison)<-c("C2-C1","C3-C1","C4-C1","C3-C2","C4-C2","C4-C3")


test <- groupComparison(contrast.matrix=comparison, data=quant)
Skyline.intensity.comparison.result <- test$ComparisonResult

write.csv(Skyline.intensity.comparison.result, file='Supplementary Table ST4.csv')

##############################
## 2.4 Visualization for testing result

groupComparisonPlots(Skyline.intensity.comparison.result, type="VolcanoPlot", address="Skyline_intensity_new_")

groupComparisonPlots(Skyline.intensity.comparison.result, type="Heatmap", address="Skyline_intensity_new_")

groupComparisonPlots(Skyline.intensity.comparison.result, type="ComparisonPlot", address="Skyline_intensity_new_")

