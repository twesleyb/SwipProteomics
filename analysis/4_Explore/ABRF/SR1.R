
##############################
## Analysis for peak intensity data
## input : SkylineIntensity_unique.tsv, which was provided to participants
##############################

##############################
#### 1. Pre-processing before analysis in MSstats
##############################

library(reshape2)

##############################
## 1.1 Read data and modify the column name and class

raw <- read.csv(file="SkylineIntensities_unique.tsv", sep="\t")

head(raw)

colnames(raw) <- c("PeptideSequence", "ProteinName", "Run", "Precursor.Mz", "PrecursorCharge", "Product.Mz",
                   "ProductCharge", "FragmentIon", "Retention.Time", "Intensity", "Probability", "QValue")

## class of intensity is factor, change it as numeric
class(raw$Intensity)
raw$Intensity <- as.numeric(as.character(raw$Intensity))

##############################
## 1.2 remove duplicated rows

##############################
### 1.2.1 check unique rows : remove peptide sequences which are not unique for protein
tmp <- raw
unique(tmp$FragmentIon)

## new column for unique feature (peptide  * charge)
tmp$fea <- paste(tmp$PeptideSequence, tmp$PrecursorCharge, sep="_")

## there are 12 MS runs
length(unique(tmp$Run))

## there are 29575 unique featuares
length(unique(tmp$fea))

## therefore 1064700 rows should be in the data (3 monoisotopic peaks per feature)
length(unique(tmp$fea)) * length(unique(tmp$Run)) * 3

## however there are more rows, which means there are duplicated rows
dim(tmp)


##############
### 1.2.2 multiple Precursor.Mz

## check unique rows
tmptmp <- xtabs(~fea, tmp)
unique(tmptmp) # there are some features which has more than one measurement
unique(tmptmp)/36

## get feature id which has more than one measurement per run
featureid <- names(tmptmp[tmptmp > 36])
featureid ## 1086 features

## look details for example
tmptmp[names(tmptmp) == "YYPQQAPMPAAAPQQAYYGTAPSTSK_3"]
temptemp <- tmp[tmp$fea == "YYPQQAPMPAAAPQQAYYGTAPSTSK_3" & tmp$Run == '1A', ]
dim(temptemp)
temptemp

## In this case, there are total 72 meausurements for one peptide and charge. (it should be 36 measurements)
## Because there are two different Precursor.Mz values.
## need to remove one of Precursor.MZ information

unid <- unique(tmp[, c("fea","Precursor.Mz")])
dim(unid)
head(unid)
unitemp <- xtabs(~ fea, unid)
unique(unitemp) # 1,2,3,4 precursor.mz per feature are available

## get feature id which has more than one measurement per run
featureid <- names(unitemp[unitemp > 1])
featureid ## 1084 features


## save the feature and precursor.mz that is needed to keep
selectfeature <- NULL

## get the combination of feature and precursor.mz which need to be kept

for(i in 1:length(featureid)){

  sub <- tmp[tmp$fea == featureid[i], ]
  sub$fea.precursor.mz <- paste(sub$fea, sub$Precursor.Mz, sep="_")

  # 1. get any m/z which have less NA Qvalue,
  sub <- sub[!is.na(sub$QValue), ]
  count <- xtabs(~ fea.precursor.mz, sub)

  if ( any(count[-1] != count[1]) ) {
    selectfeature <- c(selectfeature, names(count[which.max(count)]))
  } else {

    # 2. if the number of NA Q values are the same, use precursor.mz with more Q values less than 0.10
    sub <- sub[sub$QValue < 0.1, ]
    count <- xtabs(~ fea.precursor.mz, sub)

    if ( any(count[-1] != count[1]) ) {
      selectfeature <- c(selectfeature, names(count[which.max(count)]))
    } else {

      # 3. if the number of NA and less than 0.10 Q values are the same, use precursor.mz with highest mean intensity.
      meanfeature <- aggregate(Intensity ~ fea.precursor.mz, data=sub, function(x) median(x, na.rm=TRUE))
      meanfeature <- meanfeature[order(meanfeature$Intensity, decreasing=T), ]

      ## choose top n
      maxfeature <- meanfeature$fea.precursor.mz[1]

      selectfeature <- c(selectfeature, as.character(maxfeature))
    }
  }
  print(i)
}


raw$fea <- paste(raw$PeptideSequence, raw$PrecursorCharge, sep="_")

## keep features that has no problem
raw1 <- raw[-which(raw$fea %in% featureid), ]
dim(raw1)

## get subset which has features that have more than one precursor.mz per feature
raw2 <- raw[which(raw$fea %in% featureid), ]
dim(raw2)

## keep one precursor.mz per feature
raw2$fea.precursor.mz <- paste(raw2$fea, raw2$Precursor.Mz, sep="_")

raw2 <- raw2[which(raw2$fea.precursor.mz %in% selectfeature), ]
dim(raw2)
raw2 <- raw2[, -which(colnames(raw2) %in% "fea.precursor.mz")]

raw <- rbind(raw1, raw2)
dim(raw)


##############
### 1.2.3 multiple measurements

## check unique rows again
tmptmp <- xtabs(~fea, raw)
unique(tmptmp) # Still there are some features which have multiple measurements

## get feature id which has more than one measurement per run
featureid <- names(tmptmp[tmptmp > 36])
featureid ## 60 features

## look details for some example
tmptmp[names(tmptmp) == "APQAAGVIHNDLMNTFILAQVMK_3"]
temptemp <- raw[raw$fea == "APQAAGVIHNDLMNTFILAQVMK_3", ]
dim(temptemp)
temptemp

## suspect two groups of retention.Time, however, not always.
## get highest abundance if there are multiple measurements

## keep features that has no problem
raw1 <- raw[-which(raw$fea %in% featureid), ]
dim(raw1)

## get subset which has features that have more than one precursor.mz per feature
raw2 <- raw[which(raw$fea %in% featureid), ]
dim(raw2)


## new data frame including highest measurements per feature and per run
raw2new <- NULL

## loop for 60 features which have multiple measurements per run
for(i in 1:length(featureid)){

  sub <- raw2[raw2$fea == featureid[i], ]

  ## sum of mooisotopic peaks
  data_w <- dcast( PeptideSequence + ProteinName + Precursor.Mz + PrecursorCharge + Product.Mz
                   + ProductCharge + FragmentIon ~ Run, data=sub, value.var='Intensity',
                   fun.aggregate=function(x) max(x, na.rm=TRUE), fill=NULL)

  ## make long format
  newdata <- melt(data_w, id.vars=c('PeptideSequence', 'ProteinName', 'Precursor.Mz', 'PrecursorCharge',
                                    'Product.Mz', 'ProductCharge', 'FragmentIon'))
  colnames(newdata)[colnames(newdata) %in% c("variable","value")] <- c('Run','Intensity')

  raw2new <- rbind(raw2new, newdata)

  print(i)
}

## add missing infomation and rearrange data frame
raw2info <- unique(raw2[, c("fea", "Run", "Probability", 'QValue')])

raw2new$fea <- paste(raw2new$PeptideSequence, raw2new$PrecursorCharge, sep="_")

raw2new <- merge(raw2new, raw2info, by=c('fea', 'Run'))
raw2new <- raw2new[, c(4, 3, 5, 6, 7, 8, 9, 2, 10, 11, 12, 1)]

## merge raw1 and new raw2
raw1 <- raw1[, -which(colnames(raw1) %in% c("Retention.Time"))]
raw1 <- raw1[, c(2, 1, 4, 5, 6, 7, 8, 3, 9, 10, 11, 12)]

raw <- rbind(raw1, raw2new)
dim(raw)


##############################
## 1.3 Sum for monoisotopic peaks per peptide and precursor charge

## add the column for unique peptide and precursor
raw$pepprecursor <- paste(raw$PeptideSequence, raw$PrecursorCharge, sep="_")
dim(raw)

## sum of mooisotopic peaks
data_w <- dcast( pepprecursor ~ Run, data=raw, value.var='Intensity', fun.aggregate=sum, fill=NULL)

## make long format
newdata <- melt(data_w, id.vars=c('pepprecursor'))
colnames(newdata)[colnames(newdata) %in% c("variable","value")] <- c('Run','Intensity')

## assignn protein name
uniinfo <- unique(raw[, c("ProteinName", "PeptideSequence", "PrecursorCharge", "pepprecursor")])
raw <- merge(newdata, uniinfo, by="pepprecursor")
dim(raw)
head(raw)


##############################
## 1.4 assign the annotation

## read annotation information
annot <- data.frame(Run=unique(raw$Run),
                    Condition = c('Condition1', 'Condition2', 'Condition3',
                                  'Condition2', 'Condition3', 'Condition4',
                                  'Condition3', 'Condition4', 'Condition2',
                                  'Condition4', 'Condition1', 'Condition1'),
                    BioReplicate = c('1', '2', '3', '2', '3', '4', '3', '4', '2', '4', '1', '1'))

## merge it by Run
raw <- merge(raw, annot, by="Run")
head(raw)

## add other required information
raw$FragmentIon <- "sum"
raw$ProductCharge <- NA
raw$IsotopeLabelType <- "L"

unique(raw[, c("Run", "BioReplicate", "Condition")])

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

dataProcessPlots(quant, type="QCplot", ylimDown=0, width=7, height=7, which.Protein = 1, address="Skyline_intensity_")

dataProcessPlots(quant, type="Profileplot", ylimDown=0, featureName="NA", width=7, height=7, address="Skyline_intensity_")

dataProcessPlots(quant, type="Conditionplot", ylimDown=0, featureName="NA", width=7, height=7, address="Skyline_intensity_")



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

write.csv(Skyline.intensity.comparison.result, file='Supplementary Table ST1.csv')

##############################
## 2.4 Visualization for testing result

groupComparisonPlots(Skyline.intensity.comparison.result, type="VolcanoPlot", address="Skyline_intensity_")

groupComparisonPlots(Skyline.intensity.comparison.result, type="Heatmap", address="Skyline_intensity_")

groupComparisonPlots(Skyline.intensity.comparison.result, type="ComparisonPlot", address="Skyline_intensity_")
