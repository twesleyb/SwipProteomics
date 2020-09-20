
##############################
## Analysis for spectral count data
##############################

##############################
## 1. Read data and reformat
## from ftp > distro > id > tsv
## input : ID_1A.tsv, ID_1B.tsv, ID_1C.tsv, ID_2A.tsv, ID_2B.tsv, ID_2C.tsv,  ID_3A.tsv, ID_3B.tsv, ID_3C.tsv, ID_4A.tsv, ID_4B.tsv, ID_4C.tsv
##############################

r1 <- read.csv(file="ID_1A.tsv", sep="\t")
r1$Run <- '1A'
r1$Condition <- 'Condition1'
head(r1)
dim(r1) # 28862
length(unique(r1$Protein)) ## 3998

r2 <- read.csv(file="ID_1B.tsv", sep="\t")
r2$Run <- '1B'
r2$Condition <- 'Condition1'
head(r2)
dim(r2) # 29656
length(unique(r2$Protein)) ## 4055

r3 <- read.csv(file="ID_1C.tsv", sep="\t")
r3$Run <- '1C'
r3$Condition <- 'Condition1'
head(r3)
dim(r3) # 28238
length(unique(r3$Protein)) ## 4074

r4 <- read.csv(file="ID_2A.tsv", sep="\t")
r4$Run <- '2A'
r4$Condition <- 'Condition2'
head(r4)
dim(r4) # 27698
length(unique(r4$Protein)) ## 4013

r5 <- read.csv(file="ID_2B.tsv", sep="\t")
r5$Run <- '2B'
r5$Condition <- 'Condition2'
head(r5)
dim(r5) # 29436
length(unique(r5$Protein)) ## 4041

r6 <- read.csv(file="ID_2C.tsv", sep="\t")
r6$Run <- '2C'
r6$Condition <- 'Condition2'
head(r6)
dim(r6) # 28495
length(unique(r6$Protein)) ## 4039

r7 <- read.csv(file="ID_3A.tsv", sep="\t")
r7$Run <- '3A'
r7$Condition <- 'Condition3'
head(r7)
dim(r7) # 28485
length(unique(r7$Protein)) ## 4042

r8 <- read.csv(file="ID_3B.tsv", sep="\t")
r8$Run <- '3B'
r8$Condition <- 'Condition3'
head(r8)
dim(r8) # 29620
length(unique(r8$Protein)) ## 4068

r9 <- read.csv(file="ID_3C.tsv", sep="\t")
r9$Run <- '3C'
r9$Condition <- 'Condition3'
head(r9)
dim(r9) # 29128
length(unique(r9$Protein)) ## 4031

r10 <- read.csv(file="ID_4A.tsv", sep="\t")
r10$Run <- '4A'
r10$Condition <- 'Condition4'
head(r10)
dim(r10) # 27368
length(unique(r10$Protein)) ## 3989

r11 <- read.csv(file="ID_4B.tsv", sep="\t")
r11$Run <- '4B'
r11$Condition <- 'Condition4'
head(r11)
dim(r11) # 28776
length(unique(r11$Protein)) ## 3967

r12 <- read.csv(file="ID_4C.tsv", sep="\t")
r12$Run <- '4C'
r12$Condition <- 'Condition4'
head(r12)
dim(r12) # 28073
length(unique(r12$Protein)) ## 4026


raw <- rbind(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12)
head(raw)

unique(raw[, c('Run', 'Condition')])
length(unique(raw$Protein)) # 5776

pepid.Data <- raw

head(pepid.Data)

## Reformat
library(reshape2)

Y <- dcast(Protein ~ Run, data=pepid.Data)

countData <- as.matrix(Y[,-c(1)])
protName <- as.character(Y[,1])      


##############################
## 2. DESeq2: NB GLM + Wald test
##############################

library(DESeq2)

## Format the data for DESeq
## use M1 as reference in the design matrix
colData <- data.frame(M1=factor(c(rep(1,3),rep(0, 9))), 
                      M2=factor(c(rep(0, 3), rep(1,3), rep(0, 6))),
                      M3=factor(c(rep(0, 6), rep(1,3),rep(0, 3))), 
                      M4=factor(c(rep(0, 9), rep(1,3))))

## Comparisons with 1st one  

## Design model ~ M2+M3+M4
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ M2+M3+M4)

## Test for differential expression with NB distribution and exact test
dds <- DESeq(dds)
dds # object of DESeqDataSet class 
resultsNames(dds)

## block automatic independent filtering after testing
c1vs2 <- results(dds, contrast=c('M2', '0', '1'), independentFiltering = FALSE)
c1vs3 <- results(dds, contrast=c('M3', '0', '1'), independentFiltering = FALSE)
c1vs4 <- results(dds, contrast=c('M4', '0', '1'), independentFiltering = FALSE)

rm(dds)

## Comparisons with 2nd one 

## Design model ~ M1+M3+M4
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ M1+M3+M4)

## Test for differential expression with NB distribution and exact test
dds <- DESeq(dds)

c2vs3 <- results(dds, contrast=c('M3', '0', '1'), independentFiltering = FALSE)
c2vs4 <- results(dds, contrast=c('M4', '0', '1'), independentFiltering = FALSE)

rm(dds)

## Comparisons with 3rd one 

## Design model ~ M1+M2+M4
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ M1+M2+M4)

## Test for differential expression with NB distribution and exact test
dds <- DESeq(dds)

c3vs4 <- results(dds, contrast=c('M4', '0', '1'), independentFiltering = FALSE)


## Combine everything
results_count <- data.frame(rbind(c1vs2, c1vs3, c1vs4, c2vs3, c2vs4, c3vs4))
results_count$Label <- factor(c(rep(c("C2-C1"),dim(Y)[1]),rep(c("C3-C1"),
                                  dim(Y)[1]),rep(c("C4-C1"),dim(Y)[1]),
                                  rep(c("C3-C2"),dim(Y)[1]), rep(c("C4-C2"),dim(Y)[1]), 
                                  rep(c("C4-C3"),dim(Y)[1])))
results_count$protein <- Y$Protein

results_count <- results_count[, c(8,7,1,2,3,4,5,6)]

## Save the result as csv file
write.csv(results_count, file='Supplementary Table ST2.csv')

