##  SwipProteomics

## Comparison: 'Control-Mutant'

#### Method = `MSstatsTMT`

| FDR < 0.05|
|----------:|
|        657|


#### Method = `edgeR + Sum + IRS`
| FDR < 0.05|
|----------:|
|        686|

#### MSstatsTMT Top 5 significant proteins:

|Symbol |     log2FC| PValue| PAdjust| Tstatistic|        SE|       DF|
|:------|----------:|------:|-------:|----------:|---------:|--------:|
|Washc5 | -1.3012495|      0|       0|  -29.97081| 0.0434172| 34.00003|
|Washc1 | -1.8735855|      0|       0|  -22.25587| 0.0841839| 34.00001|
|Washc4 | -1.5169269|      0|       0|  -21.75455| 0.0697292| 34.00001|
|Washc2 | -1.0677318|      0|       0|  -21.67245| 0.0492668| 34.00003|
|Rab21  | -0.5450661|      0|       0|  -16.81386| 0.0324177| 34.00004|



## Comparison: 'Intra-fraction'

#### Method = `MSstats`

| F10| F4| F5| F6| F7| F8| F9|
|---:|--:|--:|--:|--:|--:|--:|
|  22| 26| 24| 20| 27| 30| 34|

Total number of unique proteins (FDR < 0.05): 78

#### Method = `edgeR`

| F10| F4| F5| F6| F7| F8| F9|
|---:|--:|--:|--:|--:|--:|--:|
|   7| 12| 21| 20| 25| 36| 19|

Total number of unique proteins (FDR < 0.05): <NA>


Commonly significant sig prots: 

|Gene   |UniProt|
|:------|:------|
|Rab21  |P35282 |
|Washc4 |Q3UMB9 |
|Washc2 |Q6PGL7 |
|Washc5 |Q8C2E7 |
|Washc1 |Q8VDD8 |


## Rank-Correlation of PValues
```
cor(edger,msstats,method="spearman",use="pairwise.complete")
```

| Control.F#-Mutant.F#| Control-Mutant|
|--------------------:|--------------:|
|            0.6424383|      0.7921496|

## Rank-Correlation of log2fc:

| Control.F#-Mutant.F#| Control-Mutant|
|--------------------:|--------------:|
|            0.8554378|       0.929861|
