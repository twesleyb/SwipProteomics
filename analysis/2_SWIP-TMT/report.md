Loading SwipProteomics
Dividing work into 100 chunks...

Total:21 s


|Protein |Symbol | Entrez|   Mixture|  Genotype| BioFraction| Residuals|
|:-------|:------|------:|---------:|---------:|-----------:|---------:|
|Q8C2E7  |Washc5 | 223593| 0.0004987| 0.9130491|   0.0442397| 0.0422125|
|Q8VDD8  |Washc1 |  68767| 0.0070136| 0.8808449|   0.0436028| 0.0685388|
|Q3UMB9  |Washc4 | 319277| 0.0133109| 0.8722170|   0.0489876| 0.0654845|
|Q6PGL7  |Washc2 |  28006| 0.0000000| 0.7646015|   0.1685491| 0.0668494|
|Q9CR27  |Washc3 |  67282| 0.0205066| 0.6701031|   0.0640885| 0.2453017|
|Q19LI2  |A1bg   | 117586| 0.2186346| 0.5170321|   0.0000000| 0.2643333|

Evaluating Nakagawa goodness-of-fit, refitting modules...


|Protein |Symbol | Entrez|   Mixture|  Genotype| BioFraction| Residuals|  R2.fixef|  R2.total|
|:-------|:------|------:|---------:|---------:|-----------:|---------:|---------:|---------:|
|Q8C2E7  |Washc5 | 223593| 0.0004987| 0.9130491|   0.0442397| 0.0422125| 0.9744804| 0.9762336|
|Q8VDD8  |Washc1 |  68767| 0.0070136| 0.8808449|   0.0436028| 0.0685388| 0.9232585| 0.9298346|
|Q3UMB9  |Washc4 | 319277| 0.0133109| 0.8722170|   0.0489876| 0.0654845| 0.9353344| 0.9494330|
|Q6PGL7  |Washc2 |  28006| 0.0000000| 0.7646015|   0.1685491| 0.0668494| 0.9409087| 0.9409087|
|Q9CR27  |Washc3 |  67282| 0.0205066| 0.6701031|   0.0640885| 0.2453017| 0.7341464| 0.7521275|
|Q19LI2  |A1bg   | 117586| 0.2186346| 0.5170321|   0.0000000| 0.2643333| 0.4804040| 0.7027026|


|Protein |Symbol | Entrez|   Mixture|  Genotype| BioFraction| Residuals|  R2.fixef|  R2.total|
|:-------|:------|------:|---------:|---------:|-----------:|---------:|---------:|---------:|
|P62908  |Rps3   |  27050| 0.0005153| 0.0001460|   0.9971205| 0.0022183| 0.9974401| 0.9979517|
|Q8VDJ3  |Hdlbp  | 110611| 0.0001679| 0.0000000|   0.9968848| 0.0029474| 0.9965136| 0.9966430|
|Q9D8E6  |Rpl4   |  67891| 0.0004759| 0.0003769|   0.9968231| 0.0023241| 0.9970994| 0.9975514|
|Q6ZWN5  |Rps9   |  76846| 0.0000000| 0.0000000|   0.9968041| 0.0031959| 0.9967923| 0.9967923|
|P62281  |Rps11  |  27207| 0.0001360| 0.0001584|   0.9968025| 0.0029032| 0.9970040| 0.9971343|
|Q6ZQ08  |Cnot1  | 234594| 0.0001539| 0.0001600|   0.9967507| 0.0029354| 0.9964280| 0.9965353|

R2 threshold: 0.7


| r2_threshold| out| percent| total| final|
|------------:|---:|-------:|-----:|-----:|
|          0.7| 791|   0.114|  6910|  6119|

Number of proteins with poor fit: 791

WASHC* protein goodness-of-fit statistics:


|Protein |Symbol | Entrez|   Mixture|  Genotype| BioFraction| Residuals|  R2.fixef|  R2.total|
|:-------|:------|------:|---------:|---------:|-----------:|---------:|---------:|---------:|
|Q8C2E7  |Washc5 | 223593| 0.0004987| 0.9130491|   0.0442397| 0.0422125| 0.9744804| 0.9762336|
|Q8VDD8  |Washc1 |  68767| 0.0070136| 0.8808449|   0.0436028| 0.0685388| 0.9232585| 0.9298346|
|Q3UMB9  |Washc4 | 319277| 0.0133109| 0.8722170|   0.0489876| 0.0654845| 0.9353344| 0.9494330|
|Q6PGL7  |Washc2 |  28006| 0.0000000| 0.7646015|   0.1685491| 0.0668494| 0.9409087| 0.9409087|
|Q9CR27  |Washc3 |  67282| 0.0205066| 0.6701031|   0.0640885| 0.2453017| 0.7341464| 0.7521275|
Loading SwipProteomics
Warning message:
Removing 791 proteins with poor fit before building network. 


| samples| proteins|
|-------:|--------:|
|      42|     6119|

Generating protein co-variation network.

Performing network enhancement.

Creating protein-protein interaction network.

PPI graph:


|Edges  |Nodes |
|:------|:-----|
|93,573 |6,119 |

Saved adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved ne_adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved ppi_adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved ppi_adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved ppi_adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved ppi_adjm.rda in /home/twesleyb/projects/SwipProteomics/rdata

Saved norm_prot.rda in /home/twesleyb/projects/SwipProteomics/data
Performing Leidenalg clustering utilizing the SurpriseVertexPartition method to find optimal partition(s).
... Initial partition: Clustering with 6119 elements and 67 clusters.

Splitting 22 modules that contain more than 100 nodes.
... Final partition: Clustering with 6119 elements and 535 clusters.
Loading SwipProteomics

Removing modules that contain less than 5 nodes.

Module statistic(s) used to evaluate module  preservation:
avg.weight, avg.cor, avg.contrib.

Criterion for module preservation: strong.


Evaluating  preservation of Swip modules in the Swip network...
... 295 of 334 Swip modules are preserved in the Swip network.

Loading SwipProteomics


|nProts |kModules |pClustered |medSize |
|:------|:--------|:----------|:-------|
|6,119  |295      |0.893      |15      |

lmer fit to WASH complex (Washc1, Washc2, Washc3, Washc4, Washc5) proteins:
>>>	Abundance ~ 0 + Genotype:BioFraction + (1 | Mixture) + (1 | Protein)


|Term                   | Estimate|    SE|    DF| Tvalue|Pvalue    |
|:----------------------|--------:|-----:|-----:|------:|:---------|
|Control:BioFractionF4  |    6.883| 0.151| 6.891| 45.636|2.91e-09  |
|Control:BioFractionF5  |    7.167| 0.151| 6.891| 47.520|8.249e-10 |
|Control:BioFractionF6  |    7.465| 0.151| 6.891| 49.494|2.289e-09 |
|Control:BioFractionF7  |    7.495| 0.151| 6.891| 49.692|6.248e-10 |
|Control:BioFractionF8  |    7.327| 0.151| 6.891| 48.580|2.017e-09 |
|Control:BioFractionF9  |    7.138| 0.151| 6.891| 47.328|4.722e-10 |
|Control:BioFractionF10 |    7.756| 0.151| 6.891| 51.424|1.931e-09 |
|Mutant:BioFractionF4   |    5.729| 0.151| 6.891| 37.983|4.595e-10 |
|Mutant:BioFractionF5   |    5.933| 0.151| 6.891| 39.334|2.303e-09 |
|Mutant:BioFractionF6   |    6.043| 0.151| 6.891| 40.065|5.369e-10 |
|Mutant:BioFractionF7   |    6.082| 0.151| 6.891| 40.322|2.384e-09 |
|Mutant:BioFractionF8   |    5.927| 0.151| 6.891| 39.299|6.424e-10 |
|Mutant:BioFractionF9   |    5.897| 0.151| 6.891| 39.101|1.991e-09 |
|Mutant:BioFractionF10  |    6.054| 0.151| 6.891| 40.142|3.631e-10 |
R2m: Marginal; variation explained by fixed effects.
R2c: Conditional; total variation explained by the model.


|       R2m|       R2c|
|---------:|---------:|
| 0.7620866| 0.8928053|


|Contrast       |    log2FC| percentControl| Pvalue| Tstatistic|        SE|  DF| nProteins|
|:--------------|---------:|--------------:|------:|----------:|---------:|---:|---------:|
|Mutant-Control | -1.366663|      0.3877871|      0|  -36.94728| 0.0369896| 190|         5|

Assessing module-level contrasts with lmerTest.

Time to analyze 295 modules:
Time difference of 5.356059 secs
Warning message:
0 modules with singular fits will be removed. 

Final number of modules : 295

Final percent clustered : 0.893

Final Median module size: 15

Washc4 assigned to module: M7


|Module | Size|Contrast       |     log2FC| percentControl|    Pvalue| Tstatistic|        SE|        DF|       FDR|   Padjust|
|:------|----:|:--------------|----------:|--------------:|---------:|----------:|---------:|---------:|---------:|---------:|
|M7     |   56|Mutant-Control | -0.2742921|      0.8268559| 0.0000000| -15.873873| 0.0172795| 2281.0000| 0.0000000| 0.0000000|
|M11    |   22|Mutant-Control |  0.2736723|      1.2088811| 0.0000000|  14.870168| 0.0184041|  887.0000| 0.0000000| 0.0000000|
|M77    |   55|Mutant-Control |  0.1030176|      1.0740176| 0.0000000|  11.670933| 0.0088268| 2240.0000| 0.0000000| 0.0000000|
|M192   |   14|Mutant-Control | -0.3763005|      0.7704106| 0.0000000| -11.384577| 0.0330535|  559.0003| 0.0000000| 0.0000000|
|M198   |   11|Mutant-Control | -0.4774924|      0.7182249| 0.0000000| -11.451096| 0.0416984|  435.9998| 0.0000000| 0.0000000|
|M197   |   11|Mutant-Control | -0.1994586|      0.8708773| 0.0000000| -11.263643| 0.0177082|  438.0000| 0.0000000| 0.0000000|
|M272   |   14|Mutant-Control |  0.2888011|      1.2216247| 0.0000000|  10.453953| 0.0276260|  559.0000| 0.0000000| 0.0000000|
|M190   |   19|Mutant-Control | -0.1597226|      0.8951972| 0.0000000|  -9.917399| 0.0161053|  764.0000| 0.0000000| 0.0000000|
|M268   |   30|Mutant-Control |  0.1631997|      1.1197679| 0.0000000|   9.015929| 0.0181013| 1214.9995| 0.0000000| 0.0000000|
|M211   |   33|Mutant-Control |  0.0662855|      1.0470175| 0.0000000|   8.836857| 0.0075010| 1338.0000| 0.0000000| 0.0000000|
|M189   |   21|Mutant-Control | -0.1132031|      0.9245331| 0.0000000|  -8.798891| 0.0128656|  846.0004| 0.0000000| 0.0000000|
|M142   |   17|Mutant-Control |  0.1936104|      1.1436221| 0.0000000|   8.778039| 0.0220562|  682.0000| 0.0000000| 0.0000000|
|M18    |   37|Mutant-Control | -0.0735388|      0.9503041| 0.0000000|  -8.540361| 0.0086107| 1501.9990| 0.0000000| 0.0000000|
|M32    |   20|Mutant-Control | -0.1019988|      0.9317412| 0.0000000|  -8.540805| 0.0119425|  804.9999| 0.0000000| 0.0000000|
|M161   |   11|Mutant-Control | -0.1631785|      0.8930554| 0.0000000|  -8.698029| 0.0187604|  436.0000| 0.0000000| 0.0000000|
|M281   |   15|Mutant-Control |  0.1367084|      1.0993939| 0.0000000|   8.503607| 0.0160765|  600.0000| 0.0000000| 0.0000000|
|M291   |   18|Mutant-Control |  0.0913319|      1.0653533| 0.0000000|   8.241355| 0.0110821|  723.0000| 0.0000000| 0.0000000|
|M134   |   10|Mutant-Control |  0.1666064|      1.1224152| 0.0000000|   8.287066| 0.0201044|  395.0000| 0.0000000| 0.0000000|
|M228   |   11|Mutant-Control | -0.1322447|      0.9124107| 0.0000000|  -8.210687| 0.0161064|  436.0000| 0.0000000| 0.0000000|
|M202   |   23|Mutant-Control | -0.0924249|      0.9379449| 0.0000000|  -7.920189| 0.0116695|  928.0001| 0.0000000| 0.0000000|
|M29    |   22|Mutant-Control |  0.0987445|      1.0708412| 0.0000000|   7.351494| 0.0134319|  887.0000| 0.0000000| 0.0000000|
|M147   |   12|Mutant-Control | -0.1249699|      0.9170232| 0.0000000|  -7.114749| 0.0175649|  477.0000| 0.0000000| 0.0000000|
|M154   |   22|Mutant-Control | -0.0902188|      0.9393803| 0.0000000|  -7.015441| 0.0128600|  887.0000| 0.0000000| 0.0000000|
|M22    |   30|Mutant-Control |  0.0845672|      1.0603696| 0.0000000|   6.935346| 0.0121937| 1215.0000| 0.0000000| 0.0000000|
|M129   |   20|Mutant-Control |  0.1016411|      1.0729933| 0.0000000|   6.917460| 0.0146934|  805.0000| 0.0000000| 0.0000000|
|M284   |   10|Mutant-Control |  0.2258668|      1.1694797| 0.0000000|   7.017590| 0.0321858|  395.0000| 0.0000000| 0.0000000|
|M117   |   15|Mutant-Control | -0.0660552|      0.9552464| 0.0000000|  -6.814190| 0.0096938|  600.0000| 0.0000000| 0.0000000|
|M227   |   12|Mutant-Control | -0.1409830|      0.9069010| 0.0000000|  -6.807196| 0.0207109|  479.0000| 0.0000000| 0.0000000|
|M223   |   21|Mutant-Control |  0.1007160|      1.0723055| 0.0000000|   6.718123| 0.0149917|  846.0001| 0.0000000| 0.0000000|
|M111   |   26|Mutant-Control | -0.0835975|      0.9437015| 0.0000000|  -6.655153| 0.0125613| 1051.0001| 0.0000000| 0.0000000|
|M102   |   17|Mutant-Control |  0.0790571|      1.0563274| 0.0000000|   6.564923| 0.0120423|  684.0000| 0.0000000| 0.0000000|
|M135   |   10|Mutant-Control |  0.1085429|      1.0781388| 0.0000000|   6.463465| 0.0167933|  394.9999| 0.0000000| 0.0000001|
|M283   |   14|Mutant-Control |  0.1043699|      1.0750248| 0.0000000|   6.149421| 0.0169723|  559.0003| 0.0000000| 0.0000004|
|M196   |   12|Mutant-Control | -0.0803814|      0.9458076| 0.0000000|  -6.165493| 0.0130373|  477.0000| 0.0000000| 0.0000004|
|M126   |   23|Mutant-Control | -0.0865374|      0.9417804| 0.0000000|  -5.943666| 0.0145596|  927.9999| 0.0000000| 0.0000012|
|M195   |   12|Mutant-Control | -0.0935300|      0.9372268| 0.0000000|  -5.981504| 0.0156365|  477.0000| 0.0000000| 0.0000013|
|M200   |   25|Mutant-Control |  0.0639949|      1.0453564| 0.0000000|   5.884404| 0.0108753| 1009.9999| 0.0000000| 0.0000016|
|M4     |   71|Mutant-Control | -0.0395735|      0.9729425| 0.0000000|  -5.767548| 0.0068614| 2896.0001| 0.0000001| 0.0000026|
|M199   |   10|Mutant-Control | -0.1463641|      0.9035247| 0.0000000|  -5.856844| 0.0249903|  395.0002| 0.0000001| 0.0000029|
|M6     |   69|Mutant-Control | -0.0685942|      0.9535667| 0.0000001|  -5.462365| 0.0125576| 2814.0001| 0.0000004| 0.0000151|
|M193   |   13|Mutant-Control | -0.0933278|      0.9373581| 0.0000002|  -5.303853| 0.0175962|  518.0000| 0.0000012| 0.0000496|
|M170   |    6|Mutant-Control | -0.1517568|      0.9001537| 0.0000004|  -5.230647| 0.0290130|  231.0000| 0.0000027| 0.0001116|
|M203   |   23|Mutant-Control | -0.0433783|      0.9703800| 0.0000007|  -4.997512| 0.0086800|  928.0001| 0.0000048| 0.0002047|
|M140   |   22|Mutant-Control | -0.0544643|      0.9629520| 0.0000010|  -4.921234| 0.0110672|  887.0000| 0.0000069| 0.0003024|
|M152   |    9|Mutant-Control | -0.1060666|      0.9291178| 0.0000023|  -4.801975| 0.0220881|  354.0000| 0.0000152| 0.0006852|
|M97    |   11|Mutant-Control |  0.0589365|      1.0416976| 0.0000025|   4.773507| 0.0123466|  436.0000| 0.0000159| 0.0007302|
|M156   |   19|Mutant-Control |  0.0542887|      1.0383470| 0.0000026|   4.736345| 0.0114622|  764.0000| 0.0000163| 0.0007653|
|M212   |   33|Mutant-Control |  0.0397198|      1.0279142| 0.0000045|   4.607038| 0.0086215| 1338.0002| 0.0000275| 0.0013197|
|M112   |   21|Mutant-Control | -0.0360643|      0.9753120| 0.0000053|  -4.583179| 0.0078688|  848.0000| 0.0000317| 0.0015540|
|M149   |   12|Mutant-Control | -0.0617257|      0.9581173| 0.0000064|  -4.563888| 0.0135248|  477.0000| 0.0000377| 0.0018869|
|M295   |    7|Mutant-Control |  0.0936763|      1.0670859| 0.0000073|   4.573833| 0.0204809|  272.0001| 0.0000421| 0.0021482|
|M162   |   11|Mutant-Control | -0.1108000|      0.9260744| 0.0000081|  -4.516445| 0.0245326|  436.0001| 0.0000460| 0.0023912|
|M86    |   15|Mutant-Control | -0.0842771|      0.9432571| 0.0000085|  -4.490707| 0.0187670|  600.0000| 0.0000474| 0.0025121|
|M148   |   12|Mutant-Control |  0.0804122|      1.0573201| 0.0000097|   4.472579| 0.0179789|  477.0001| 0.0000528| 0.0028528|
|M35    |   12|Mutant-Control | -0.0469925|      0.9679521| 0.0000141|  -4.387404| 0.0107108|  479.0000| 0.0000757| 0.0041651|
|M21    |   30|Mutant-Control | -0.0442885|      0.9697680| 0.0000168|  -4.321384| 0.0102487| 1215.0000| 0.0000873| 0.0049477|
|M191   |   18|Mutant-Control | -0.0656675|      0.9555031| 0.0000170|  -4.330149| 0.0151652|  723.0001| 0.0000873| 0.0050164|
|M3     |   75|Mutant-Control | -0.0398810|      0.9727352| 0.0000172|  -4.305734| 0.0092623| 3060.0001| 0.0000873| 0.0050640|
|M289   |    6|Mutant-Control |  0.1085978|      1.0781798| 0.0000269|   4.284718| 0.0253454|  231.0000| 0.0001343| 0.0079246|
|M39    |   89|Mutant-Control | -0.0463964|      0.9683521| 0.0000331|  -4.156348| 0.0111628| 3634.0007| 0.0001627| 0.0097595|
|M286   |    8|Mutant-Control |  0.0598141|      1.0423315| 0.0000352|   4.197231| 0.0142509|  313.0000| 0.0001705| 0.0103977|
|M10    |   24|Mutant-Control |  0.0928730|      1.0664919| 0.0000373|   4.142961| 0.0224170|  969.0001| 0.0001774| 0.0109974|
|M222   |   22|Mutant-Control |  0.0542679|      1.0383320| 0.0000407|   4.124312| 0.0131580|  887.0005| 0.0001905| 0.0120004|
|M171   |   36|Mutant-Control |  0.0297809|      1.0208570| 0.0000563|   4.039842| 0.0073718| 1461.0000| 0.0002594| 0.0165990|
|M229   |   10|Mutant-Control |  0.0727991|      1.0517553| 0.0000724|   4.010358| 0.0181528|  397.0000| 0.0003288| 0.0213705|
|M178   |   12|Mutant-Control |  0.0518544|      1.0365964| 0.0000923|   3.943658| 0.0131488|  477.0000| 0.0004125| 0.0272247|
|M79    |   36|Mutant-Control |  0.0495771|      1.0349615| 0.0001118|   3.874006| 0.0127974| 1461.0004| 0.0004922| 0.0329799|
|M256   |   17|Mutant-Control |  0.0613920|      1.0434720| 0.0001537|   3.806329| 0.0161289|  682.0011| 0.0006669| 0.0453517|
|M274   |    9|Mutant-Control |  0.0991000|      1.0711051| 0.0001590|   3.817753| 0.0259577|  354.0000| 0.0006796| 0.0468909|

Modules with greater than 10% change:


|Module | Size|Contrast       |     log2FC| percentControl| Pvalue| Tstatistic|        SE|        DF|     FDR|   Padjust|
|:------|----:|:--------------|----------:|--------------:|------:|----------:|---------:|---------:|-------:|---------:|
|M7     |   56|Mutant-Control | -0.2742921|      0.8268559|  0e+00| -15.873873| 0.0172795| 2281.0000| 0.0e+00| 0.0000000|
|M11    |   22|Mutant-Control |  0.2736723|      1.2088811|  0e+00|  14.870168| 0.0184041|  887.0000| 0.0e+00| 0.0000000|
|M192   |   14|Mutant-Control | -0.3763005|      0.7704106|  0e+00| -11.384577| 0.0330535|  559.0003| 0.0e+00| 0.0000000|
|M198   |   11|Mutant-Control | -0.4774924|      0.7182249|  0e+00| -11.451096| 0.0416984|  435.9998| 0.0e+00| 0.0000000|
|M197   |   11|Mutant-Control | -0.1994586|      0.8708773|  0e+00| -11.263643| 0.0177082|  438.0000| 0.0e+00| 0.0000000|
|M272   |   14|Mutant-Control |  0.2888011|      1.2216247|  0e+00|  10.453953| 0.0276260|  559.0000| 0.0e+00| 0.0000000|
|M190   |   19|Mutant-Control | -0.1597226|      0.8951972|  0e+00|  -9.917399| 0.0161053|  764.0000| 0.0e+00| 0.0000000|
|M268   |   30|Mutant-Control |  0.1631997|      1.1197679|  0e+00|   9.015929| 0.0181013| 1214.9995| 0.0e+00| 0.0000000|
|M142   |   17|Mutant-Control |  0.1936104|      1.1436221|  0e+00|   8.778039| 0.0220562|  682.0000| 0.0e+00| 0.0000000|
|M161   |   11|Mutant-Control | -0.1631785|      0.8930554|  0e+00|  -8.698029| 0.0187604|  436.0000| 0.0e+00| 0.0000000|
|M134   |   10|Mutant-Control |  0.1666064|      1.1224152|  0e+00|   8.287066| 0.0201044|  395.0000| 0.0e+00| 0.0000000|
|M284   |   10|Mutant-Control |  0.2258668|      1.1694797|  0e+00|   7.017590| 0.0321858|  395.0000| 0.0e+00| 0.0000000|
|M227   |   12|Mutant-Control | -0.1409830|      0.9069010|  0e+00|  -6.807196| 0.0207109|  479.0000| 0.0e+00| 0.0000000|
|M199   |   10|Mutant-Control | -0.1463641|      0.9035247|  0e+00|  -5.856844| 0.0249903|  395.0002| 1.0e-07| 0.0000029|
|M170   |    6|Mutant-Control | -0.1517568|      0.9001537|  4e-07|  -5.230647| 0.0290130|  231.0000| 2.7e-06| 0.0001116|
Number of significant modules (Bonferroni<0.05): 69
Loading SwipProteomics

Evaluating goodness-of-fit of modules.
  |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%boundary (singular) fit: see ?isSingular
  |                                                                              |===                                                                   |   4%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%boundary (singular) fit: see ?isSingular
  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |======                                                                |   8%boundary (singular) fit: see ?isSingular
  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%boundary (singular) fit: see ?isSingular
  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%boundary (singular) fit: see ?isSingular
  |                                                                              |========                                                              |  12%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%boundary (singular) fit: see ?isSingular
  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%boundary (singular) fit: see ?isSingular
  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  17%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |============                                                          |  18%boundary (singular) fit: see ?isSingular
  |                                                                              |=============                                                         |  18%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=============                                                         |  19%boundary (singular) fit: see ?isSingular
  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%boundary (singular) fit: see ?isSingular
  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |================                                                      |  22%boundary (singular) fit: see ?isSingular
  |                                                                              |================                                                      |  23%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=================                                                     |  24%boundary (singular) fit: see ?isSingular
  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%boundary (singular) fit: see ?isSingular
  |                                                                              |=====================                                                 |  29%boundary (singular) fit: see ?isSingular
  |                                                                              |=====================                                                 |  30%boundary (singular) fit: see ?isSingular
  |                                                                              |=====================                                                 |  31%boundary (singular) fit: see ?isSingular
  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%boundary (singular) fit: see ?isSingular
  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |========================                                              |  35%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |==========================                                            |  37%boundary (singular) fit: see ?isSingular
  |                                                                              |==========================                                            |  38%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===========================                                           |  38%boundary (singular) fit: see ?isSingular
  |                                                                              |===========================                                           |  39%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |============================                                          |  41%boundary (singular) fit: see ?isSingular
  |                                                                              |=============================                                         |  41%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=============================                                         |  42%boundary (singular) fit: see ?isSingular
  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%boundary (singular) fit: see ?isSingular
  |                                                                              |===============================                                       |  44%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%boundary (singular) fit: see ?isSingular
  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%boundary (singular) fit: see ?isSingular
  |                                                                              |==================================                                    |  49%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%boundary (singular) fit: see ?isSingular
  |                                                                              |=======================================                               |  55%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=======================================                               |  56%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |========================================                              |  57%boundary (singular) fit: see ?isSingular
  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=========================================                             |  59%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |==========================================                            |  59%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |==========================================                            |  60%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |==========================================                            |  61%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%boundary (singular) fit: see ?isSingular
  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%boundary (singular) fit: see ?isSingular
  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%boundary (singular) fit: see ?isSingular
  |                                                                              |==================================================                    |  71%boundary (singular) fit: see ?isSingular
  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |======================================================                |  77%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |======================================================                |  78%boundary (singular) fit: see ?isSingular
  |                                                                              |=======================================================               |  78%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=======================================================               |  79%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%boundary (singular) fit: see ?isSingular
  |                                                                              |========================================================              |  81%boundary (singular) fit: see ?isSingular
  |                                                                              |=========================================================             |  81%boundary (singular) fit: see ?isSingular
  |                                                                              |=========================================================             |  82%boundary (singular) fit: see ?isSingular
  |                                                                              |==========================================================            |  82%boundary (singular) fit: see ?isSingular
  |                                                                              |==========================================================            |  83%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===========================================================           |  84%boundary (singular) fit: see ?isSingular
  |                                                                              |===========================================================           |  85%boundary (singular) fit: see ?isSingular
  |                                                                              |============================================================          |  85%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=============================================================         |  88%boundary (singular) fit: see ?isSingular
  |                                                                              |==============================================================        |  88%boundary (singular) fit: see ?isSingular
  |                                                                              |==============================================================        |  89%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===============================================================       |  89%boundary (singular) fit: see ?isSingular
  |                                                                              |===============================================================       |  90%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%boundary (singular) fit: see ?isSingular
  |                                                                              |=================================================================     |  93%boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
boundary (singular) fit: see ?isSingular
  |                                                                              |=================================================================     |  94%boundary (singular) fit: see ?isSingular
  |                                                                              |==================================================================    |  94%boundary (singular) fit: see ?isSingular
  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  97%boundary (singular) fit: see ?isSingular
  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%boundary (singular) fit: see ?isSingular
  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%
There were problems fitting 0 models.
Loading SwipProteomics

Performing GSE analysis for all modules:
  |                                                                              |                                                                      |   0%  |                                                                              |                                                                      |   1%  |                                                                              |=                                                                     |   1%  |                                                                              |=                                                                     |   2%  |                                                                              |==                                                                    |   2%  |                                                                              |==                                                                    |   3%  |                                                                              |==                                                                    |   4%  |                                                                              |===                                                                   |   4%  |                                                                              |===                                                                   |   5%  |                                                                              |====                                                                  |   5%  |                                                                              |====                                                                  |   6%  |                                                                              |=====                                                                 |   6%  |                                                                              |=====                                                                 |   7%  |                                                                              |=====                                                                 |   8%  |                                                                              |======                                                                |   8%  |                                                                              |======                                                                |   9%  |                                                                              |=======                                                               |   9%  |                                                                              |=======                                                               |  10%  |                                                                              |=======                                                               |  11%  |                                                                              |========                                                              |  11%  |                                                                              |========                                                              |  12%  |                                                                              |=========                                                             |  12%  |                                                                              |=========                                                             |  13%  |                                                                              |=========                                                             |  14%  |                                                                              |==========                                                            |  14%  |                                                                              |==========                                                            |  15%  |                                                                              |===========                                                           |  15%  |                                                                              |===========                                                           |  16%  |                                                                              |============                                                          |  16%  |                                                                              |============                                                          |  17%  |                                                                              |============                                                          |  18%  |                                                                              |=============                                                         |  18%  |                                                                              |=============                                                         |  19%  |                                                                              |==============                                                        |  19%  |                                                                              |==============                                                        |  20%  |                                                                              |==============                                                        |  21%  |                                                                              |===============                                                       |  21%  |                                                                              |===============                                                       |  22%  |                                                                              |================                                                      |  22%  |                                                                              |================                                                      |  23%  |                                                                              |================                                                      |  24%  |                                                                              |=================                                                     |  24%  |                                                                              |=================                                                     |  25%  |                                                                              |==================                                                    |  25%  |                                                                              |==================                                                    |  26%  |                                                                              |===================                                                   |  26%  |                                                                              |===================                                                   |  27%  |                                                                              |===================                                                   |  28%  |                                                                              |====================                                                  |  28%  |                                                                              |====================                                                  |  29%  |                                                                              |=====================                                                 |  29%  |                                                                              |=====================                                                 |  30%  |                                                                              |=====================                                                 |  31%  |                                                                              |======================                                                |  31%  |                                                                              |======================                                                |  32%  |                                                                              |=======================                                               |  32%  |                                                                              |=======================                                               |  33%  |                                                                              |=======================                                               |  34%  |                                                                              |========================                                              |  34%  |                                                                              |========================                                              |  35%  |                                                                              |=========================                                             |  35%  |                                                                              |=========================                                             |  36%  |                                                                              |==========================                                            |  36%  |                                                                              |==========================                                            |  37%  |                                                                              |==========================                                            |  38%  |                                                                              |===========================                                           |  38%  |                                                                              |===========================                                           |  39%  |                                                                              |============================                                          |  39%  |                                                                              |============================                                          |  40%  |                                                                              |============================                                          |  41%  |                                                                              |=============================                                         |  41%  |                                                                              |=============================                                         |  42%  |                                                                              |==============================                                        |  42%  |                                                                              |==============================                                        |  43%  |                                                                              |==============================                                        |  44%  |                                                                              |===============================                                       |  44%  |                                                                              |===============================                                       |  45%  |                                                                              |================================                                      |  45%  |                                                                              |================================                                      |  46%  |                                                                              |=================================                                     |  46%  |                                                                              |=================================                                     |  47%  |                                                                              |=================================                                     |  48%  |                                                                              |==================================                                    |  48%  |                                                                              |==================================                                    |  49%  |                                                                              |===================================                                   |  49%  |                                                                              |===================================                                   |  50%  |                                                                              |===================================                                   |  51%  |                                                                              |====================================                                  |  51%  |                                                                              |====================================                                  |  52%  |                                                                              |=====================================                                 |  52%  |                                                                              |=====================================                                 |  53%  |                                                                              |=====================================                                 |  54%  |                                                                              |======================================                                |  54%  |                                                                              |======================================                                |  55%  |                                                                              |=======================================                               |  55%  |                                                                              |=======================================                               |  56%  |                                                                              |========================================                              |  56%  |                                                                              |========================================                              |  57%  |                                                                              |========================================                              |  58%  |                                                                              |=========================================                             |  58%  |                                                                              |=========================================                             |  59%  |                                                                              |==========================================                            |  59%  |                                                                              |==========================================                            |  60%  |                                                                              |==========================================                            |  61%  |                                                                              |===========================================                           |  61%  |                                                                              |===========================================                           |  62%  |                                                                              |============================================                          |  62%  |                                                                              |============================================                          |  63%  |                                                                              |============================================                          |  64%  |                                                                              |=============================================                         |  64%  |                                                                              |=============================================                         |  65%  |                                                                              |==============================================                        |  65%  |                                                                              |==============================================                        |  66%  |                                                                              |===============================================                       |  66%  |                                                                              |===============================================                       |  67%  |                                                                              |===============================================                       |  68%  |                                                                              |================================================                      |  68%  |                                                                              |================================================                      |  69%  |                                                                              |=================================================                     |  69%  |                                                                              |=================================================                     |  70%  |                                                                              |=================================================                     |  71%  |                                                                              |==================================================                    |  71%  |                                                                              |==================================================                    |  72%  |                                                                              |===================================================                   |  72%  |                                                                              |===================================================                   |  73%  |                                                                              |===================================================                   |  74%  |                                                                              |====================================================                  |  74%  |                                                                              |====================================================                  |  75%  |                                                                              |=====================================================                 |  75%  |                                                                              |=====================================================                 |  76%  |                                                                              |======================================================                |  76%  |                                                                              |======================================================                |  77%  |                                                                              |======================================================                |  78%  |                                                                              |=======================================================               |  78%  |                                                                              |=======================================================               |  79%  |                                                                              |========================================================              |  79%  |                                                                              |========================================================              |  80%  |                                                                              |========================================================              |  81%  |                                                                              |=========================================================             |  81%  |                                                                              |=========================================================             |  82%  |                                                                              |==========================================================            |  82%  |                                                                              |==========================================================            |  83%  |                                                                              |==========================================================            |  84%  |                                                                              |===========================================================           |  84%  |                                                                              |===========================================================           |  85%  |                                                                              |============================================================          |  85%  |                                                                              |============================================================          |  86%  |                                                                              |=============================================================         |  86%  |                                                                              |=============================================================         |  87%  |                                                                              |=============================================================         |  88%  |                                                                              |==============================================================        |  88%  |                                                                              |==============================================================        |  89%  |                                                                              |===============================================================       |  89%  |                                                                              |===============================================================       |  90%  |                                                                              |===============================================================       |  91%  |                                                                              |================================================================      |  91%  |                                                                              |================================================================      |  92%  |                                                                              |=================================================================     |  92%  |                                                                              |=================================================================     |  93%  |                                                                              |=================================================================     |  94%  |                                                                              |==================================================================    |  94%  |                                                                              |==================================================================    |  95%  |                                                                              |===================================================================   |  95%  |                                                                              |===================================================================   |  96%  |                                                                              |====================================================================  |  96%  |                                                                              |====================================================================  |  97%  |                                                                              |====================================================================  |  98%  |                                                                              |===================================================================== |  98%  |                                                                              |===================================================================== |  99%  |                                                                              |======================================================================|  99%  |                                                                              |======================================================================| 100%There were 50 or more warnings (use warnings() to see the first 50)


Number of modules with something interesting going on: 117

Significant Modules with significant gse:
(29 of 69 significant modules.)


|Module |TopPathway                                                    |         FE|   Padjust|
|:------|:-------------------------------------------------------------|----------:|---------:|
|M10    |CORUM: Class C Vps complex (hVPS11, hVPS18, hVPS16, rVPS33a ) | 227.750000| 0.0000001|
|M11    |CORUM: Fibrinogen complex                                     | 248.454545| 0.0000167|
|M117   |CORUM: Arp2/3 protein complex                                 | 156.171429| 0.0001716|
|M134   |LopitDC: LYSOSOME                                             |  27.440000| 0.0416553|
|M135   |LopitDC: LYSOSOME                                             |  36.586667| 0.0007617|
|M140   |CORUM: Protein-sorting complex (Stam1, Hgs, Eps15)            | 165.666667| 0.0136492|
|M154   |LopitDC: PM                                                   |   6.708007| 0.0032901|
|M162   |Uezu et al., 2016: ePSD                                       |  21.858054| 0.0004395|
|M170   |Uezu et al., 2016: ePSD                                       |  24.043860| 0.0493791|
|M171   |CORUM: Ribosome, cytoplasmic                                  |  52.797068| 0.0000000|
|M178   |LopitDC: RIBOSOME                                             |  22.333333| 0.0000272|
|M18    |LopitDC: CYTOSOL                                              |   4.629225| 0.0036394|
|M189   |CORUM: IGF1R-CXCR4-GNAI2-GNB1 complex                         | 130.166667| 0.0247622|
|M191   |CORUM: GABA-A receptor (GABRA1, GABRB2, GABRD)                | 202.444444| 0.0090481|
|M196   |CORUM: AMPA receptor complex (anti-GluA1-a)                   |  68.375000| 0.0026574|
|M211   |LopitDC: ER                                                   |  12.153959| 0.0000000|
|M212   |CORUM: SPG3A-SPG4-SPG31 complex                               | 110.424242| 0.0311675|
|M222   |CORUM: SNARE complex (Stx5, Gosr2, Sec22b, Bet1)              | 124.250000| 0.0272317|
|M227   |CORUM: SNARE complex (Snap25, Vamp3, Vamp2, Napa, Stx12)      | 182.266667| 0.0129786|
|M274   |CORUM: HOPS complex                                           | 202.444444| 0.0106292|
|M286   |LopitDC: LYSOSOME                                             |  57.166667| 0.0000021|
|M29    |LopitDC: CYTOSOL                                              |   6.369967| 0.0011682|
|M3     |Uezu et al., 2016: ePSD                                       |   6.411696| 0.0007998|
|M39    |LopitDC: MITOCHONDRION                                        |   9.363271| 0.0000000|
|M4     |CORUM: ESCRT-II complex                                       |  76.985916| 0.0006198|
|M6     |LopitDC: MITOCHONDRION                                        |   4.326184| 0.0000001|
|M7     |CORUM: CCC-Wash (WASH1, FAM21C) complex                       |  76.733418| 0.0000000|
|M77    |LopitDC: ER                                                   |  12.559091| 0.0000000|
|M79    |LopitDC: ER                                                   |   9.284274| 0.0000000|

Modules with significant LopitDC gse:


|Pathway                    |TopModule |        FE| Padjust|
|:--------------------------|:---------|---------:|-------:|
|LopitDC: CYTOSOL           |M16       |  8.953343| 0.0e+00|
|LopitDC: ER                |M77       | 12.559091| 0.0e+00|
|LopitDC: GA                |M82       | 45.185185| 0.0e+00|
|LopitDC: LYSOSOME          |M286      | 57.166667| 2.1e-06|
|LopitDC: MITOCHONDRION     |M39       |  9.363271| 0.0e+00|
|LopitDC: NUCLEUS/CHROMATIN |M59       |  3.401832| 1.0e-07|
|LopitDC: PEROXISOME        |M145      | 88.452525| 0.0e+00|
|LopitDC: PM                |M155      |  9.662725| 6.0e-07|
|LopitDC: PROTEASOME        |M89       | 58.101535| 0.0e+00|
|LopitDC: RIBOSOME          |M171      | 31.018518| 0.0e+00|


|Module |TopPathway                                                                                                                        |         FE|       FDR|
|:------|:---------------------------------------------------------------------------------------------------------------------------------|----------:|---------:|
|M169   |CORUM: LIN2-LIN7-SAP97 complex                                                                                                    | 520.571429| 0.0012436|
|M167   |CORUM: Xin-Cdh2-Ctnnb1-Ctnnd1 complex                                                                                             | 390.500000| 0.0024847|
|M151   |CORUM: Dystrobrevin-syntrophin complex, brain-derived                                                                             | 372.681818| 0.0000071|
|M262   |CORUM: NELF complex (Negative elongation factor complex)                                                                          | 341.687500| 0.0033121|
|M287   |CORUM: Phosphatidylinositol 4-kinase complex                                                                                      | 341.687500| 0.0033121|
|M218   |CORUM: SNARE complex (RINT1, ZW10, p31, Stx18)                                                                                    | 341.625000| 0.0000095|
|M106   |CORUM: Coatomer complex                                                                                                           | 325.476190| 0.0000000|
|M107   |CORUM: MutS-alpha-PK-zeta complex                                                                                                 | 303.722222| 0.0039045|
|M96    |CORUM: CAND1-CUL2-RBX1 complex                                                                                                    | 303.722222| 0.0039045|
|M279   |CORUM: Gamma-tubulin complex                                                                                                      | 303.666667| 0.0000001|
|M172   |CORUM: eIF3 complex (EIF3S6, EIF3S5, EIF3S4, EIF3S3, EIF3S6IP, EIF3S2, EIF3S9, EIF3S12,  EIF3S10, EIF3S8,  EIF3S1, EIF3S7, PCID1) | 265.603239| 0.0000000|
|M109   |CORUM: COG1-COG8-COG5-COG6-COG7 subcomplex                                                                                        | 218.640000| 0.0088620|
|M163   |CORUM: Exocyst complex                                                                                                            | 205.087500| 0.0000724|
|M99    |CORUM: Multisynthetase complex                                                                                                    | 202.444444| 0.0000000|
|M253   |CORUM: AFF4 super elongation complex (SEC)                                                                                        | 199.054545| 0.0107958|
|M61    |CORUM: p54(nrb)-PSF-matrin3 complex                                                                                               | 195.214286| 0.0000355|
|M184   |CORUM: DGCR8 multiprotein complex                                                                                                 | 182.233333| 0.0066376|
|M290   |CORUM: NCKAP1-WASF3-RAC1 complex                                                                                                  | 173.523809| 0.0124144|
|M104   |CORUM: XFIM complex                                                                                                               | 168.246154| 0.0153327|
|M187   |CORUM: GNAI1-GNB2-GNG12 complex                                                                                                   | 158.434783| 0.0149527|
|M174   |CORUM: DNA-PK-Ku-eIF2-NF90-NF45 complex                                                                                           | 136.700000| 0.0002737|
|M240   |CORUM: PU.1-Sin3A-Hdac-MeCP2 complex                                                                                              | 130.166667| 0.0247622|
|M53    |CORUM: MIB complex                                                                                                                | 124.227273| 0.0302169|
|M232   |CORUM: GATOR2 complex                                                                                                             | 115.073684| 0.0335648|
|M63    |CORUM: LAS1L-PELP1-TEX10-WDR18-NOL9-SENP3  complex                                                                                | 113.875000| 0.0004352|
|M94    |CORUM: Arp2/3 protein complex                                                                                                     | 111.551020| 0.0187345|
|M113   |CORUM: BBS1-BBS4-BBS5-PKD1-TTC8 complex                                                                                           | 109.340000| 0.0372669|
|M60    |CORUM: C1q complex                                                                                                                | 101.222222| 0.0371748|
|M259   |CORUM: Set1A complex                                                                                                              |  99.454545| 0.0483858|
|M276   |CORUM: hASC-1 complex                                                                                                             |  97.625000| 0.0444956|
|M254   |CORUM: 40S ribosomal subunit, cytoplasmic                                                                                         |  97.526754| 0.0000000|
|M68    |CORUM: 17S U2 snRNP                                                                                                               |  96.529412| 0.0000000|
|M67    |CORUM: CCR4-NOT-CNOT8-CNOT6 complex                                                                                               |  95.947368| 0.0008697|
|M62    |CORUM: TRAPP complex                                                                                                              |  93.754286| 0.0000000|
|M277   |Takamori et al., 2006: Transporter / Channel                                                                                      |  90.530827| 0.0000000|
|M145   |LopitDC: PEROXISOME                                                                                                               |  88.452525| 0.0000000|
|M83    |CORUM: Oligosaccharyltransferase complex (Stt3B variant)                                                                          |  82.195489| 0.0350970|
|M114   |Takamori et al., 2006: Endocytosis- related proteins                                                                              |  82.005000| 0.0014586|
|M90    |CORUM: CCT-Pdcl complex                                                                                                           |  68.431925| 0.0000000|
|M143   |LopitDC: PEROXISOME                                                                                                               |  68.290553| 0.0000000|
|M89    |CORUM: PA700-20S-PA28 complex                                                                                                     |  59.166667| 0.0000000|
|M82    |LopitDC: GA                                                                                                                       |  45.185185| 0.0000000|
|M27    |CORUM: Kinase maturation complex 1                                                                                                |  42.734375| 0.0118131|
|M186   |CORUM: Spliceosome                                                                                                                |  34.306250| 0.0000682|
|M73    |CORUM: Spliceosome                                                                                                                |  32.162109| 0.0000000|
|M74    |CORUM: Parvulin-associated pre-rRNP complex                                                                                       |  30.932203| 0.0070343|
|M92    |Takamori et al., 2006: Metabolic enzymes                                                                                          |  28.489583| 0.0027867|
|M2     |CORUM: AMPA receptor complex (anti-GluA2-a)                                                                                       |  23.819686| 0.0003561|
|M179   |CORUM: Spliceosome                                                                                                                |  19.492188| 0.0001110|
|M173   |CORUM: Spliceosome                                                                                                                |  17.657629| 0.0000038|
|M175   |CORUM: Spliceosome                                                                                                                |  15.315290| 0.0003977|
|M59    |CORUM: Spliceosome                                                                                                                |  13.907939| 0.0000000|
|M65    |CORUM: Spliceosome                                                                                                                |  13.644531| 0.0000239|
|M176   |LopitDC: RIBOSOME                                                                                                                 |  12.761905| 0.0148451|
|M78    |LopitDC: ER                                                                                                                       |  12.646687| 0.0000000|
|M70    |CORUM: Spliceosome                                                                                                                |  12.612592| 0.0010329|
|M58    |LopitDC: MITOCHONDRION                                                                                                            |   9.950222| 0.0037998|
|M155   |LopitDC: PM                                                                                                                       |   9.662725| 0.0000006|
|M41    |LopitDC: MITOCHONDRION                                                                                                            |   9.567521| 0.0000000|
|M241   |CORUM: Spliceosome                                                                                                                |   9.529514| 0.0179654|
|M66    |LopitDC: RIBOSOME                                                                                                                 |   9.403509| 0.0351344|
|M1     |Takamori et al., 2006: Small GTPases and related proteins                                                                         |   9.249493| 0.0258322|
|M16    |LopitDC: CYTOSOL                                                                                                                  |   8.953343| 0.0000000|
|M81    |LopitDC: ER                                                                                                                       |   8.625390| 0.0000002|
|M48    |LopitDC: MITOCHONDRION                                                                                                            |   8.291852| 0.0000000|
|M55    |LopitDC: MITOCHONDRION                                                                                                            |   8.291852| 0.0088951|
|M56    |LopitDC: MITOCHONDRION                                                                                                            |   8.291852| 0.0088951|
|M17    |LopitDC: CYTOSOL                                                                                                                  |   7.995935| 0.0000000|
|M36    |LopitDC: CYTOSOL                                                                                                                  |   7.785515| 0.0010979|
|M52    |LopitDC: MITOCHONDRION                                                                                                            |   7.773611| 0.0031625|
|M158   |LopitDC: PM                                                                                                                       |   7.595832| 0.0015332|
|M278   |Takamori et al., 2006: All                                                                                                        |   7.134052| 0.0151220|
|M49    |LopitDC: MITOCHONDRION                                                                                                            |   7.107302| 0.0000804|
|M23    |LopitDC: CYTOSOL                                                                                                                  |   6.980117| 0.0000007|
|M131   |LopitDC: PM                                                                                                                       |   6.917632| 0.0059877|
|M34    |LopitDC: CYTOSOL                                                                                                                  |   6.812326| 0.0007864|
|M45    |LopitDC: MITOCHONDRION                                                                                                            |   6.784242| 0.0000006|
|M42    |LopitDC: MITOCHONDRION                                                                                                            |   6.633482| 0.0000000|
|M160   |LopitDC: PM                                                                                                                       |   6.588221| 0.0222836|
|M40    |LopitDC: MITOCHONDRION                                                                                                            |   6.546199| 0.0000000|
|M46    |LopitDC: MITOCHONDRION                                                                                                            |   6.546199| 0.0000144|
|M26    |LopitDC: CYTOSOL                                                                                                                  |   6.487929| 0.0000451|
|M128   |LopitDC: PM                                                                                                                       |   6.456457| 0.0039790|
|M44    |LopitDC: MITOCHONDRION                                                                                                            |   6.433333| 0.0000000|
|M25    |LopitDC: CYTOSOL                                                                                                                  |   6.228412| 0.0000608|
|M43    |LopitDC: MITOCHONDRION                                                                                                            |   6.218889| 0.0000001|
|M47    |LopitDC: MITOCHONDRION                                                                                                            |   6.218889| 0.0000804|
|M132   |LopitDC: PM                                                                                                                       |   6.149007| 0.0255644|
|M20    |LopitDC: CYTOSOL                                                                                                                  |   6.027496| 0.0000140|
|M213   |LopitDC: ER                                                                                                                       |   5.999070| 0.0052312|
|M24    |LopitDC: CYTOSOL                                                                                                                  |   5.906253| 0.0000451|
|M50    |LopitDC: MITOCHONDRION                                                                                                            |   5.330476| 0.0082747|
|M125   |LopitDC: PM                                                                                                                       |   5.270577| 0.0048037|
|M130   |LopitDC: PM                                                                                                                       |   5.124172| 0.0490439|
|M19    |LopitDC: CYTOSOL                                                                                                                  |   5.037686| 0.0001612|
|M127   |LopitDC: PM                                                                                                                       |   5.031005| 0.0255644|
|M80    |LopitDC: ER                                                                                                                       |   4.951613| 0.0062747|
|M182   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   4.315466| 0.0074460|
|M71    |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   4.027769| 0.0000531|
|M258   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   3.485569| 0.0041294|
|M260   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   3.356474| 0.0376029|
|M69    |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   3.257754| 0.0020297|
|M245   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   3.203907| 0.0305501|
|M64    |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   3.064607| 0.0006082|
|M5     |LopitDC: PM                                                                                                                       |   2.940829| 0.0287597|
|M230   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   2.517355| 0.0107404|
|M231   |LopitDC: NUCLEUS/CHROMATIN                                                                                                        |   2.407905| 0.0339811|
