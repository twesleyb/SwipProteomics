Loading required package: Matrix

Attaching package: ‘lmerTest’

The following object is masked from ‘package:lme4’:

    lmer

The following object is masked from ‘package:stats’:

    step

Loading SwipProteomics


|Coefficient          | Estimate| Std. Error|       df|  t value|p value   |
|:--------------------|--------:|----------:|--------:|--------:|:---------|
|ConditionControl.F10 | 7.618697|   0.121126| 17.30594| 62.89894|7.045e-22 |
|ConditionControl.F4  | 6.710959|   0.121126| 17.30594| 55.40477|6.264e-21 |
|ConditionControl.F5  | 6.945583|   0.121126| 17.30594| 57.34180|3.467e-21 |
|ConditionControl.F6  | 7.240081|   0.121126| 17.30594| 59.77313|1.696e-21 |
|ConditionControl.F7  | 7.321074|   0.121126| 17.30594| 60.44180|1.4e-21   |
|ConditionControl.F8  | 7.129632|   0.121126| 17.30594| 58.86129|2.21e-21  |
|ConditionControl.F9  | 6.954472|   0.121126| 17.30594| 57.41518|3.391e-21 |
|ConditionMutant.F10  | 5.784403|   0.121126| 17.30594| 47.75525|8.066e-20 |
|ConditionMutant.F4   | 5.404300|   0.121126| 17.30594| 44.61718|2.592e-19 |
|ConditionMutant.F5   | 5.567441|   0.121126| 17.30594| 45.96405|1.555e-19 |
|ConditionMutant.F6   | 5.640188|   0.121126| 17.30594| 46.56463|1.245e-19 |
|ConditionMutant.F7   | 5.631680|   0.121126| 17.30594| 46.49440|1.277e-19 |
|ConditionMutant.F8   | 5.492772|   0.121126| 17.30594| 45.34759|1.961e-19 |
|ConditionMutant.F9   | 5.781022|   0.121126| 17.30594| 47.72734|8.148e-20 |
Abundance ~ 0 + Condition + (1|Mixture)


|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|Mutant-Control | -1.516956|      0.3494225| 0.0572533|  -26.49552|      0| 26|FALSE      |


|Coefficient         |   Estimate| Std. Error|       df|    t value|p value   |
|:-------------------|----------:|----------:|--------:|----------:|:---------|
|(Intercept)         |  7.6186967|  0.1211260| 17.30594|  62.898935|7.045e-22 |
|ConditionControl.F4 | -0.9077381|  0.1514779| 26.00000|  -5.992543|2.51e-06  |
|ConditionControl.F5 | -0.6731139|  0.1514779| 26.00000|  -4.443643|0.000146  |
|ConditionControl.F6 | -0.3786161|  0.1514779| 26.00000|  -2.499480|0.01908   |
|ConditionControl.F7 | -0.2976228|  0.1514779| 26.00000|  -1.964793|0.06021   |
|ConditionControl.F8 | -0.4890643|  0.1514779| 26.00000|  -3.228617|0.003356  |
|ConditionControl.F9 | -0.6642249|  0.1514779| 26.00000|  -4.384961|0.0001704 |
|ConditionMutant.F10 | -1.8342939|  0.1514779| 26.00000| -12.109313|3.427e-12 |
|ConditionMutant.F4  | -2.2143966|  0.1514779| 26.00000| -14.618607|4.701e-14 |
|ConditionMutant.F5  | -2.0512553|  0.1514779| 26.00000| -13.541610|2.751e-13 |
|ConditionMutant.F6  | -1.9785088|  0.1514779| 26.00000| -13.061366|6.264e-13 |
|ConditionMutant.F7  | -1.9870162|  0.1514779| 26.00000| -13.117528|5.683e-13 |
|ConditionMutant.F8  | -2.1259246|  0.1514779| 26.00000| -14.034549|1.21e-13  |
|ConditionMutant.F9  | -1.8376747|  0.1514779| 26.00000| -12.131632|3.29e-12  |
Abundance ~ 1 + Condition + (1|Mixture)


|Contrast       |    log2FC| percentControl|        SE| Tstatistic| Pvalue| DF|isSingular |
|:--------------|---------:|--------------:|---------:|----------:|------:|--:|:----------|
|Mutant-Control | -1.516956|      0.3494225| 0.0572533|  -26.49552|      0| 26|FALSE      |
