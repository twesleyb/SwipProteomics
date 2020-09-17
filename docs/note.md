### `tidy_peptide` - the Swip TMT dataset, raw peptides

 | Name                       | Class| Values |
 |----------------------------|-----|---------------------------------|
 | Confidence                 | chr |"High" "High" "High" "High" ...        |
 | Accession                  | chr |"A1L317" "A2A432" "A2A432" "A2A432" ...  |
 | Description                | chr |"Keratin, type I cytoskeletal 24" ...  |
 | Sequence                   | chr |"QTMQDLNDR" "AQIHQFREDSLDSVLFLK" "DIMIQFK" ...  |
 | Modifications              | chr |"1xOxidation [M3]; 1xTMTpro [N-Term]" ... |
 | Retention.Time             | num |16.9 76.9 70.6 78.7 19.5 ... |
 | Ion.Score                  | int |12 7 30 34 18 55 32 22 42 27 ... |
 | q-Value                    | num |5.91e-03 8.84e-06 1.34e-03 1.09e-04 5.08e-04 |
 | Sample                     | chr |Abundance.F1.126.Sample.53395.Control4.Exp1" |
 | Intensity                  | num |NA 71.6 36.2 254 NA ... |
 | Date of Brain Fractionation| int |11062019 11062019 11062019 11062019 ... |
 | Proteomics ID              | chr |"Exp1Control4" "Exp1Control4" ...  |
 | Mouse ID                   | int |23733 23733 23733 23733 23733 23733 23733 ...    |
 | Genotype                   | chr |"WT" "WT" "WT" "WT" ... |
 | Age (mo)                   | int |10 10 10 10 10 10 10 10 10 10 ... |
 | DOB                        | int |1062019 1062019 1062019 1062019 ... |
 | Sex                        | chr |"Female" "Female" "Female" "Female" ... |
 | Treatment                  | chr |"Control" "Control" "Control" "Control" ...  |
 | Experiment                 | chr |"Exp1" "Exp1" "Exp1" "Exp1" ... |
 | Channel                    | chr |"126" "126" "126" "126" ... |
 | Fraction                   | chr |"F4" "F4" "F4" "F4" ... |
 | Cfg Force (xg)             | chr |"5,000" "5,000" "5,000" "5,000" ... |


### `input.pd` the dataset included with the MSstatsTMT package

 | Name                       | Class| Values |
 |----------------------------|-----|---------------------------------|
 | ProteinName    | Factor|w/ 10 levels "P04406","P06576" ... |
 | PeptideSequence| Factor|w/ 228 levels "[K].aAGASDVVLYk.[I]" ... |
 | Charge         | int   | 2 3 3 3 3 3 3 2 3 3 ... |
 | PSM            | Factor|w/ 305 levels "[K].aAGASDVVLYk.[I]_2", ... |
 | Mixture        | Factor|w/ 5 levels "Mixture1","Mixture2",..: 1 1 1 1 1 1 1 1 1 1 ... |
 | TechRepMixture | Factor|w/ 3 levels "1","2","3": 1 1 1 1 1 1 1 1 1 1 ... |
 | Run            | Factor|w/ 15 levels "161117_SILAC_HeLa_UPS1_TMT10_Mixture1_01.raw" ... |
 | Channel        | Factor|w/ 10 levels "126","127C","127N", ... |
 | Condition      | Factor|w/ 5 levels "0.125","0.5",..: 5 5 5 5 5 5 5 5 ... |
 | BioReplicate   | Factor|w/ 25 levels "Mixture1_0.125",..: 5 5 5 5 5 5 ... |
 | Intensity      | num   |8348 28327 1275011 80590 2231 ... |
