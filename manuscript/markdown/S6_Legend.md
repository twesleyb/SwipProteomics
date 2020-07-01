# Supplementary File 6

Cytoscape subnetworks for all __255__ modules, arranged in numerical order. See
also Supplementary File 9 for the Cytoscape file containing all module subnetworks.

Vizualization of covariation graphs is challenging because every pair of nodes is connected.
In order to resolve the underlying structure of each module, a dynamic threshold was 
set for each module to remove weak edges. The edge weight threshold was set as the maximum edge weight at which the graph is a single connected component.
Edges weaker than the threshold were removed. This strategy enables visualization of the strongest paths in a graph.

#### Example Plot:
![plot](../../figs/Proteins/S6_Example.png)

#### Network Attributes:
| Attribute | Description |
| --------- | ----------- |
| _Node Color_  | Significant proteins¹ were assigned their modules' color. |
| _Node Size_   | Node weighted degree centrality. Node size indicates the importance of a node in its module. |
| _Edge Color_  | Edge color cooresponds to the strength of the correlation between two proteins (gray=weak, dark-red=strong). Black edges indicate PPIs compiled from the HitPredict database. |
| _Layout_ | Cytoscape's perfuse force-directed layout; weight = edge weight. |

__¹__ Significant Proteins: WT v MUT FDR < 0.05 (n=686; see Supplementary File 5).
