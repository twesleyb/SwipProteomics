#!/usr/bin/env Rscript

#' ### 2.2 dataProcessPlotsTMT()	
#' 	
#' Visualization for explanatory data analysis. To illustrate the quantitative data after data-preprocessing and quality control of TMT runs, dataProcessPlotsTMT takes the quantitative data from converter functions (`PDtoMSstatsTMTFormat`, `MQtoMSstatsTMTFormat` and `SpectroMinetoMSstatsTMTFormat`) and summarized data from function `proteinSummarization` as input. It generates two types of figures in pdf files as output :	
#' 	
#' (1) profile plot (specify "ProfilePlot" in option type), to identify the potential sources of variation for each protein;	
#' 	
#' (2) quality control plot (specify "QCPlot" in option type), to evaluate the systematic bias between MS runs.	
#' 	
#' 	
#' #### Arguments	
#' 	
#' * `data.peptide` : name of the data with peptide-level, which can be the output of converter functions (`PDtoMSstatsTMTFormat`, `MQtoMSstatsTMTFormat` and `SpectroMinetoMSstatsTMTFormat`).	
#' * `data.summarization` : name of the data with protein-level, which can be the output of `proteinSummarization` function.	
#' * `type` : choice of visualization. "ProfilePlot" represents profile plot of log intensities across MS runs.	
#' "QCPlot" represents quality control plot of log intensities across MS runs.	
#' * `ylimUp` : upper limit for y-axis in the log scale.	
#' FALSE(Default) for Profile Plot and QC Plot use the upper limit as rounded off maximum of log2(intensities) after normalization + 3.	
#' * `ylimDown` : lower limit for y-axis in the log scale. FALSE(Default) for Profile Plot and QC Plot is 0.	
#' * `x.axis.size` : size of x-axis labeling for "Run" and "channel" in Profile Plot and QC Plot.	
#' * `y.axis.size` : size of y-axis labels. Default is 10.	
#' * `text.size` : size of labels represented each condition at the top of graph in Profile Plot and QC plot. Default is 4.	
#' * `text.angle` : angle of labels represented each condition at the top of graph in Profile Plot and QC plot. Default is 0.	
#' * `legend.size` : size of legend above graph in Profile Plot. Default is 7.	
#' * `dot.size.profile` : size of dots in profile plot. Default is 2.	
#' * `ncol.guide` : number of columns for legends at the top of plot. Default is 5.	
#' * `width` : width of the saved file. Default is 10.	
#' * `height` : height of the saved file. Default is 10.	
#' * `which.Protein` : Protein list to draw plots. List can be names of Proteins or order numbers of Proteins.	
#' Default is "all", which generates all plots for each protein. For QC plot, "allonly" will generate one QC plot with all proteins.	
#' * `originalPlot` : TRUE(default) draws original profile plots, without normalization.	
#' * `summaryPlot` : TRUE(default) draws profile plots with protein summarization for each channel and MS run.	
#' * `address` : the name of folder that will store the results. Default folder is the current working directory.	
#' The other assigned folder has to be existed under the current working directory.	
#' An output pdf file is automatically created with the default name of "ProfilePlot.pdf" or "QCplot.pdf".	
#' The command address can help to specify where to store the file as well as how to modify the beginning of the file name.	
#' If address=FALSE, plot will be not saved as pdf file but showed in window.	
#' 	
#' #### Example	
#' 	
## Profile plot without norm channnels and empty channels	
dataProcessPlotsTMT(data.peptide = input.pd,	
                     data.summarization = quant.msstats,	
                     type = 'ProfilePlot',	
                     width = 21, # adjust the figure width since there are 15 TMT runs.	
                     height = 7)	
	
# ## Profile plot with all the channels	
# quant.msstats.all <- proteinSummarization(input.pd,	
#                                       method="msstats",	
#                                       normalization=TRUE,	
#                                       remove_norm_channel=FALSE,	
#                                       remove_empty_channel=FALSE)	
# 	
# dataProcessPlotsTMT(data.peptide = input.pd,	
#                      data.summarization = quant.msstats.all,	
#                      type = 'ProfilePlot',	
#                      width = 21, # adjust the figure width since there are 15 TMT runs.	
#                      height = 7)	
	
## Quality control plot 	
# dataProcessPlotsTMT(data.peptide=input.pd,	
                     # data.summarization=quant.msstats, 	
                     # type='QCPlot',	
                     # width = 21, # adjust the figure width since there are 15 TMT runs. 	
                     # height = 7)	
