#!/usr/bin/env Rscript

## ---- prepare the renv
renv::load("~/projects/SwipProteomics")

## ---- RCy3/Cytoscape info

library(RCy3)

packageVersion("RCy3")
# [1] ‘2.10.2’

RCy3::cytoscapePing() 
# You are connected to Cytoscape!

cytoscapeVersionInfo()
#      apiVersion cytoscapeVersion 
#            "v1"          "3.8.2" 


# make an example graph
g <- igraph::make_graph("Zachary")

igraph::write_graph(g, "graph.gml", format="gml")

file.exists("graph.gml")

RCy3::importNetworkFromFile("graph.gml")


myfile <- file.path(getwd(), "graph.gml")
win_file <- gsub("/mnt/d/", "D:/", myfile)

RCy3::importNetworkFromFile(win_file)

#Failed to execute: http://localhost:1234/v1/commands/network/load%20file
#Error: File 'galFiltered.sif' not found:
#Execution halted


## --- example

# download galFiltered.sif
URL <- "https://raw.githubusercontent.com/twesleyb/RCy3/master/inst/extdata/galFiltered.sif"
netw <- basename(URL)
download.file(URL,destfile=netw, quiet=TRUE) # downloaded 6822 bytes


file.exists(netw)
# [1] TRUE

importNetworkFromFile(netw)

#Failed to execute: http://localhost:1234/v1/commands/network/load%20file
#Error: File 'galFiltered.sif' not found:
#Execution halted


#I've been using `RCy3::importNetworkFromFile` to import graphs saved as `gml` files into Cytoscape. The following illustrates my problem. 
#
#```R
#library(RCy3)
#
#packageVersion("RCy3")
## [1] ‘2.10.2’
#
#RCy3::cytoscapePing() 
## You are connected to Cytoscape!
#
#cytoscapeVersionInfo()
##      apiVersion cytoscapeVersion 
##            "v1"          "3.8.2" 
#
## create an example graph
#g <- igraph::make_graph("Zachary")
#igraph::write_graph(g, "graph.gml", format="gml")
#
#file.exists("graph.gml")
## TRUE
#
## load into Cytoscape
#RCy3::importNetworkFromFile("graph.gml")
#```
#I encounter the error:
#> Failed to execute: http://localhost:1234/v1/commands/network/load%20file
#>Error: File 'graph.gml' not found:
#>Execution halted
#
#The following also reproduces my problem with the example `galFiltered.sif` network.
#
#### Reproducible Example
#
#```R
## download 'galFiltered.sif'
#URL <- "https://raw.githubusercontent.com/twesleyb/RCy3/master/inst/extdata/galFiltered.sif"
#netw <- basename(URL)
#
#download.file(URL, destfile=netw, quiet=TRUE)
#
#file.exists(netw)
## [1] TRUE
#
#importNetworkFromFile(netw)
#```
#> Failed to execute: http://localhost:1234/v1/commands/network/load%20file
#> Error: File 'galFiltered.sif' not found:
#> Execution halted
#
#This method worked about a month ago, I've updated RCy3 and Cytoscape, but I have not tried going back in time to see if a previous version of RCy3 works.
#
#I'm working on the Window Subsytem for Linux. Although Cytoscape is installed on the Windows side of things and I am using R on the linux side of things, previously, this didn't seem to matter much because--I'm guessing because the communication between the API intermediate. I only encountered problems when trying to reference paths with Cytoscape. Since its a Windows program you need to make sure that the file-path is a Windows file path. The following worked for me. 
#
#```R
#myfile
#[1] "/mnt/d/projects/SwipProteomics/graph.gml"
#
#win_file <- gsub("/mnt/d/", "D:/", myfile)
#
#win_file
#[1] "D:/projects/SwipProteomics/graph.gml"
#```
