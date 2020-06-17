#' createCytoscapeGraph
createCytoscapeGraph <- function(netw_g, ppi_g, nodes, module_name,
				 n_cutoffs=5000, netwdir=getwd(),
				 imgsdir=getwd(), netw_layout=NULL) {

	# Imports
	suppressPackageStartupMessages({
		library(RCy3)
	})

	# Default layout.
	if (is.null(netw_layout)) { 
		netw_layout <- 'force-directed edgeAttribute=weight' 
	}

	# Check that we are connected to Cytoscape.
	response <- cytoscapePing()
	if (response != "You are connected to Cytoscape!") { 
		stop("Start Cytoscape first!") 
	}

	# Subset graph: keep nodes in module.
	graph <- netw_g
	idx <- match(nodes, names(V(graph)))
	subg <- induced_subgraph(graph, vids = V(graph)[idx])

	# Prune weak edges.
	# Seq from min(edge.weight) to max to generate cutoffs.
	n_edges <- length(E(subg))
	min_weight <- min(E(subg)$weight)
	max_weight <- max(E(subg)$weight)
	cutoffs <- seq(min_weight, max_weight, length.out = n_cutoffs)

	# Define a function that checks if graph is connnected at a 
	# given edge weight threshold.
	is_connected <- function(graph,threshold) {
		filt_graph <- delete.edges(graph, 
				       which(E(graph)$weight <= threshold))
		return(is.connected(filt_graph))
	}

	# Check if graph is connected or not at various thresholds.
	message("Thresholding graph...")
	checks <- sapply(cutoffs, function(threshold) {
				 is_connected(subg,threshold)
				       })

	# Limit is max(cutoff) at which the graph is still connected.
	limit <- cutoffs[max(which(checks==TRUE))]
	if (all(checks)) { stop("Error thesholding graph.") }

	# Prune edges. NOTE: This removes all edge types.
	g <- delete.edges(subg, which(E(subg)$weight <= limit))
	n_edges <- length(E(g))

	# Write graph to file this is faster than sending to cytoscape.
	myfile <- file.path(netwdir, paste0(module_name, ".gml"))
	write_graph(g, myfile, format = "gml")

	# Send to Cytoscape.
	## FIXME: underscores from edge weight attributes are removed!
	winfile <- gsub("/mnt/d/", "D:/", myfile)
	cysnetw <- importNetworkFromFile(winfile)
	Sys.sleep(2); unlink(myfile)

	####################################################################
	## VISUAL STYLE DEFAULTS:
	####################################################################

	# Create a visual style.
	style.name <- paste(module_name, "style", sep = "-")

	# VISUAL STYLE DEFAULTS:
	defaults <- list(
		 NODE_FILL_COLOR = col2hex("gray"),
		 NODE_TRANSPARENCY = 200,
		 NODE_SIZE = 35,
		 NODE_SHAPE = "ellipse",
		 NODE_LABEL_TRANSPARENCY = 255,
		 NODE_LABEL_FONT_SIZE = 12,
		 NODE_LABEL_COLOR = col2hex("black"),
		 NODE_BORDER_TRANSPARENCY = 200,
		 NODE_BORDER_WIDTH = 4,
		 NODE_BORDER_PAINT = col2hex("black"),
		 NODE_TRANSPARENCY = 200,
		 EDGE_STROKE_UNSELECTED_PAINT = col2hex("black"),
		 EDGE_WIDTH = 2,
		 NETWORK_BACKGROUND_PAINT = col2hex("white")
		)

	####################################################################
	## VISUAL STYLE MAPPINGS:
	####################################################################

# MAPPED PROPERTIES:
min_max_weight <- c(min(E(g)$weight), max(E(g)$weight))
edge_colors <- c(col2hex("gray"), col2hex("dark red"))
mappings <- list(
		 NODE_FILL_COLOR = mapVisualProperty("node fill color","Color","p"),
		 NODE_LABEL = mapVisualProperty("node label", 
							"symbol", "p"),
EDGE_TRANSPARENCY = mapVisualProperty("edge transparency", "weight", 
				      "c", min_max_weight, c(155, 255)),
EDGE_STROKE_UNSELECTED_PAINT = mapVisualProperty("edge stroke unselected paint",
						 "weight", "c",min_max_weight,
						 edge_colors) )

	# Create a visual style.
	createVisualStyle(style.name, defaults = defaults, 
			  mappings = mappings)
	# Apply to graph.
	setVisualStyle(style.name)
	Sys.sleep(3)
	# Collect PPI edges.
	idx <- match(nodes, names(V(ppi_g)))
		subg <- induced_subgraph(ppi_g, 
					 vids = V(ppi_g)[idx])
	edge_list <- apply(as_edgelist(subg, names = TRUE), 1, as.list)
	# If edge list is only of length 1, unnest it to avoid problems.
	if (length(edge_list) == 1) {
		edge_list <- unlist(edge_list, recursive = FALSE)
	}
	# Add PPI edges to Cytoscape graph.
	if (length(edge_list) > 0) {
		ppi_edges <- addCyEdges(edge_list)
		# Add PPIs and set to black.
		selected_edges <- selectEdges(ppi_edges, by.col = "SUID")
		# Set to black with edge bend.
		namen <- "EDGE_STROKE_UNSELECTED_PAINT"
		setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = col2hex("black"), 
				      visual.property = namen,
				      bypass = TRUE
				      )
		setEdgePropertyBypass(
				      edge.names = selected_edges$edges,
				      new.values = TRUE,
				      visual.property = "EDGE_BEND",
				      bypass = TRUE
				      )
	} # Ends IF statement.
	#  Clean-up.
	clearSelection()
	Sys.sleep(2) 
	# Apply layout.
	layoutNetwork(netw_layout)
	Sys.sleep(2)
	fitContent()
	# Save Image.
	netw_image <- file.path(imgsdir, module_name)
	winfile <- gsub("/mnt/d/", "D:/", netw_image)
	exportImage(winfile, "svg")
	# Free up some memory.
	cytoscapeFreeMemory()
}
