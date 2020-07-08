rgb2col <- function(r,g,b) {

       	#create colour name vs. rgb mapping table 
	colourMap = data.frame(colourNames = colours(),t(col2rgb(colours())))

	#input test colours
	r = 155
	g 
	testDF = data.frame(colourNames="testCol",red = r,green = g,blue = b)

	#combine both tables
	combDF = rbind(testDF,colourMap)

	#convert in matrix representation 
	combMat= (as.matrix(combDF[,-1]))

	#add row labels as colour names
	rownames(combMat) = combDF[,1]

	#compute euclidean distance between test rgb vector and all the colours
	#from mapping table 
	#using dist function compute distance matrix,retain only upper matrix
	#find minimum distance point from test vector

	#find closest matching colour name
	approxMatchCol = which.min(as.matrix(dist(combMat,upper=TRUE))[1,][-1])

	#compare test colour with approximate matching colour
	scales::show_col(c(rgb2hex(r,g,b),rgb2hex(colourMap[approxMatchCol,2:4])))

	#return colour name
	return(approxMatchCol)

}
