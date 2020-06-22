#!/usr/bin/env Rscript

cmd <- c("curl http://colormind.io/api/ --data-binary","{'model':'default'}")
cmd <- paste(cmd,collapse=" ")

curl 'http://colormind.io/api/' --data-binary '{"model":"default"}'

curl http://colormind.io/api/ --data-binary '{"model":"default"}'
					
# {"result":[[214,78,69],[247,242,163],[201,216,147],[57,141,112],[62,80,64]]}
system(cmd, intern = TRUE)

curl 'http://colormind.io/api/' --data-binary '{"input":[[44,43,44],[90,83,82],"N","N","N"],"model":"default"}'
					
