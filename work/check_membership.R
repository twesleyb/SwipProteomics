#!/usr/bin/env Rscript

root = "~/projects/SwipProteomics"
renv::load(root, quiet=TRUE)
devtools::load_all(root, quiet=TRUE)


data(ne_surprise2_partition)
subpart = partition

data(ne_surprise_partition)
# partition

M21 = names(partition[partition==21])
M23 = names(partition[partition==23])
M33 = names(partition[partition==33])

message("M21:")
table(subpart[M21])

message("M23:")
table(subpart[M23])

message("M33:")
table(subpart[M33])
