#!/usr/bin/env Rscript

root <- "~/projects/SwipProteomics"
devtools::load_all(root)

data(alt_partition)
alt_part <- partition

data(swip_partition)
swip_part <- partition

data(msstats_partition)
msstats_part <- partition

prots <- Reduce(intersect, list(names(swip_part),names(alt_part),names(msstats_part)))

length(prots)

calcPartSim(swip_part[prots], alt_part[prots])

calcPartSim(swip_part[prots], msstats_part[prots])
