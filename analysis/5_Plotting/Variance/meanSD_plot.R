#!/usr/bin/env Rscript


renv::load("~/projects/SynaptopathyProteomics")

root <- "~/projects/SwipProteomics"
devtools::load_all(root)

data(msstats_prot)

myfile <- file.path(root,"rdata","msstats_psm.rda")
load(myfile)

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# set plotting theme and font
ggtheme(); set_font("Arial",file.path(root,"figs"))

# cast the data into a matrix
raw_dm <- msstats_psm %>% reshape2::dcast(PSM ~ Mixture + Channel + Condition, 
				       value.var = "Intensity") %>%
                       as.data.table() %>% as.matrix(rownames="PSM")

norm_dm <- msstats_prot %>% reshape2::dcast(Protein ~ Mixture + Condition, 
				       value.var = "norm_Abundance") %>%
                       as.data.table() %>% as.matrix(rownames="Protein")

# generate a mean-variance plot
plot_data <- vsn::meanSdPlot(raw_dm)
plot <- plot_data$gg
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(colour = "black", fill=NA))
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))

plot_data <- vsn::meanSdPlot(norm_dm)
plot <- plot_data$gg
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(panel.border = element_rect(colour = "black", fill=NA))
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))

# save
myfile <- file.path(root,"figs","Variance","meanSD.pdf")
ggsave(myfile,plot)
