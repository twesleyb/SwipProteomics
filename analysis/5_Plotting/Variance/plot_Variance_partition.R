#!/usr/bin/env Rscript

# prepare the env
root <- "~/projects/SwipProteomics"
renv::load(root)
devtools::load_all(root)

# output directory for figures
figsdir <- file.path(root,"figs","Variance")
if (!dir.exists(figsdir)) { dir.create(figsdir) }

# set plotting theme
ggtheme(); set_font("Arial",font_path=file.path(root,"fonts"))


# load the data
data(swip)
data(partition)
data(protein_gof)
data(msstats_prot) 

# imports
suppressPackageStartupMessages({
  library(dplyr)
  library(gridExtra)
  library(data.table)
  library(variancePartition)
})


# prepare the data
tidy_df <- reshape2::melt(protein_gof,id.vars=c("Protein","Symbol","Entrez"),
		     value.name="Variance", variable.name="Term") %>% 
           filter(Term %notin% c("R2.fixef","R2.total"))

# summary
summary_df <- tidy_df %>% group_by(Term) %>% 
	summarize(Mean=mean(Variance), SD = sd(Variance),.groups="drop")


# generate the plot
df <- tidy_df %>%
      filter(Term=="Genotype") %>% arrange(desc(Variance)) %>% 
      mutate(x = as.numeric(factor(Protein,levels=Protein))) %>% 
      select(Protein,x) %>% left_join(tidy_df,by="Protein")
plot <- ggplot(df)
plot <- plot + aes(x=x)
plot <- plot + aes(y=Variance)
plot <- plot + aes(group=Term)
plot <- plot + aes(stat=Variance)
plot <- plot + aes(fill=Term)
plot <- plot + geom_area()
plot <- plot + xlab("Protein")
plot <- plot + ylab("Percentage of Variance")
#plot <- plot + theme(axis.line.x=element_line())
#plot <- plot + theme(axis.line.y=element_line())
plot <- plot + theme(panel.background = element_blank())
plot <- plot + theme(axis.text.x = element_text(color="black",size=11))
plot <- plot + theme(axis.text.x = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(axis.text.y = element_text(color="black",size=11))
plot <- plot + theme(axis.text.y = element_text(angle=0,hjust=1,family="Arial"))
plot <- plot + theme(panel.border = element_rect(colour = "black", fill=NA))
plot <- plot + scale_x_continuous(expand=c(0,0))
plot <- plot + scale_y_continuous(expand=c(0,0))


# save plot
myfile <- file.path(root,"figs","Variance","Variance_partition.pdf")
ggsave(myfile, plot,height=4.5,width=4.5)

# create and save table
tt <- ttheme_default(base_size = 11, 
		     core = list(bg_params = list(fill = "white")))
tab <- tableGrob(summary_df, rows = NULL, theme = tt)
myfile <- file.path(figsdir,"Variance_table.pdf")
ggsaveTable(tab,myfile)
