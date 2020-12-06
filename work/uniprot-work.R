#!/usr/bin/env Rscript

## ---- input 
root <- "~/projects/SwipProteomics"

# downloaded on 12/06/2020
zipfile <-  file.path(root,"downloads","uniprot-download.tab.gz")

## ---- prepare the renv

renv::load(root, quiet=TRUE)

# library(SwipProteomics)
devtools::load_all(root, quiet=TRUE)


## ---- load the data

data(swip)
data(gene_map)
data(washc_prots)
data(wash_interactome)
data(ne_surprise_partition)

## load uniprot 'Subcellular location [CC]' annotations
# parse the uniprot data

parse_uniprot <- function(zipfile) {
  df <- data.table::fread(zipfile)
  colnames(df)[1] <- "Protein"
  colnames(df)[colnames(df) == "Subcellular location [CC]"] <- "Subcell_loc"
  df = df %>% dplyr::select(Protein,"Subcell_loc") 
  df = df %>% tidyr::separate_rows(Subcell_loc,sep="\\. ") # split rows (annotations) delimited by '.'
  df = df %>% subset(Subcell_loc != "") # remove empty
  df = df %>% dplyr::mutate(Subcell_loc = gsub("SUBCELLULAR LOCATION: ", "",Subcell_loc))
  df = df %>% dplyr::mutate(Subcell_loc = gsub("\\.", "",Subcell_loc)) # remove '.'
  df = df %>% subset(!grepl("Note=",Subcell_loc)) # drop notes (NOTE: rows with just text remain)
  df = df %>% subset(grepl("\\{*\\}",Subcell_loc)) # keep evidence
  df = df %>% tidyr::separate_rows(Subcell_loc,sep="\\; ") # split rows (annotations) delimited by ';'
  # this isnt even for R!
  # https://howtodoinjava.com/java/regex/start-end-of-string/
  #grepl("^[\\{]",x) # match start with special
  #grepl("[\\}]$",x) # match end with special
  df = df %>% subset(!grepl("^[\\{]", Subcell_loc)) # drop rows that are just evidence (from end of Note= usually)
  df = df %>% subset(grepl("\\{*\\}", Subcell_loc)) # drop rows without evidence
  # final push
  mylist <- lapply(df$Subcell_loc, function(x) setNames(gsub("\\}","", strsplit(x," \\{")[[1]]),nm=c("Annotation","Evidence")))
  names(mylist) <- df$Protein
  anno_df = dplyr::bind_rows(mylist, .id="Protein")
  return(anno_df)
} #EOF


subcell_df  <- parse_uniprot(zipfile)

## ---- save

uniprot_subcell <- subcell_df
myfile <- "uniprot_subcell.rda"
save(uniprot_subcell, file=myfile,version=2)

head(subcell_df)
nrow(subcell_df)

# annot for perccent part
sum(names(partition) %in% unique(subcell_df$Protein))/length(partition)

df = subcell_df %>% filter(grepl("Endosome|endosome",Annotation))
endosome <- unique(df$Protein)

myfile <- "endosome.rda"
save(endosome,file=myfile,version=2)

endosome %in% names(partition)

data(ne_surprise_partition)
m23 <- names(which(partition==23))
sum(m23 %in% endosome)


data(ne_surprise2_partition)
m501 = names(which(partition==501))

washc_prots %in% endosome

sum(m501 %in% endosome)

