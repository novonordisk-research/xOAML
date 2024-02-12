# Code to prepare input data to generate polygenic risk scores (PRS) at the gene and pathway levels
# Nielsen et al, Nat Comms 2024

################# Load libraries ###############
library(msigdf)
library(dplyr)
library(rtracklayer)
library(tibble)
library(data.table)

################# Load annotation ###############

# Data from https://www.gencodegenes.org/human/release_43lift37.html
gtf <- rtracklayer::import("input_data_path/gencode.v43lift37.annotation.gtf")
gtf_df=as.data.frame(gtf)
head(gtf_df)

################# Extract protein-coding genes names ###############

coding_genes <- gtf_df %>%
  filter(transcript_type == "protein_coding") %>%
  select(gene_name) %>%
  distinct()  %>%
  mutate(pathway_data = paste(paste0("GRS_",gene_name), gene_name)) %>%
  select(pathway_data) 
head(coding_genes)

################# Save file ###############

write.table(coding_genes,
            "input_data_path/coding_genes.txt", 
            col.names = F, 
            row.names = F, 
            quote=F, 
            sep=" ")


################# Generate a list of genes per pathways ###############

# Pathway data from https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2022.1.Hs/

list_msigdb <- list()

setnames_list <- c("c2.cp.reactome",
                   "c2.cp.kegg",
                   "c2.cp.biocarta",
                   "c3.tft",
                   "c5.go.bp",
                   "h.all")

# Data from 
for (set_names in setnames_list) {
  
  msigdb_data <- GSA::GSA.read.gmt(paste0("input_data_path/",set_names,".v2022.1.Hs.symbols.gmt"))
  
  list_msigdb[[set_names]] <- cbind(msigdb_data$geneset.names,msigdb_data$genesets) %>% 
    as.data.frame() %>%
    magrittr::set_colnames(c("geneset","symbol")) %>%
    mutate(symbol = lapply(symbol,paste,collapse = " ")) %>%
    mutate(set_name = set_names) %>%
    dplyr::select(set_name, everything()) %>%
    mutate(pathway_data = paste(paste0("PRS_",set_name,"_",geneset), symbol))
  
  
}

msigdb_all <- do.call(rbind,list_msigdb)
head(msigdb_all)
table(msigdb_all$set_name)
dim(msigdb_all)


################# Save file per pathway set ###############

for (set_names in setnames_list){
  
  msigdb_all %>%
    filter(set_name == set_names) %>%
    select(pathway_data) %>%
    write.table(.,
                paste0("input_data_path/",set_names,"_for_PRS.txt"), 
                col.names = F, 
                row.names = F, 
                quote=F, 
                sep=" ")
  
}

################# Save file including all selected pathway sets ###############

write.table(msigdb_all$pathway_data,
            "input_data_path/msigdb_subset_for_PRS.txt", 
            col.names = F, 
            row.names = F, 
            quote=F, 
            sep=" ")



