library(data.table)


run_all <- function(args){
  output_file <- args[1]
  input_files <- args[2:length(args)]
  
  tab_to_print <- rbindlist(lapply(input_files,fread))
  
  ###### extract chrom,position,reference,alternative from varname


  
  tab_to_print[,"chromosome"]<-sapply(strsplit(as.character(tab_to_print$var_name), "_"), "[[", 1)
  tab_to_print[,"start"]<-sapply(strsplit(as.character(tab_to_print$var_name), "_"), "[[", 2)
  splitted<-sapply(strsplit(as.character(tab_to_print$var_name), "_"), "[[", 3)
  tab_to_print[,"allele"]<- splitted
  tab_to_print[,"reference"]<-sapply(strsplit(as.character(splitted), "/"), "[[", 1)
  #tab_to_print[,"alternative"]<-sapply(strsplit(as.character(splitted), "/"), "[[", 2)
  tab_to_print[,start := as.integer(start)]
  tab_to_print[,end := as.integer(start) + nchar(reference) - 1L]
  tab_to_print[,strand := "+"]

  # VEP > WARNING: Alleles look like an insertion (-/T) but coordinates are not start = end + 1, PROTO:
  tab_to_print[,reference := gsub("/.*","",allele)]
  tab_to_print[reference == "-",end := start - 1L]
  
  tab_to_print <- tab_to_print[,.(chromosome,start,end,allele,strand)]

  setkey(tab_to_print)
  
  fwrite(tab_to_print,file = output_file,col.names = F,row.names = F,sep = "\t",quote = F)
}


# develop and test
#args <- c("annotate/all_variants.tsv","somatic_seq_results/Pca_4.variants.tsv","somatic_seq_results/Pca_5.variants.tsv","somatic_seq_results/Pca_6.variants.tsv","somatic_seq_results/Pca_7.variants.tsv","somatic_seq_results/Pca_8.variants.tsv","somatic_seq_results/Pca_10.variants.tsv","somatic_seq_results/Pca_11.variants.tsv","somatic_seq_results/Pca_12.variants.tsv","somatic_seq_results/Pca_13.variants.tsv","somatic_seq_results/Pca_17.variants.tsv","somatic_seq_results/Pca_18.variants.tsv","somatic_seq_results/Pca_19.variants.tsv","somatic_seq_results/Pca_20.variants.tsv","somatic_seq_results/Pca_21.variants.tsv","somatic_seq_results/Pca_22.variants.tsv","somatic_seq_results/Pca_23.variants.tsv","somatic_seq_results/Pca_24.variants.tsv","somatic_seq_results/Pca_25.variants.tsv","somatic_seq_results/Pca_26.variants.tsv")
#setwd("/mnt/ssd/ssd_1/snakemake/stage359_PC.seq_A/somatic_variant_calling")

#develop and test rj WES 86
#args <- c("annotate/all_variants.tsv", "somatic_seq_results/CK1151.variants.tsv", "somatic_seq_results/JD2011.variants.tsv", "somatic_seq_results/KM1755.variants.tsv")
#setwd("/mnt/ssd/ssd_1/snakemake/stage362_solid_tumors_children.WES_86/somatic_variant_calling")


#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)