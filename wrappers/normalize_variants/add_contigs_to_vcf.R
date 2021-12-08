library(data.table)

args <- commandArgs(trailingOnly = T)

add_contigs_to_vcf <- function(vcf_file,dict_file){
    lines <- readLines(vcf_file)
    header_non_contig_idx <- grep("^##[^c][^o][^n][^t][^i]",lines)
    non_header_idx <- grep("^##",lines,invert = T)
    
    contig_tab <- fread(dict_file,skip = "@SQ",header = F)
    
    contig_lines <- paste0("##contig=<ID=",gsub("SN\\:","",contig_tab$V2),",length=",gsub("LN\\:","",contig_tab$V3),">")
    lines <- paste(c(lines[header_non_contig_idx],contig_lines,lines[non_header_idx]),collapse = "\n")
    cat(lines,"\n",file = vcf_file)
    
}

add_contigs_to_vcf(args[1],args[2])