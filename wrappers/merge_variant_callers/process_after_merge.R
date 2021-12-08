suppressMessages(library(data.table))
suppressMessages(library(VennDiagram))

run_all <- function(args){
  var_file <- args[1]
  output_file <- args[2]
  variant_type <- args[3]

  vcf <- vcfR::read.vcfR(var_file,verbose = F)
  
  if(variant_type == "germline"){
    var_tab <- data.table(chrom = vcf@fix[,"CHROM"]
                          ,position = as.integer(vcf@fix[,"POS"])
                          ,reference = vcf@fix[,"REF"]
                          ,alternative = vcf@fix[,"ALT"]
                          ,filter = vcf@fix[,"FILTER"]
                          ,genotype = vcfR::extract.gt(vcf)[,1]
                          ,variant_freq = vcfR::extract.gt(vcf,element = "AD",as.numeric = T)[,1]
                          ,coverage_depth = vcfR::extract.gt(vcf,element = "DP",as.numeric = T)[,1])

    var_tab[is.na(coverage_depth),coverage_depth := vcfR::extract.gt(vcf,element = "DPI",as.numeric = T)[,1][!is.na(vcfR::extract.gt(vcf,element = "DPI",as.numeric = T)[,1])]]

    var_tab[,index := seq_along(var_tab$chrom)]
    var_tab[,variant_freq := round((coverage_depth - variant_freq) / coverage_depth,5)]
    
    # ADD INFO ABOUT MULTIPLE CALLERS and GET UNIQUE VARIANTS
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "GQ")[,1]),caller := "haplotypecaller"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "ALD")[,1]),caller := "vardict"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "GQX")[,1]),caller := "strelka"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "ABQ")[,1]),caller := "varscan"]
    var_tab[caller == "haplotypecaller",filter := "PASS"]
    var_tab[caller == "varscan",variant_freq := 1 - variant_freq]
    
    if(any(var_tab$caller == "varscan",na.rm = T)){
      ADF <- vcfR::extract.gt(vcf,element = "ADF",as.numeric = T)[var_tab$caller == "varscan"]
      RDF <- vcfR::extract.gt(vcf,element = "RDF",as.numeric = T)[var_tab$caller == "varscan"]
      ADR <- vcfR::extract.gt(vcf,element = "ADR",as.numeric = T)[var_tab$caller == "varscan"]
      RDR <- vcfR::extract.gt(vcf,element = "RDR",as.numeric = T)[var_tab$caller == "varscan"]
      
      var_tab[var_tab$caller == "varscan",strand_bias := boot::inv.logit(ADF * (ADR + RDR) / (ADR * (ADF + RDF)))] 
                
    } else {
      var_tab[,strand_bias := NA]
    }
    
    var_tab <- var_tab[alternative != "*"]

    #my filter
    var_tab <- var_tab[variant_freq * coverage_depth > 4]

    # GET CALLING STATISTICS
    stat_tab <- copy(var_tab)
    
    #FILTER CALLS
    var_tab <- var_tab[filter == "PASS"]
    var_tab <- var_tab[coverage_depth > 0]
    var_tab[variant_freq > 1,variant_freq := 1]
    
    if(nrow(var_tab) > 0){
      # GET CALLING STATISTICS
      stat_tab[,is_pass := filter == "PASS"]
      stat_tab[,by_caller_sum := .N,by = c("caller")]
      
      stat_tab_out1 <- stat_tab[,list(variants = .N),by = c("caller")]
      stat_tab_out1[,c("is_pass","percent_of_caller") := NA]
      stat_tab_out2 <- stat_tab[,list(variants = .N,percent_of_caller = round(.N / by_caller_sum[1] * 100,1)),by = c("caller","is_pass")]
      setcolorder(stat_tab_out1,names(stat_tab_out2))
      setorder(stat_tab_out1)
      setorder(stat_tab_out2)
      
      write.table(rbind(stat_tab_out1,stat_tab_out2),file = gsub(".tsv$",".stats_tab.tsv",output_file),sep = "\t",row.names = F,col.names = T,quote = F,na = "")
      
      caller_types <- unique(stat_tab$caller)
      venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
      is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_all.tiff",output_file)
                            ,imagetype = "tiff"
                            ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                            ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                            ,cex = 1.5
                            ,cat.fontface=2
                            ,category.names=caller_types,main = "Germline varints all")
      
      
      stat_tab <- copy(var_tab)
      venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
      is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_pass.tiff",output_file)
                            ,imagetype = "tiff"
                            ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                            ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                            ,cex = 1.5
                            ,cat.fontface=2
                            ,category.names=caller_types,main = "Germline varints PASS")
    }

    var_tab[,call_info := paste0(round(variant_freq,2),",",coverage_depth)]
    
    
    setorder(var_tab,-variant_freq)
    
  } else {
    # READING VCF INTO A UNIFIED TABLE UNIFICATION OF READ COUNTS FROM DIFFERENT CALLERS
    var_tab <- data.table(type = vcfR::extract.info(vcf,"TYPE")
                          ,chrom = vcf@fix[,"CHROM"]
                          ,position = as.integer(vcf@fix[,"POS"])
                          ,reference = vcf@fix[,"REF"]
                          ,alternative = vcf@fix[,"ALT"]
                          ,filter = vcf@fix[,"FILTER"]
                          ,tumor_variant_freq = vcfR::extract.gt(vcf,element = "AD",as.numeric = T)[,2]
                          ,tumor_variant_freq_text = vcfR::extract.gt(vcf,element = "AD")[,2]
                          ,tumor_A = vcfR::extract.gt(vcf,element = "AU",as.numeric = T)[,2]
                          ,tumor_C = vcfR::extract.gt(vcf,element = "CU",as.numeric = T)[,2]
                          ,tumor_G = vcfR::extract.gt(vcf,element = "GU",as.numeric = T)[,2]
                          ,tumor_T = vcfR::extract.gt(vcf,element = "TU",as.numeric = T)[,2]
                          ,tumor_INDEL = vcfR::extract.gt(vcf,element = "TIR",as.numeric = T)[,2]
                          ,tumor_depth = vcfR::extract.gt(vcf,element = "DP",as.numeric = T)[,2]
                          ,normal_variant_freq = vcfR::extract.gt(vcf,element = "AD",as.numeric = T)[,1]
                          ,normal_variant_freq_text = vcfR::extract.gt(vcf,element = "AD")[,1]
                          ,normal_A = vcfR::extract.gt(vcf,element = "AU",as.numeric = T)[,1]
                          ,normal_C = vcfR::extract.gt(vcf,element = "CU",as.numeric = T)[,1]
                          ,normal_G = vcfR::extract.gt(vcf,element = "GU",as.numeric = T)[,1]
                          ,normal_T = vcfR::extract.gt(vcf,element = "TU",as.numeric = T)[,1]
                          ,normal_INDEL = vcfR::extract.gt(vcf,element = "TIR",as.numeric = T)[,1]
                          ,normal_depth = vcfR::extract.gt(vcf,element = "DP",as.numeric = T)[,1]
                          ,status = vcfR::extract.info(vcf,"STATUS"))

    var_tab[,index := seq_along(var_tab$chrom)]
    var_tab[!is.na(tumor_A) & alternative == "A",tumor_variant_freq := tumor_A]
    var_tab[!is.na(tumor_C) & alternative == "C",tumor_variant_freq := tumor_C]
    var_tab[!is.na(tumor_G) & alternative == "G",tumor_variant_freq := tumor_G]
    var_tab[!is.na(tumor_T) & alternative == "T",tumor_variant_freq := tumor_T]
    var_tab[!is.na(tumor_INDEL),tumor_variant_freq := tumor_INDEL]
    var_tab[!is.na(normal_A) & alternative == "A",normal_variant_freq := normal_A]
    var_tab[!is.na(normal_C) & alternative == "C",normal_variant_freq := normal_C]
    var_tab[!is.na(normal_G) & alternative == "G",normal_variant_freq := normal_G]
    var_tab[!is.na(normal_T) & alternative == "T",normal_variant_freq := normal_T]
    var_tab[!is.na(normal_INDEL),normal_variant_freq := normal_INDEL]

    var_tab[,c("tumor_A", "tumor_C","tumor_G", "tumor_T","tumor_INDEL","normal_A", "normal_C","normal_G", "normal_T","normal_INDEL") := NULL]

    var_tab[!is.na(vcfR::extract.gt(vcf,element = "F1R2")[,1]),normal_variant_freq := as.numeric(gsub(".*,","",normal_variant_freq_text))]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "F1R2")[,1]),tumor_variant_freq := as.numeric(gsub(".*,","",tumor_variant_freq_text))]
    var_tab[,tumor_variant_freq := round(tumor_variant_freq / tumor_depth,5)]
    var_tab[,normal_variant_freq := round(normal_variant_freq / normal_depth,5)]

    # ADD INFO ABOUT MULTIPLE CALLERS and GET UNIQUE VARIANTS
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "QSTD")[,1]),caller := "vardict"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "SUBDP")[,1]),caller := "strelka"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "SUBDP50")[,1]),caller := "strelka"]
    var_tab[!is.na(vcfR::extract.gt(vcf,element = "F1R2")[,1]),caller := "mutect2"]



    # GET CALLING STATISTICS
    stat_tab <- copy(var_tab)

    #FILTER CALLS
    var_tab <- var_tab[caller != "vardict" | grepl("Somatic",status),]
    var_tab <- var_tab[filter == "PASS"]
    var_tab <- var_tab[tumor_depth > 0 & tumor_variant_freq > 0 & normal_depth > 0]
    var_tab[tumor_variant_freq > 1,tumor_variant_freq := 1]

    if(nrow(var_tab) > 0){

      # GET CALLING STATISTICS
      stat_tab[,is_pass := filter == "PASS"]
      stat_tab[!(is.na(status) | status == "LikelySomatic" | status == "StrongSomatic"),status := "NotSomatic"]
      stat_tab[,status := factor(status,levels = c("NotSomatic","LikelySomatic","StrongSomatic"))]
      stat_tab[,by_caller_sum := .N,by = c("caller")]

      stat_tab_out1 <- stat_tab[,list(variants = .N),by = c("caller")]
      stat_tab_out1[,c("is_pass","status","percent_of_caller") := NA]
      stat_tab_out2 <- stat_tab[,list(variants = .N,percent_of_caller = round(.N / by_caller_sum[1] * 100,1)),by = c("caller","is_pass","status")]
      setcolorder(stat_tab_out1,names(stat_tab_out2))
      setorder(stat_tab_out1)
      setorder(stat_tab_out2)

      write.table(rbind(stat_tab_out1,stat_tab_out2),file = gsub(".tsv$",".stats_tab.tsv",output_file),sep = "\t",row.names = F,col.names = T,quote = F,na = "")

      caller_types <- unique(stat_tab$caller)
      venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
      is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_all.tiff",output_file)
                   ,imagetype = "tiff"
                   ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                   ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                   ,cex = 1.5
                   ,cat.fontface=2
                   ,category.names=caller_types,main = "Somatic variants all")

      stat_tab <- copy(var_tab)

      venn_list <- lapply(caller_types,function(x) stat_tab[caller == x,paste(chrom,position,sep = "_")])
      is_ok <- venn.diagram(venn_list,gsub(".tsv$",".stats_ven_pass.tiff",output_file)
                            ,imagetype = "tiff"
                            ,fill=c("red", "green","blue","yellow")[seq_along(caller_types)]
                            ,alpha=c(0.5,0.5,0.5,0.5)[seq_along(caller_types)]
                            ,cex = 1.5
                            ,cat.fontface=2
                            ,category.names=caller_types,main = "Somatic variants PASS")

      var_tab[,status := gsub("(.).*","\\1",status)]
      var_tab[caller == "vardict" ,caller := paste0(caller,"_",status)]
    }

    var_tab[,call_info := paste0("T:",round(tumor_variant_freq,2),",",tumor_depth,",N:",round(normal_variant_freq,2),",",normal_depth)]

    setorder(var_tab,-tumor_variant_freq)

  }
  
  var_tab <- var_tab[,c("caller_count", "callers","all_callers_info") := list(.N,paste(caller,collapse = ","),paste(call_info,collapse = ";")),by = c("chrom","position","reference","alternative")]
  var_tab <- unique(var_tab,by = c("chrom","position","reference","alternative"))
  
  if(nrow(var_tab) > 0){
    var_tab[,new_ref := reference]
    var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_ref := "-"]
    var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_ref := stringi::stri_sub(new_ref,from = nchar(alternative) + 1)]
    
    var_tab[,new_alt := alternative]
    var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_alt := "-"]
    var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_alt := stringi::stri_sub(new_alt,from = nchar(reference) + 1)]
    
    var_tab[,new_pos := position]
    var_tab[sapply(seq_along(reference),function(x) grepl(reference[x],alternative[x])),new_pos := new_pos + nchar(reference)]
    var_tab[sapply(seq_along(reference),function(x) grepl(alternative[x],reference[x])),new_pos := new_pos + nchar(alternative)]
    
    var_tab[,var_name := paste0(chrom,"_",new_pos,"_",new_ref,"/",new_alt)]
  } else {
    var_tab[,var_name := chrom]
  }

  
  vcf_out <- vcf[var_tab$index]
  vcfR::write.vcf(vcf_out,gsub(".tsv$",".vcf",output_file))
  
  if(variant_type == "germline"){
    var_tab <- var_tab[,list(var_name,genotype,variant_freq,coverage_depth,caller_count,callers,strand_bias)]
  } else {
    var_tab <- var_tab[,list(var_name,tumor_variant_freq,tumor_depth,normal_variant_freq,normal_depth,status,caller_count,callers,all_callers_info)]
  }
  
  write.table(var_tab,file = output_file,sep = "\t",row.names = F,col.names = T,quote = F,na = "")


}

# develop and test
# args <- character(3)
# args[1] <- "/mnt/ssd/ssd_1/snakemake/stage265_solid_tumors_children.reanalyze_30_36_WES/somatic_variant_calling/merged/AB1561.somatic.not_filtered.vcf"
# args[2] <- "/mnt/ssd/ssd_1/snakemake/stage265_solid_tumors_children.reanalyze_30_36_WES/somatic_variant_calling/merged/AB1561.somatic.tsv"
# args[3] <- "somatic"

# develop and test
# args <- character(3)
# args[1] <- "/mnt/ssd/ssd_1/snakemake/MedGen/sequencing_results/projects/Cosimo/pub_WES_Crescenzo/somatic_variant_calling/merged/GPS27_pos.somatic.not_filtered.vcf"
# args[2] <- "/mnt/ssd/ssd_1/snakemake/MedGen/sequencing_results/projects/Cosimo/pub_WES_Crescenzo/somatic_variant_calling/merged/GPS27_pos.somatic.tsv"
# args[3] <- "somatic"

# develop and test
# args <- character(3)
# args[1] <- "/mnt/ssd/ssd_1/snakemake/stage269_BRONCO.33/germline_variant_calling/merged/BR-0739.germline.not_filtered.vcf"
# args[2] <- "/mnt/ssd/ssd_1/snakemake/stage269_BRONCO.33/germline_variant_calling/merged/BR-0739.germline.tsv"
# args[3] <- "germline"

#run as Rscript

args <- commandArgs(trailingOnly = T)
run_all(args)
