suppressMessages(library(data.table))
Sys.setenv("R_ZIPCMD" = "zip")

#ZMENENO 7.12.2020
fread_vector_of_files <- function(file_list,regex,add_column = "sample"){
  list_of_tabs <- lapply(seq_along(file_list),function(x){
    res <- fread(file_list[x])
    res[,(add_column) := gsub(regex,"\\1",file_list[x])]
    if(nrow(res) > 0){
      return(res)
    } else {
      return(NULL)
    }
  })  
  if(any(sapply(list_of_tabs,is.null))){
    empty_sample_names <<- gsub(regex,"\\1",file_list[sapply(list_of_tabs,is.null)])
  }  
  return(rbindlist(list_of_tabs))
}


annotate_with_intervals <- function(var_tab,annot_tab,annotate_cols_names = tail(names(annot_tab),1)){
  location_tab <- var_tab[,list(chrom,start = pos,end = pos)]
  location_tab <- unique(location_tab)
  setnames(annot_tab,names(annot_tab)[1:3],c("chrom", "start", "end"))
  setkey(annot_tab, chrom, start, end)
  res <- foverlaps(location_tab, annot_tab, type="any")
  res <- res[,c("chrom","i.start",annotate_cols_names),with = F]
  setnames(res,names(res)[1:2],c("chrom", "pos"))
  var_tab <- merge(var_tab,res,by = c("chrom","pos"))
  return(var_tab)
}

run_all <- function(args){
  #read all params as variables
  annot_file <- args[1]
  output_file <- args[2]
  per_sample_results_dir <- args[3]
  format_file <- args[4]
  VF_threshold <- as.numeric(args[5]) / 100
  cohort_data_filename <- args[6]
  var_files <- args[7:length(args)]


  cohort_data <- fread(cohort_data_filename)

  #load format config file
  full_format_configs <- readLines(format_file)
  global_format_configs <- data.table(gsub("=.*","",full_format_configs[1:(which(full_format_configs == "") - 1)[1]]),
                                      gsub(".*?=(.*)","\\1",full_format_configs[1:(which(full_format_configs == "") - 1)[1]]))
  
  col_config <- fread(format_file,skip = "orig_name")
  col_config[,orig_name := sub(".*::","",orig_name)]
  
  #load processed annotated vars
  annot_tab <- fread(annot_file)
  
  #load and filter vars from all samples
  sample_filename_pattern <- paste0(".*\\/(.*).processed.tsv")
  all_var_tab <- fread_vector_of_files(var_files,regex = sample_filename_pattern)
  all_var_tab <- filter_variants(all_var_tab,VF_threshold = VF_threshold)
  
  final_unformated_tab <- merge(all_var_tab,annot_tab,by = "var_name",allow.cartesian=TRUE)
  
  #keep only cols in config which are in the final table
  col_config <- col_config[orig_name %in% names(final_unformated_tab)]

  final_formated_tab <- format_final_var_table(final_unformated_tab,global_format_configs,col_config)
  
  write_out_per_sample_vars(final_formated_tab,per_sample_results_dir,full_format_configs,output_file,col_config)
  
  #print cast var_table
  final_formated_tab[,is_in := 1]
  var_sample_presence <- data.table::dcast.data.table(final_formated_tab,formula = var_name ~ sample,fun.aggregate = identity,fill = 0,value.var = "is_in")
  final_formated_tab[,is_in := NULL]
  unique_var_tab <- unique(final_formated_tab,by = c("var_name","Gene_symbol"))
  unique_var_tab[,sample := NULL]
  cast_var_tab <- merge(var_sample_presence,unique_var_tab,by = "var_name")
  setcolorder(cast_var_tab,c("var_name","Gene_symbol"))
  
  fwrite(cast_var_tab,file = output_file,sep = "\t")
  openxlsx::write.xlsx(list(all_vars_cast=cast_var_tab),file = gsub(".tsv",".xlsx",output_file))
  
}

filter_variants <- function(all_var_tab,VF_threshold = 0,coverage_alarm = c(1,10,40),single_transcript = T){
  
  all_var_tab <- all_var_tab[variant_freq > VF_threshold,]
  all_var_tab[,alarm := ""]
  all_var_tab[coverage_depth < coverage_alarm[3],alarm := "Low coverage"]
  all_var_tab[coverage_depth < coverage_alarm[2],alarm := "Very low coverage"]
  all_var_tab[coverage_depth < coverage_alarm[1],alarm := "No coverage"]
  
  return(all_var_tab)
}


compute_and_write_mut_load  <- function(variant_tab,mut_load_output_file,global_format_configs){
  variant_tab <- unique(variant_tab,by = c("sample","var_name"))
  mut_load_config <- as.data.table(tstrsplit(global_format_configs[V1 == "mut_load"]$V2,split = "::"))
  
  mut_load_res_tab <- data.table(sample = unique(variant_tab$sample))
  for(index in seq_along(mut_load_config$V1)){
    
    filter_text <- trimws(mut_load_config[index,]$V3)
    filtered_var_table <- eval(parse(text = paste0("variant_tab[",filter_text,"]")))
    intervals <- fread(mut_load_config[index,]$V2)
    intervals[,is_in := "x"]
    
    filtered_var_table <- annotate_with_intervals(filtered_var_table,intervals,annotate_cols_names = "is_in")
    filtered_var_table <- filtered_var_table[!is.na(is_in)]
    tab <- filtered_var_table[,list(mutation_load = round(.N * 10^6 / sum(intervals$end - intervals$start + 1),2)),by = sample]
    mut_load_res_tab <- merge(mut_load_res_tab,tab,by = "sample",all.x = T)
    mut_load_res_tab[is.na(mutation_load),mutation_load := 0]
    setnames(mut_load_res_tab,"mutation_load",mut_load_config[index,]$V1)
    
  }
  openxlsx::write.xlsx(mut_load_res_tab,file = mut_load_output_file)
}



format_final_var_table  <- function(variant_tab,global_format_configs,col_config){
  
  #apply format specific rounding
  rounding <- as.numeric(global_format_configs[V1 == "rounding"]$V2)
  cols <- names(variant_tab)[sapply(variant_tab,is.numeric)]
  variant_tab[,(cols) := round(.SD,rounding), .SDcols=cols]
  
  #config specific sorting
  order_by_config <- strsplit(global_format_configs[V1 == "sort_by"]$V2,",")[[1]]
  order_direction <- sign(-as.numeric(grepl("^-",order_by_config)) + 0.5)
  order_by_config <- gsub("^-","",order_by_config)
  setorderv(variant_tab,cols = order_by_config[order_by_config %in% names(variant_tab)],order = order_direction[order_by_config %in% names(variant_tab)],na.last = T)
  
  #ZMENENO 7.12.2020
  orig_names_vec <- c("var_name",col_config$orig_name)
  new_names_vec <- c("var_name",col_config$new_name)  
  # orig_names_vec <- c("var_name","user_annotation","comment",col_config$orig_name,"DB_upload_info_projectsample_id")
  # new_names_vec <- c("var_name","user_annotation","comment",col_config$new_name,"DB_upload_info_projectsample_id")
  
  if(!any("sample" == col_config$orig_name)){
    orig_names_vec <- c("sample",orig_names_vec)
    new_names_vec <- c("sample",new_names_vec)
  }
  
  variant_tab <- variant_tab[,orig_names_vec,with = F]
  setnames(variant_tab,orig_names_vec,new_names_vec)
  
  if(any(names(variant_tab) == "1000g_EUR_AF")){
    variant_tab <- suppressWarnings(variant_tab[,`1000g_EUR_AF` := as.numeric(`1000g_EUR_AF`)])
    variant_tab <- variant_tab[is.na(`1000g_EUR_AF`),`1000g_EUR_AF` := 0]
  }
  
  if(any(names(variant_tab) == "gnomAD_NFE")){
    variant_tab <- suppressWarnings(variant_tab[,gnomAD_NFE := as.numeric(gnomAD_NFE)])
    variant_tab <- variant_tab[is.na(gnomAD_NFE),gnomAD_NFE := 0]
  }
  
  return(variant_tab)
}

write_out_per_sample_vars  <- function(variant_tab,per_sample_results_dir,full_format_configs,output_file,col_config){
  
  if(!dir.exists(per_sample_results_dir)){
    dir.create(per_sample_results_dir)
  }
  
  if(!dir.exists(paste0(per_sample_results_dir,"/tsv_formated/"))){
    dir.create(paste0(per_sample_results_dir,"/tsv_formated/"))
  }
  
  if(any(grepl("^filtered_res:",full_format_configs))){
    filtered_config <- full_format_configs[(grep("^filtered_res:",full_format_configs) + 1):length(full_format_configs)]
    filtered_config <- filtered_config[1:(which(filtered_config == "") - 1)[1]]
    filtered_config <- as.data.table(tstrsplit(filtered_config,":"))
    
  } else {
    filtered_config <- NULL
  }
  

  for(my_sample in unique(variant_tab$sample)){
    sample_tab <- variant_tab[sample == my_sample]
    if(!any("sample" == col_config$orig_name)){
      sample_tab[,sample := NULL]
    }

    
    fwrite(sample_tab,file = paste0(per_sample_results_dir,"/tsv_formated/",my_sample,".variants.tsv"),sep = "\t")
    
    if(!is.null(filtered_config)){
      
      res_list <- list(all = sample_tab)
      res_list_add <- lapply(filtered_config$V2,function(x){
        eval(parse(text = paste0("sample_tab[",x,"]")))
      })
      names(res_list_add) <-  filtered_config$V1
      
      openxlsx::write.xlsx(c(res_list,res_list_add),file = paste0(per_sample_results_dir,"/",my_sample,".variants.xlsx"))
      
    } else {
      openxlsx::write.xlsx(sample_tab,file = paste0(per_sample_results_dir,"/",my_sample,".variants.xlsx"))
    }
  }
  
  if(length(empty_sample_names) > 0){
    sample_tab <- variant_tab[0,]
    if(!any("sample" == col_config$orig_name)){
      sample_tab[,sample := NULL]
    }

    
    for(empty_sample_name in empty_sample_names){
      openxlsx::write.xlsx(sample_tab,file = paste0(per_sample_results_dir,"/",my_sample,".variants.xlsx"))
      fwrite(unformated_sample_tab,file = paste0(per_sample_results_dir,"/tsv_formated/",my_sample,".variants.tsv"),sep = "\t")
    }
  }

}


empty_sample_names <<- character()

# develop and test
# setwd("/mnt/ssd/ssd_1/snakemake/stage359_PC.seq_A/somatic_variant_calling")
# setwd("/mnt/ssd/ssd_1/sequia/211206__germline_small_var_call__1332")
# args <- c("annotate/all_variants.annotated.processed.tsv","final_variant_table.tsv","mutation_loads.xlsx","per_sample_final_var_tabs","~/BioRoots/germline_small_var_call/resources/formats/default.txt","25","merged/BR-1052.processed.tsv","merged/BR-1053.processed.tsv","merged/BR-1054.processed.tsv","merged/BR-1055.processed.tsv","merged/BR-1056.processed.tsv","merged/BR-1060.processed.tsv","merged/BR-1061.processed.tsv","merged/BR-1062.processed.tsv","merged/BR-1063.processed.tsv","merged/BR-1064.processed.tsv","merged/BR-1065.processed.tsv","merged/BR-1066.processed.tsv","merged/H-0381.processed.tsv","merged/H-0388.processed.tsv","merged/H-0389.processed.tsv","merged/H-0397.processed.tsv")

#run as Rscript
args <- commandArgs(trailingOnly = T)
run_all(args)



