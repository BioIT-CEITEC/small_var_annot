suppressMessages(library(data.table))
suppressMessages(library(parallel))

Sys.setenv("R_ZIPCMD" = "zip")


#DEFINE EXTRACTED DBs from Existing_variation vector
extracted_DBs <<- c("snpDB","COSMIC","HGMD","NHLBI_ESP")
extracted_DBs_IDs <<- c("rs","COS[VM]","^[CHB][IMSGXDRP][0-9]+$","ESP")

aa_names_tab <<- data.table(three = c("Ala" ,"Arg", "Asn" ,"Asp", "Cys", "Glu", "Gln", "Gly" ,"His" ,"Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser" ,"Thr", "Trp" ,"Tyr" ,"Val",".fs.*","Ter"),
                            one = c("A", "R","N" ,"D" ,"C","E", "Q", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V","fs","X"))

gene_trans_spec <<- unlist(list(
  BRCA1 = "ENST00000357654",
  BRCA2 = "ENST00000544455",
  TP53  = "NM_000546"
))


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


extract_existing_variations_to_separate_cols <- function(Existing_variation_vector){
  
  split_variant_list <- strsplit(Existing_variation_vector,",")
  
  res_mat <- sapply(split_variant_list,function(x) {
    if(x[1] != "-"){
      return(sapply(extracted_DBs_IDs,function(y) paste(x[grep(y,x)],collapse = ",")))
    } else {
      return(character(length(extracted_DBs_IDs)))
    }
  })
  res_mat <- t(res_mat)
  res_tab <- as.data.table(res_mat)
  names(res_tab) <- extracted_DBs
  return(res_tab)
}


add_custom_DB <- function(tab,custom_DBs,custom_DB_folder){
  location_tab <- as.data.table(tstrsplit(tab$var_name,"_",names = c("chrom","pos","t"),type.convert = T)[1:2])
  location_tab[,type := "b_var"][,index := seq_along(type)]
  setorder(location_tab)
  
  
  for(DB_text in custom_DBs){
    DB_name_split <- strsplit(DB_text,":")[[1]]
    DB_name <- DB_name_split[1]
    DB_filename <- list.files(custom_DB_folder,pattern = DB_name,full.names = T)
    if(length(DB_filename) == 1){
      orig_DB <- fread(DB_filename,sep = "\t",quote="")
      DB <- orig_DB[, c(T,lapply(orig_DB[,-1,with = FALSE], is.numeric) == TRUE), with = FALSE]
      orig_DB <- unique(orig_DB,by = names(DB)[1:3])
      setorder(orig_DB)
      DB <- unique(DB,by = names(DB)[1:3])
      if(length(names(DB)) > 2 && all(DB[[3]] - DB[[2]] >= 0)){
        DB_start <- DB[,c(1,2),with = F][,type := "a_start"][,index := seq_along(type)]
        setnames(DB_start,c(1:2),c("chrom","pos"))
        DB_end <- DB[,c(1,3),with = F][,type := "end"][,index := 0]
        setnames(DB_end,c(1:2),c("chrom","pos"))
        DB <- rbind(DB_start,DB_end)
        DB <- rbind(DB,location_tab)
        setorder(DB)
        indeces <- mapply(seq, which(DB$type == "a_start") + 1, which(DB$type == "end"))
        indeces <- lapply(indeces,function(x) x[-1])
        indeces_var <- unlist(indeces)
        indeces_start <- rep(seq_along(indeces),sapply(indeces,length))
        select_DB <- DB[indeces_var - 1][,start_idx := indeces_start]
        select_DB <- unique(select_DB[type == "b_var"],by = c("chrom","pos"))
        if(length(DB_name_split) > 1){
          cols <- strsplit(DB_name_split[2],"+",fixed = T)[[1]]
          cols <- cols[cols %in% names(orig_DB)]
          cols_with_DB_name <- paste(DB_name,cols,sep = "__")
          tab[,(cols_with_DB_name) := "."]
          tab[select_DB$index,(cols_with_DB_name) := orig_DB[select_DB$start_idx,cols,with = F]]
          for(col in cols_with_DB_name){
            tab[tab[[col]] == "",(col) := "."]
          }
          
        } else {
          tab[,(DB_name) := "."]
          tab[select_DB$index,(DB_name) := "yes"]
        }
      }
    }
  }
}


load_and_process_annot_tab <- function(annot_file,ref_name,col_config = NULL,resources_dir,custom_DB_folder,gtf_file){
  
  #READ ANNOTATION TABLE
  annot_tabFULL <- fread(annot_file,sep = "\t",header = T,skip = "#Uploaded_variation",verbose = F,showProgress = T)
  setnames(annot_tabFULL,"#Uploaded_variation","var_name")
  
  
  # SPLIT into chunks
  cores <- 100
  chunks <- parallel::splitIndices(nrow(annot_tabFULL), ncl = cores)
  ##############################################################################
  annot_tab_NEW <- data.table()
  
  for(i in 1:cores){
    # get partial chunk
    annot_tab <- annot_tabFULL[chunks[[i]],]
    
    annot_tab_extra_parse <- annot_tab$Extra
    annot_tab_extra_parse <- strsplit(annot_tab_extra_parse,";")
    annot_tab_extra_names_order <- order(unlist(lapply(annot_tab_extra_parse,seq_along)))
    annot_extra_names <- unlist(lapply(annot_tab_extra_parse,function(x) gsub("(.*)=.*","\\1",x)))
    annot_extra_value <- unlist(lapply(annot_tab_extra_parse,function(x) gsub(".*=(.*)","\\1",x)))
    annot_extra_index <- sapply(annot_tab_extra_parse,length)
    annot_extra_index <- rep(seq_along(annot_extra_index),annot_extra_index)
    annot_tab_extra <- data.table(index = annot_extra_index,names = annot_extra_names,value = annot_extra_value)
    annot_tab_extra <- dcast.data.table(annot_tab_extra,formula = index ~ names,fill = NA,value.var = "value")
    
    annot_tab <- cbind(annot_tab,annot_tab_extra)
    remove(annot_tab_extra)
    remove(annot_tab_extra_parse)
    remove(annot_extra_names)
    remove(annot_extra_value)
    remove(annot_extra_index)
    remove(annot_tab_extra_names_order)
    
    annot_tab[,Extra := NULL]
    annot_tab[,index := NULL]
    
    annot_tab <- unique(annot_tab,by = c("var_name","Location","Allele","Gene","Feature","STRAND"))
    annot_tab <- annot_tab[,var_name := gsub("(.*/[^/]*)/.*","\\1",var_name)]
    
    if(any(names(annot_tab) == "HGVSc")){
      annot_tab[,HGVSc := gsub(".*\\:c","c",HGVSc)]
    }
    if(any(names(annot_tab) == "HGVSp")){
      annot_tab[,HGVSp := gsub(".*\\:p","p",HGVSp)]
    }
    
    if(any(names(annot_tab) == "HGVSp")){
      annot_tab[grep("synonymous_variant",Consequence),HGVSp := gsub("p\\.([^0-9]+)([0-9]+)\\%3D","p\\.=(p\\.\\1\\2\\1)",HGVSp)]
      for(idx in 1:nrow(aa_names_tab)){
        annot_tab[,HGVSp := gsub(aa_names_tab$three[idx],aa_names_tab$one[idx],HGVSp)]
      }
    }
    
    # if(any(names(annot_tab) == "EUR_MAF")){
    #   annot_tab[,EUR_MAF := gsub("[ACGTN-]+:([01][^,]*).*","\\1",EUR_MAF)]
    # }
    
    if(any(names(annot_tab) == "Consequence")){
      annot_tab[,gene_region := "exon"]
      annot_tab[grep("5_prime_UTR|splice|intron|downstream|upstream|3_prime_UTR|intergenic|mature_miRNA|non_coding|regulatory_region|TF_binding_site",Consequence)
                ,gene_region := gsub(".*(5_prime_UTR|splice|intron|downstream|upstream|3_prime_UTR|intergenic|mature_miRNA|non_coding|regulatory_region|TF_binding_site).*","\\1",Consequence)]
    }
    
    if(any(names(annot_tab) == "Existing_variation")){
      add_tab <- cbind(Existing_variation = unique(annot_tab$Existing_variation)
                       ,extract_existing_variations_to_separate_cols(unique(annot_tab$Existing_variation)))
      annot_tab <- merge(annot_tab,add_tab,by = "Existing_variation")
    }
    
    if(any(names(annot_tab) == "CLIN_SIG")){
      clinvar_tab <- fread(paste0(custom_DB_folder,"/clinvar_annot.tsv"))
      annot_tab <- merge(annot_tab,clinvar_tab,by = "var_name",all.x = T)
      annot_tab[,CLIN_SIG := NULL]
    }
    
    if(!any(names(annot_tab) == "SYMBOL")){
      if(file.exists(gtf_file)){
        gtf <- fread(gtf_file)
        SYMBOL_annot <- unique(grep("gene_id .*; gene_name.*gene_biotype",gtf$V9,value = T))
        SYMBOL_annot <- unique(gsub(".*gene_id \"(.*)\"; gene_name \"(.*)\";.*gene_biotype \"([^ ;\"]*).*","\\1;\\2;\\3",SYMBOL_annot))
        SYMBOL_annot_tab <- tstrsplit(SYMBOL_annot,";")
        SYMBOL_annot_tab <- data.table(Gene = SYMBOL_annot_tab[[1]],SYMBOL = SYMBOL_annot_tab[[2]],BIOTYPE = SYMBOL_annot_tab[[3]])
        annot_tab <- merge(annot_tab,SYMBOL_annot_tab,by = "Gene",all.x = T)
        annot_tab[is.na(SYMBOL),SYMBOL := Gene]
        annot_tab[is.na(BIOTYPE),BIOTYPE := "not_available"]
      } else {
        annot_tab[,SYMBOL := Gene]
        annot_tab[,BIOTYPE := "not_available"]
      }
    }
    
    annot_tab[,c("chrom","pos","substitution") := as.data.table(tstrsplit(var_name,"_"))]
    annot_tab[,c("reference","alternative") := as.data.table(tstrsplit(substitution,"/"))]
    annot_tab[,pos := as.integer(pos)]
    
    
    for(col in names(annot_tab)) {
      if(is.character(annot_tab[[col]])){
        set(annot_tab, i=which(annot_tab[[col]]=="" | is.na(annot_tab[[col]])), j=col, value=".")
      }
    }
    
    if(any(names(annot_tab) == "EXON")){
      annot_tab[EXON != ".",exon_intron_number := paste0("exon",gsub("\\/.*","",EXON))]
      if(any(names(annot_tab) == "INTRON")){
        annot_tab[INTRON != ".",exon_intron_number := paste0("intron",gsub("\\/.*","",INTRON))]
      }
      annot_tab[is.na(exon_intron_number), exon_intron_number := "."]
      annot_tab[,full_annot_name := paste(SYMBOL,Feature,exon_intron_number,HGVSc,HGVSp,sep = ":")]
    } else {
      annot_tab[,full_annot_name := ""]
    }
    
    consequence_tab <- fread(paste0(resources_dir,"/consequences_tab.tsv"),header = F)
    consequence_tab[,consequence_index := seq_along(V1)]
    consequence_map <- data.table(Consequence = unique(annot_tab$Consequence),V1 = sapply(strsplit(unique(annot_tab$Consequence),","),head,1))
    consequence_map <- merge(consequence_map,consequence_tab[,list(V1,consequence_index)],by = "V1")
    annot_tab <- merge(annot_tab,consequence_map[,list(Consequence,consequence_index)],by = "Consequence")
    annot_tab[,specific_transcript := Feature %in% gene_trans_spec]
    
    annot_tab[,is_protein_coding := BIOTYPE == "protein_coding"]
    
    non_protein_to_remove <- function(is_protein_coding){
      if(any(is_protein_coding)){
        return(!is_protein_coding)
      } 
      else {
        return(FALSE)
      }
    }
    
    annot_tab[,non_protein_coding_to_remove := non_protein_to_remove(is_protein_coding),by = c("var_name")]
    annot_tab <- annot_tab[non_protein_coding_to_remove == F]
    
    columns_to_order <- c("specific_transcript","is_protein_coding","SOURCE","consequence_index","CANONICAL")
    order_direction <- c(-1,-1,-1,1,-1)
    order_direction <- order_direction[columns_to_order %in% names(annot_tab)]
    columns_to_order <- columns_to_order[columns_to_order %in% names(annot_tab)]
    setorderv(annot_tab,cols = columns_to_order,order = order_direction,na.last = T)
    
    annot_tab[,all_full_annot_name := paste(unique(full_annot_name),collapse = ","),by = c("var_name","SYMBOL")]
    annot_tab <- unique(annot_tab,by = c("var_name","SYMBOL"))
    
    setorder(annot_tab,chrom,pos)
    
    if(!is.null(col_config)){
      if(any(grepl("custom_DB::",col_config$orig_name))){
        custom_DBs <- col_config[grepl("custom_DB::",orig_name)]
        custom_DBs[,orig_name := sub("custom_DB::","",orig_name)]
        custom_DBs[,DB := gsub("__.*","",orig_name)] 
        custom_DBs[,cols := ""]
        custom_DBs[grepl("__",orig_name),cols := paste(gsub(".*__","",orig_name),collapse = "+"),by = DB]
        custom_DBs[,input := paste(DB,cols,sep = ":")]
        custom_DBs[,input := gsub("\\:$","",input)]
        add_custom_DB(annot_tab,unique(custom_DBs$input),custom_DB_folder)
      } 
      
      if(any(col_config$orig_name == "exon_dist") && file.exists(gtf_file)){
        gtf <- fread(gtf_file,header = F)
        gtf <- gtf[V3 == "exon",.(chrom = V1,start = V4,end = V5,gene_name = gsub(".*gene_name .([^;]*).;.*","\\1",V9),exon = T)]
        gtf <- rbind(gtf[,.(chrom = chrom,pos = end,gene_name,exon = T)],gtf[,.(chrom = chrom,pos = start,gene_name,exon = T)])
        gtf <- unique(gtf)
        
        pos_tab <- annot_tab[gene_region != "exon",.(chrom,pos,gene_name = SYMBOL,exon = F)]
        pos_tab <- unique(pos_tab)
        
        exon_pos_tab <- rbind(gtf,pos_tab)
        exon_pos_tab[,pos := as.integer(pos)]
        exon_pos_tab[,any_var := sum(!exon),by = c("chrom","gene_name")]
        exon_pos_tab <- exon_pos_tab[any_var > 0]
        setkey(exon_pos_tab)
        exon_pos_tab[,any_var := NULL]
        exon_pos_tab[,diff := abs(c(diff(exon),0)) + abs(c(0,diff(exon))), by = c("chrom","gene_name")]
        exon_pos_tab <- exon_pos_tab[exon == F | diff != 0,]
        exon_pos_tab[,exon_dist := sapply(pos, function(x) as.numeric(min(abs(x - pos[exon == T])))),by = c("chrom","gene_name")]
        setnames(exon_pos_tab,"gene_name","SYMBOL")
        
        annot_tab <- merge(annot_tab,exon_pos_tab[exon == F,.(chrom,pos,SYMBOL,exon_dist)],by = c("chrom","pos","SYMBOL"),all.x = T)
        annot_tab[is.na(exon_dist),exon_dist := 0]
        annot_tab[exon_dist == Inf,exon_dist := NA]
      }
    } 
    
    #var_gen_coord == Chr5(GRCh37):g.112176756T>A
    annot_tab[,var_gen_coord := paste0(ref_name,"____",var_name)]
    annot_tab[,var_gen_coord := gsub("^(.*)____([0-9XY]+)_([0-9]+)_([ACGT-]+)\\/([ACGT-]+)$","Chr\\2\\(\\1\\)\\:g\\.\\3\\4\\>\\5",var_gen_coord)]
    
    
    #correct renaming
    # final_tab[,full_trans_name := gsub("\\:c.*","",HGVSc)]
    
    
    annot_tab_NEW <- rbindlist(list(annot_tab_NEW, annot_tab), fill=TRUE)
    
    message(sprintf("%i out of %s iterations DONE", i, cores ))
    
    rm(annot_tab)
    suppressMessages(gc())
    
  }
  ##############################################################################
  # FOR LOOP END
  return(annot_tab_NEW)
  
}

run_all <- function(args){
  #read all params as variables
  annot_file <- args[1]
  output_file <- args[2]
  ref_name <- args[3]
  gtf_file <- args[4]
  resources_dir <- args[5]
  custom_DB_folder <- args[6]
  format_file <- args[7]
  
  
#   ref_name <- gsub("\\-.*","",ref_name)
  
  col_config <- fread(format_file,skip = "orig_name")
  # add gnomad WGS column to format
  x <- which(col_config$orig_name=="gnomADe_NFE_AF")
  if(!length(x)){x <- which(col_config$orig_name=="gnomAD_NFE_AF")}
  col_config <- as.data.table(rbind(col_config[1:x,],
                                    data.frame(orig_name="gnomADg_NFE_AF", new_name="gnomAD_WGS_NFE"),
                                    col_config[(x+1):nrow(col_config),]
                                    ))
                              
                              
  #load and process annotated vars
  annot_tab <- load_and_process_annot_tab(annot_file,ref_name,col_config,resources_dir,custom_DB_folder,gtf_file)
  
  fwrite(annot_tab,file = output_file,sep = "\t")
}



# develop and test WES 86
# setwd("/mnt/ssd/ssd_1/workspace/ROBIN2/MOII_signatures_GERMLINE/")
# args <- c("annotate/all_variants.annotated.tsv",
#           "annotate/all_variants.annotated.processed.tsv",
#           "GRCh38-p10","/mnt/references/homo_sapiens/GRCh38-p10/annot/GRCh38-p10.gtf",
#           "/home/402182/PIPELINY/small_var_annot-v2.0/resources",
#           "/mnt/references/homo_sapiens/GRCh38-p10/annot/custom_new2",
#           "/mnt/references/general/germline_small_var_call_format_files/slaby_children_germ.txt")


#run as Rscript
# 
args <- commandArgs(trailingOnly = T)
run_all(args)



