######################################
# wrapper for rule: variant_annotation
######################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: variant_annotation \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("vep 2>&1 | grep \"ensembl-vep\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## VERSION: ensembl-vep  :  "+version+"\n")
f.close()

if sum(1 for line in open(snakemake.input.tsv_for_vep)) > 1:
    assembly = re.sub("\-.*","",snakemake.params.ref_name)

    if snakemake.threads > 1:
        fork_text = "--fork " + str(snakemake.threads)
    else:
        fork_text = ""

    cache_version = 0

    organism = snakemake.params.organism_name

    if organism == "homsap":
        organism = "homo_sapiens"

    if os.path.isdir(snakemake.params.vep_dir):

        vep_merged_dir = snakemake.params.vep_dir + "/" + organism + "_merged"
        if os.path.isdir(vep_merged_dir) and not snakemake.params.not_use_merged:
            merged = "--merged"
            for dir in os.listdir(vep_merged_dir):
                 if re.search("^[0-9]+_",dir):
                     num = int(re.sub("_.*","",dir))
                     if num > cache_version:
                         cache_version = num
        else:
            merged = ""
            vep_dir = snakemake.params.vep_dir + "/" + organism
            for dir in os.listdir(vep_dir):
                 if re.search("^[0-9]+_",dir):
                     num = int(re.sub("_.*","",dir))
                     if num > cache_version:
                         cache_version = num


        if organism == "homo_sapiens":
            command = "vep --dir " + snakemake.params.vep_dir +\
                           " --everything "  + merged +\
                           " --fasta "+ snakemake.params.ref +\
                           " --offline --assembly " + assembly +\
                           " --cache_version " + str(cache_version) +\
                           " --input_file " + snakemake.input.tsv_for_vep +\
                           " --output_file " + snakemake.output.annotated +\
                           " --force_overwrite " + \
                           " --dir_plugins " + snakemake.params.dir_plugins + \
                           " --plugin CADD," + snakemake.params.CADD_DB_SNVs + ","  + snakemake.params.CADD_DB_indels + " " +\
                           fork_text + " >> " + log_filename + " 2>&1"
        else:
            command = "vep --dir " + snakemake.params.vep_dir + \
                      " --everything " + merged + \
                      " --fasta " + snakemake.params.ref + \
                      " --species " + organism + \
                      " --offline --assembly " + assembly + \
                      " --cache_version " + str(cache_version) + \
                      " --input_file " + snakemake.input.tsv_for_vep + \
                      " --output_file " + snakemake.output.annotated + \
                      " --force_overwrite " + \
                      fork_text + " >> " + log_filename + " 2>&1"

        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)


    else:
        annot_dir = snakemake.params.vep_dir.replace("/vep","")
        gff_file = "no_gff_ERROR"
        for filename in os.listdir(annot_dir):
            if re.search("gff.*.gz.tbi", filename) is not None:
                gff_file = annot_dir + "/" + filename.replace(".tbi","")


        command = "vep --gff " + gff_file +\
                       " --fasta "+ snakemake.params.ref +\
                       " --input_file " + snakemake.input.tsv_for_vep +\
                       " --output_file " + snakemake.output.annotated +\
                       " --force_overwrite " +\
                       fork_text + " >> " + log_filename + " 2>&1"


        f = open(log_filename, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
else:
    command = "touch " + snakemake.output.annotated
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
