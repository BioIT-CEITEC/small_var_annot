#############################################################
# wrapper for rule: process_and_format_annot_variants
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: process_and_format_annot_variants \n##\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_and_format_annot_variants.R "+\
            snakemake.input.annotated + " " +\
            snakemake.output.all_vars_tsv + " " +\
            os.path.dirname(snakemake.output.per_sample_var_tabs[0]) + " " +\
            snakemake.input.format_file + " " +\
            snakemake.params.min_variant_frequency + " " +\
            snakemake.input.cohort_data + " " +\
            " ".join(snakemake.input.var_tabs) +\
            " >> " + log_filename + " 2>&1"


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## args <- c(\"" + "\",\"".join(command.split(" ")[2:-3]) + "\")\n")
f.close()

shell(command)

# os.makedirs(os.path.join(os.path.dirname(snakemake.output.all_vars),"user_annotations"),exist_ok = True)
#
# command = "cp " + os.path.dirname(snakemake.input.var_tabs[0]) + "/*.xlsx " + os.path.join(os.path.dirname(snakemake.output.all_vars),"user_annotations")
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
