#############################################################
# wrapper for rule: custom_annotation
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: custom_annotation \n##\n")
f.close()

command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/custom_annotation.R "+\
            snakemake.input.annotated + " " +\
            snakemake.output.custom_annotated + " " +\
            snakemake.params.reference_name + " " +\
            snakemake.params.anno_gtf + " " +\
            snakemake.params.resources_dir + " " +\
            snakemake.params.custom_DB_folder + " " +\
            snakemake.input.format_file +\
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
