#############################################################
# wrapper for rule: process_and_format_annot_variants
#############################################################
import os
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: process_and_format_annot_variants \n##\n")
f.close()


if hasattr(snakemake.input, 'cohort_data') and len(snakemake.input.cohort_data) > 0:
    cohort_data_filename = snakemake.input.cohort_data
else:
    cohort_data_filename = "no_cohort_data"

if snakemake.params.create_cohort_data:
    create_cohort_data = "cohort_data/cohort_variants.tsv"
else:
    create_cohort_data = "dont_save_cohort_data"

if snakemake.params.isWGS == "wgs" or snakemake.params.isWGS == "rna" or os.path.getsize(snakemake.input.annotated) > 5 * 10 ** 9:
    command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_and_format_annot_variants_WGS.R "+\
            snakemake.input.annotated + " " +\
            snakemake.output.all_vars_tsv + " " +\
            os.path.dirname(snakemake.output.per_sample_var_tabs[0]) + " " +\
            snakemake.input.format_file + " " +\
            snakemake.params.min_variant_frequency + " " +\
            cohort_data_filename + " " +\
            create_cohort_data + " " +\
            snakemake.params.batch_name + " " +\
            snakemake.params.ref_dir + " " +\
            snakemake.params.organism + " " +\
            snakemake.params.mut_load_output_filename + " " +\
            " ".join(snakemake.input.var_tabs) +\
            " >> " + log_filename + " 2>&1"
else:
    command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/process_and_format_annot_variants.R "+\
            snakemake.input.annotated + " " +\
            snakemake.output.all_vars_tsv + " " +\
            os.path.dirname(snakemake.output.per_sample_var_tabs[0]) + " " +\
            snakemake.input.format_file + " " +\
            snakemake.params.min_variant_frequency + " " +\
            cohort_data_filename + " " +\
            create_cohort_data + " " +\
            snakemake.params.batch_name + " " +\
            snakemake.params.ref_dir + " " +\
            snakemake.params.organism + " " +\
            snakemake.params.mut_load_output_filename + " " +\
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
