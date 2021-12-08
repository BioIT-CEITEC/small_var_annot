######################################
# wrapper for rule: normalize_variants
######################################
import re
import os
import subprocess
from snakemake.shell import shell

log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: normalize_variants - "+ snakemake.wildcards.variant_caller +" \n##\n")
f.close()

shell.executable("/bin/bash")

# os.path.abspath(os.path.dirname(__file__))+"/../../adapters_merge.fa"

if os.stat(snakemake.input.vcf).st_size > 0:

    shell.executable("/bin/bash")

    if str(snakemake.wildcards.variant_caller) == "vardict":
        with open(snakemake.input.vcf, "r") as f:
            lines = f.readlines()
        with open(snakemake.input.vcf, "w") as f:
            for line in lines:
                if re.search("^\\t\\t\.",line) == None:
                    f.write(line)

    f = open(log_filename, 'at')
    command = "bgzip -c "+ snakemake.input.vcf +"  > "+ snakemake.input.vcf +".gz"
    f.write("## COMMAND: "+command+"\n")
    shell(command)

    f = open(log_filename, 'at')
    command = "tabix -p vcf " + snakemake.input.vcf + ".gz"
    f.write("## COMMAND: "+command+"\n")
    shell(command)


    version = str(subprocess.Popen("bcftools 2>&1 | grep \"[Vv]ersion\" | cut -f 2 -d \" \"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f = open(log_filename, 'at')
    f.write("## VERSION: bcftools "+version+"\n")
    f.close()


    if str(snakemake.wildcards.variant_caller) == "vardict":
        f = open(log_filename, 'at')
        command = "bcftools norm --multiallelics '-both' --check-ref x -f "+ snakemake.input.ref + " -o /dev/stdout " + snakemake.input.vcf + ".gz | "+ \
                  "bcftools annotate -x 'INFO/END' /dev/stdin > " + snakemake.output.vcf
        f.write("## COMMAND: "+command+"\n")
        shell(command)
    else:
        f = open(log_filename, 'at')
        command = "bcftools norm --multiallelics '-both' --check-ref x -f "+ snakemake.input.ref + " -o "+ snakemake.output.vcf + " " + snakemake.input.vcf + ".gz" + \
                  " >> " + log_filename + " 2>&1"
        f.write("## COMMAND: "+command+"\n")
        shell(command)


    # adding all contigs to vcf from .dict files because of the stupid GATK
    f = open(log_filename, 'at')
    command = "Rscript "+os.path.abspath(os.path.dirname(__file__))+"/add_contigs_to_vcf.R " + snakemake.output.vcf + " " + snakemake.input.dict +\
              " >> " + log_filename + " 2>&1"
    f.write("## COMMAND: "+command+"\n")
    shell(command)


    # normalizing sample names in genotype fields to donor.tag
    command = "sed -i 's/\tFORMAT\t[^\t]*/\tFORMAT\t" + snakemake.wildcards.sample_name + "/' " + snakemake.output.vcf

    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    shell(command)

    # sort variants with reference dict
    f = open(log_filename, 'at')
    command = "gatk SortVcf -I " + snakemake.output.vcf + " -O " + snakemake.output.vcf + " -SD " + snakemake.input.dict + " >> " + log_filename + " 2>&1"
    f.write("## COMMAND: "+command+"\n")
    shell(command)

else:
    command = "touch " + snakemake.output.vcf + " >> " + log_filename + " 2>&1 "
    f = open(log_filename, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
