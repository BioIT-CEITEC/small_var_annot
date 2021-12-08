######################################
# wrapper for rule: haplotypecaller
######################################
import os
import subprocess
from snakemake.shell import shell


log_filename = str(snakemake.log)

f = open(log_filename, 'wt')
f.write("\n##\n## RULE: haplotypecaller \n##\n")
f.close()


shell.executable("/bin/bash")

# ZAKOMENTOVANO 23.11.2020 small variants analysis
# normal_cfg = snakemake.params.normal_sample
# normal_sample_name = str(normal_cfg.loc[normal_cfg.fullname == snakemake.wildcards.fullname,"sample"].min())
# if normal_sample_name == "nan":
#     normal_sample_name = normal_cfg["sample"].min()

version = str(subprocess.Popen("gatk HaplotypeCaller --version true 2>&1 | grep \"[Vv]ersion:\" | cut -f 2 -d \":\"", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'a+')
f.write("## VERSION: gatk "+version+"\n")
f.close()

if snakemake.params.lib_ROI == "wgs":
    intervals_call = ""
else:
    intervals_call = " -L " + snakemake.input.regions

command = "mkdir -p " + os.path.dirname(snakemake.params.bamout)
shell(command)


command = "gatk --java-options \"-Xmx2g\" HaplotypeCaller" + \
               " -R "+ snakemake.input.ref +\
               " -I "+ snakemake.input.bam +\
               intervals_call +\
               " --native-pair-hmm-threads "+ str(snakemake.threads) +\
               " -O " + snakemake.output.vcf +\
               " -bamout "+ snakemake.params.bamout +\
               " >> " + log_filename + " 2>&1 "


f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


