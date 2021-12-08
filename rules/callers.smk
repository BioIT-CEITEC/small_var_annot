
def bam_input(wildcards):
    if config["material"] != "RNA":
        tag = "bam"
    else:
        tag = "RNAsplit.bam"

    return expand("mapped/{input_bam}.{tag}",input_bam=wildcards.sample_name,tag=tag)[0]


rule haplotypecaller:
    input:  bam = bam_input,
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: vcf="variant_calls/{sample_name}/haplotypecaller/haplotypecaller.vcf"
    log: "logs/{sample_name}/callers/haplotypecaller.log"
    threads: 5
    resources:
        mem_mb=6000
    params: bamout="variant_calls/{sample_name}/haplotypecaller/realigned.bam",
            lib_ROI=config["lib_ROI"]
    conda: "../wrappers/haplotypecaller/env.yaml"
    script: "../wrappers/haplotypecaller/script.py"

rule vardict:
    input:  bam = bam_input,
            ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            refdict=expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
    output: vcf="variant_calls/{sample_name}/vardict/vardict.vcf"
    log: "logs/{sample_name}/callers/vardict.log"
    threads: 10
    resources:
        mem_mb=8000
    params:
        AF_threshold=config["min_variant_frequency"]
    conda: "../wrappers/vardict/env.yaml"
    script: "../wrappers/vardict/script.py"

rule strelka:
    input:  bam = bam_input,
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            regions_gz = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed.gz",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
            regions_tbi = expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed.gz.tbi",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
    output: vcf = "variant_calls/{sample_name}/strelka/strelka.vcf"
    log: "logs/{sample_name}/callers/strelka.log"
    threads: 10
    resources:
        mem_mb=6000
    params: dir="variant_calls/{sample_name}/strelka",
            lib_ROI=config["lib_ROI"],
            vcf= "variant_calls/{sample_name}/strelka/results/variants/variants.vcf.gz"
    conda: "../wrappers/strelka/env.yaml"
    script: "../wrappers/strelka/script.py"



# def mpileup_bam_input(wildcards):
#     if config["material"] != "RNA":
#         tag = "bam"
#     else:
#         tag = "RNAsplit.bam"
#     if wildcards.sample_pair == "tumor":
#         return expand("../input_files/mapped/{tumor_bam}.{tag}",tumor_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_tumor"],tag=tag)
#     else:
#         return expand("../input_files/mapped/{normal_bam}.{tag}",normal_bam=sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "sample_name_normal"],tag=tag)
#
# def sample_orig_bam_names(wildcards):
#     if config["calling_type"] == "paired":
#         return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "original_sample_name_tumor"])[0], \
#                 'normal': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "original_sample_name_normal"])[0]}
#     else:
#         return {'tumor': expand("{val}",val = sample_tab.loc[sample_tab.sample_name == wildcards.sample_name, "original_sample_name_tumor"])[0]}



# rule varscan_single:
#     input:
#         unpack(bam_inputs),
#         ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#         regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0],
#     output:
#         vcf="variant_calls/{sample_name}/varscan/VarScan2.vcf",
#     log: "logs/{sample_name}/callers/varscan.log"
#     threads: 1
#     resources:
#         mem_mb=9000
#     params:
#         tumor_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_tumor.mpileup.gz",
#         normal_pileup = "variant_calls/{sample_name}/varscan/{sample_name}_normal.mpileup.gz",
#         snp="variant_calls/{sample_name}/varscan/VarScan2.snp.vcf",
#         indel="variant_calls/{sample_name}/varscan/VarScan2.indel.vcf",
#         extra = config["varscan_extra_params"],
#         # " --strand-filter 0 --p-value 0.95 --min-coverage 50 --min-reads2 8 --min-avg-qual 25 --min-var-freq 0.0005",
#         calling_type = config["calling_type"]
#     conda: "../wrappers/varscan/env.yaml"
#     script: "../wrappers/varscan/script.py"
#
# rule germline_varscan:
#     input:  bam = bam_input,
#             ref=expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
#             regions=expand("{ref_dir}/intervals/{lib_ROI}/{lib_ROI}.bed",ref_dir=reference_directory,lib_ROI=config["lib_ROI"])[0]
#     output: vcf = ADIR+"/varscan/{full_name}.germline.vcf"
#     log:    run = ADIR+"/sample_logs/{full_name}/germline_varscan.log",
#     params: mpileup = ADIR+"/varscan/{full_name}.germline.mpileup",
#             snps = ADIR+"/varscan/{full_name}.germline.snps.vcf",
#             indels = ADIR+"/varscan/{full_name}.germline.indels.vcf"
#     threads: 10
#     conda:  "../wraps/variant_calling/germline_varscan/env.yaml"
#     script: "../wraps/variant_calling/germline_varscan/script.py"