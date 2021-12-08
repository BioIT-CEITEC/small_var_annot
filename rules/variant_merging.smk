
rule merge_variant_callers:
    input:  vcfs = lambda wildcards: expand("variant_calls/{sample_name}/{variant_caller}/{variant_caller}.norm.vcf",\
                                            sample_name=wildcards.sample_name,\
                                            variant_caller = callers),
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            dict= expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0],
    output: not_filtered_vcf = "merged/{sample_name}.raw_calls.vcf",
            vcf= "merged/{sample_name}.processed.vcf",
            tsv = "merged/{sample_name}.processed.tsv"
    log:    "logs/{sample_name}/merge_variant_callers.log"
    threads: 1
    conda:  "../wrappers/merge_variant_callers/env.yaml"
    script: "../wrappers/merge_variant_callers/script.py"



rule normalize_variants:
    input:  vcf = "variant_calls/{sample_name}/{variant_caller}/{variant_caller}.vcf",
            ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir=reference_directory,ref_name=config["reference"])[0],
            dict= expand("{ref_dir}/seq/{ref_name}.dict",ref_dir=reference_directory,ref_name=config["reference"])[0]
    output: vcf = "variant_calls/{sample_name}/{variant_caller}/{variant_caller}.norm.vcf"
    log:    "logs/{sample_name}/callers/{variant_caller}_normalization.log"
    threads: 1
    conda: "../wrappers/normalize_variants/env.yaml"
    script: "../wrappers/normalize_variants/script.py"
