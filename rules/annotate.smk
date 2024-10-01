
####################################
# Merge variants from all samples for faster annotation
#
rule merge_variants_in_samples:
    input:  var_tabs = expand("{calling_type}_varcalls/{sample_name}.final_variants.tsv", sample_name = sample_tab.sample_name, calling_type = config["calling_type"]),
    output: tsv_for_vep = "annotate/all_variants.tsv"
    log:    "logs/merge_variants_in_samples.log"
    threads: 20
    resources:
        mem_mb=8000
    conda:  "../wrappers/merge_variants_in_samples/env.yaml"
    script: "../wrappers/merge_variants_in_samples/script.py"

## ANNOTATION of VARIANTS in all SAMPLES
rule variant_annotation:
    input:  tsv_for_vep = "annotate/all_variants.tsv"
    output: annotated = "annotate/all_variants.annotated.tsv"
    log:    "logs/variant_annotation.log"
    threads: 20
    resources:
        mem_mb=8000
    params: ref = config["organism_fasta"],
            vep_dir = config["organism_vep_dir"],
            assembly = config["assembly"],
            organism_name = config["organism"],
            format = config["format"],
            not_use_merged = config["not_use_merged"],
            CADD_DB_SNVs = config["organism_cadd_db_snvs"],
            CADD_DB_indels = config["organism_cadd_db_indels"],
            dir_plugins = config["organism_vep_dir"] + "/VEP_plugins/",
    conda:  "../wrappers/variant_annotation/env.yaml"
    script: "../wrappers/variant_annotation/script.py"

rule custom_annotation:
    input:  annotated = "annotate/all_variants.annotated.tsv",
            format_file = expand(config["tooldir"] + "/variant_calling/{calling_type}_small_var_call_format_files/" + config["format"] + ".txt",calling_type = config["calling_type"])[0],
    output: custom_annotated = "annotate/all_variants.annotated.processed.tsv"
    log:    "logs/custom_annotation.log"
    threads: 10
    resources:
        mem_mb=8000
    params: resources_dir = workflow.basedir + "/resources",
            assembly= config["assembly"],
            format = config["format"],
            custom_DB_folder = config["organism_custom_DB_folder"],
            anno_gtf = config["organism_gtf"],
            isWGS=config["lib_ROI"]
    conda:  "../wrappers/custom_annotation/env.yaml"
    script: "../wrappers/custom_annotation/script.py"
