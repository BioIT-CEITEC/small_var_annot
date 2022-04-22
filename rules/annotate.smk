
####################################
# MERGE VARIANTS
#
rule merge_variants_in_samples:
    input:  var_tabs = expand("merged/{sample_name}.processed.tsv", sample_name = sample_tab.sample_name),
    output: tsv_for_vep = "annotate/all_variants.tsv"
    log:    "logs/merge_variants_in_samples.log"
    threads: 20
    resources:
        mem_mb=8000
    conda:  "../wrappers/merge_variants_in_samples/env.yaml"
    script: "../wrappers/merge_variants_in_samples/script.py"

## ANNOTATION of VARIANTS in SAMPLES
rule variant_annotation:
    input:  tsv_for_vep = "annotate/all_variants.tsv"
    output: annotated = "annotate/all_variants.annotated.tsv"
    log:    "logs/variant_annotation.log"
    threads: 20
    resources:
        mem_mb=8000
    params: ref = expand("{ref_dir}/seq/{ref_name}.fa",ref_dir = reference_directory,ref_name = config["reference"])[0],
            vep_dir = expand("{ref_dir}/annot/vep",ref_dir = reference_directory)[0],
            ref_name = config["reference"],
            organism_name = config["organism"],
            format = config["format"],
            not_use_merged = config["not_use_merged"],
            CADD_DB_SNVs = expand("{ref_dir}/annot/vep/CADD_scores_DB/whole_genome_SNVs.tsv.gz",ref_dir = reference_directory)[0],
            CADD_DB_indels = expand("{ref_dir}/annot/vep/CADD_scores_DB/gnomad.genomes.r3.0.indel.tsv.gz",ref_dir = reference_directory)[0],
            dir_plugins = expand("{ref_dir}/annot/vep/VEP_plugins",ref_dir = reference_directory)[0],
    conda:  "../wrappers/variant_annotation/env.yaml"
    script: "../wrappers/variant_annotation/script.py"

rule custom_annotation:
    input:  annotated = "annotate/all_variants.annotated.tsv",
            format_file = GLOBAL_REF_PATH + "/general/germline_small_var_call_fromat_files/" + config["format"] + ".txt",
    output: custom_annotated = "annotate/all_variants.annotated.processed.tsv"
    log:    "logs/custom_annotation.log"
    threads: 10
    resources:
        mem_mb=8000
    params: resources_dir = workflow.basedir + "/resources",
            reference_name = config["reference"],
            format = config["format"],
            custom_DB_folder = expand("{ref_dir}/annot/custom_new2",ref_dir = reference_directory)[0],
            anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"])[0]
    conda:  "../wrappers/custom_annotation/env.yaml"
    script: "../wrappers/custom_annotation/script.py"