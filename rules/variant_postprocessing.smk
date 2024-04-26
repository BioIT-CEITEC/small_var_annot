
# ####################################
# # AFTER ANNOTATION PROCESSING
# #
def cohort_data_input(wildcards):
    if config["use_cohort_data"] == True:
        return "cohort_data/cohort_variants.tsv"
    else:
        return []

rule process_and_format_annot_variants:
    input:  var_tabs = expand("{calling_type}_varcalls/{sample_name}.final_variants.tsv", sample_name = sample_tab.sample_name, calling_type = config["calling_type"]),
            annotated = "annotate/all_variants.annotated.processed.tsv",
            format_file = expand(config["tooldir"] + "/{calling_type}_small_var_call_format_files/" + config["format"] + ".txt",calling_type = config["calling_type"])[0],
            cohort_data = cohort_data_input
    output: all_vars_xlsx = "final_variant_table.xlsx",
            all_vars_tsv = "final_variant_table.tsv",
            per_sample_var_tabs = expand("per_sample_final_var_tabs/{sample_name}.variants.xlsx", sample_name = sample_tab.sample_name),
    log:    "logs/postprocess_and_format_annot_variants.log"
    threads: 10
    resources:
        mem_mb=8000
    params: reference = config["reference"],
            min_variant_frequency = str(config["min_variant_frequency"]),
            format = config["format"],
            anno_gtf = config["organism_gtf"],
            create_cohort_data = config["create_cohort_data"],
            batch_name = config["entity_name"],
            ref_dir= reference_directory,
            organism=config["organism"],
            mut_load_output_filename= "mutation_loads.xlsx",
            isWGS=config["lib_ROI"]
    conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
    script: "../wrappers/process_and_format_annot_variants/script.py"


# rule cohort_related_processing:
#     input:  var_tabs = expand("somatic_seq_results/{sample_name}.variants.tsv", sample_name = sample_tab.sample_name),
#             annotated = "annotate/all_variants.annotated.tsv"
#     output: all_vars = "cohort_data"
#     log:    "logs/cohort_related_processing.log"
#     threads: 10
#     params: reference = config["reference"],
#             min_variant_frequency = str(config["min_variant_frequency"]),
#             format = config["format"],
#             anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"])
#     conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
#     script: "../wrappers/process_and_format_annot_variants/script.py"
#
#
# rule cohort_sample_check:
#     input:  var_tabs = expand("somatic_seq_results/{sample_name}.variants.tsv", sample_name = sample_tab.sample_name),
#             annotated = "annotate/all_variants.annotated.tsv"
#     output: all_vars = "full_variant_table.xlsx"
#     log:    "logs/postprocess_and_format_annot_variants.log"
#     threads: 10
#     params: reference = config["reference"],
#             min_variant_frequency = str(config["min_variant_frequency"]),
#             format = config["format"],
#             anno_gtf = expand("{ref_dir}/annot/{ref_name}.gtf",ref_dir = reference_directory,ref_name = config["reference"])
#     conda:  "../wrappers/process_and_format_annot_variants/env.yaml"
#     script: "../wrappers/process_and_format_annot_variants/script.py"