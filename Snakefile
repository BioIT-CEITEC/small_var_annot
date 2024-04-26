import os
import pandas as pd
import json
from snakemake.utils import min_version

configfile: "config.json"
min_version("5.18.0")
GLOBAL_REF_PATH = config["globalResources"]


##### BioRoot utilities #####
module BR:
    snakefile: gitlab("bioroots/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*

tool_dir = BR.load_tooldir()

config = BR.load_organism()

########################################################################################################################
##### Config processing #####
#conversion from new json
if config["calling_type"] == "somatic":
    sample_tab = BR.load_sample()
    if config["tumor_normal_paired"]:
        sample_tab.drop('sample_name', axis=1, inplace=True)
        sample_tab.rename(columns={'donor': 'sample_name'},inplace=True)
        sample_tab.drop_duplicates(subset="sample_name",inplace=True)
else:
    sample_tab = BR.load_sample()
    config["format"] = config["germline_format"]

#######################################################################################################################

# DEFAULT VALUES
if not "format" in config:
    config["format"] = "default"
if not "not_use_merged" in config:
    config["not_use_merged"] = False
if not "min_variant_frequency" in config:
    config["min_variant_frequency"] = 0

wildcard_constraints:
    vartype = "snvs|indels",
    sample = "|".join(sample_tab.sample_name),
    lib_name = "[^\.\/]+",
    read_pair_tag = "(_R.)?"


####################################
# SEPARATE RULES
include: "rules/annotate.smk"
include: "rules/variant_postprocessing.smk"

####################################
# RULE ALL
rule all:
    input:
        all_vars_xlsx = "final_variant_table.xlsx",
        all_vars_tsv = "final_variant_table.tsv",
        per_sample_var_tabs = expand("per_sample_final_var_tabs/{sample_name}.variants.xlsx", sample_name = sample_tab.sample_name),

