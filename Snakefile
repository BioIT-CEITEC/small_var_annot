import os
import pandas as pd
import json
from snakemake.utils import min_version

configfile: "config.json"
min_version("5.18.0")
GLOBAL_REF_PATH = config["globalResources"]

# Reference processing
#

if config["lib_ROI"] != "wgs" and config["lib_ROI"] != "rna":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

#### Setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference2.json"),)
reference_dict = json.load(f)
f.close()
config["species_name"] = [organism_name for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]
config["organism"] = config["species_name"].split(" (")[0].lower().replace(" ","_")
if len(config["species_name"].split(" (")) > 1:
    config["species"] = config["species_name"].split(" (")[1].replace(")","")



##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
########################################################################################################################
##### Config processing #####
#conversion from new json
if config["calling_type"] == "somatic":
    sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
    if config["tumor_normal_paired"]:
        sample_tab.drop('sample_name', axis=1, inplace=True)
        sample_tab.rename(columns={'donor': 'sample_name'},inplace=True)
        sample_tab.drop_duplicates(subset="sample_name",inplace=True)
else:
    sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
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

