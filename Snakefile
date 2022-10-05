import os
import pandas as pd
import json
from snakemake.utils import min_version

configfile: "config.json"

min_version("5.18.0")
GLOBAL_REF_PATH = config["globalResources"]

# Reference processing
#
if config["lib_ROI"] != "wgs" and config["lib_ROI"] != "RNA":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]
else:
    config["lib_ROI"] = "wgs"

# setting organism from reference
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
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")

if "tumor_normal" in sample_tab:
    sample_tab.drop(sample_tab[sample_tab.tumor_normal == "normal"].index)

if "donor" in sample_tab:
    sample_tab["sample_name"] = sample_tab["donor"]

if config["calling_type"] == "germline":
    config["format"] = config["germline_format"]
else:
    config["format"] = config["somatic_format"]


# if not config["is_paired"]:
#     read_pair_tags = [""]
#     paired = "SE"
# else:
#     read_pair_tags = ["_R1","_R2"]
#     paired = "PE"



# callers = config["callers"].split(';')


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

