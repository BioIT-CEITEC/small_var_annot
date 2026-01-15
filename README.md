# BioIT-CEITEC/small_var_annot
<p align="justify">
This Snakemake workflow is designed for the annotation of small genetic variants, including both somatic and germline variants. The pipeline processes pre-called small variant files and merges them across all analyzed samples (optional). The variants are then annotated using VEP, with support for adding custom annotations. Finally, a post-processing step transforms the annotated output into the selected format. The final results are provided both as a single aggregated variant table and as separate tables for each analyzed sample.
</p>

## Requirements
- Linux environment
- Snakemake ≥ 5.18.0
- Conda / Mamba
- Python

All remaining dependencies are handled by Snakemake using Conda environments defined per rule.
## Parameters
### Required parameters
The workflow uses a local *config.json* file, which is further extended using shared *BioIT-CEITEC/bioroots_utilities*. The required parameters need to be specified in config file. This workflow is primarily designed so that the configuration is generated and managed via an internal GUI.

- `organism`  
  Target organism used for annotation (e.g. `homo_sapiens`).
  
- `assembly`  
  Genome assembly version.

- `release`  
  Annotation release version.

- `lib_ROI`  
  Indicates whether the input data originate from whole-genome sequencing
  or targeted regions.

- `tumor_normal_paired`  
  Specifies whether tumor–normal paired samples are used.

- `calling_type`  
  Variant calling type, supported values include `somatic` and `germline`.

- `format`  
  Output formatting scheme used for the final variant tables.

### Optional parameters
- `use_cohort_data`  
  If enabled, variants from previous experiments are loaded for cohort-level annotation, requires the file `cohort_data/cohort_variants.tsv`.

- `create_cohort_data`  
  If enabled, cohort-level variant information is generated from the current run
  and stored in `cohort_data/cohort_variants.tsv`.

- `min_variant_frequency`  
  Minimum variant frequency threshold used during variant filtering, default value is 0.

- `not_use_merged`

### Sample parameters
- `sample_name`  
  Sample identification.  

- `entity_name`  
  Entity identification.

- `donor`  
  Used when tumor–normal paired samples are enabled.

## Usage
The workflow is executed using Snakemake and requires a prepared configuration file
and pre-called variant files. From the root directory of the repository, run:

```bash
snakemake --use-conda --cores <N>
```

### Required inputs

- `{calling_type}_varcalls/{sample_name}.final_variants.tsv`

  TSV files containing variant calls produced by an upstream variant calling pipeline.  
  The expected input file structure depends on the selected `calling_type` (e.g. somatic or germline).

### Optional inputs

- `cohort_data/cohort_variants.tsv` 

  If enabled via the `use_cohort_data` parameter, variants from previous runs
  are loaded and used for cohort-level annotation.

## Output
### Main outputs

- `final_variant_table.xlsx`  
  Aggregated and annotated variant table containing all processed variants.

- `final_variant_table.tsv`  
  Tab-delimited version of the aggregated variant table.

- `per_sample_final_var_tabs/`  
  Directory containing per-sample annotated variant tables in Excel format.

### Additional outputs

- `annotate/`  
  Intermediate and final annotation files generated during the annotation step.

- `config.json/`  
  Snapshot of the configuration file used for the run, stored for reproducibility.

- `logs/`  
  Log files produced by individual workflow steps.

### Optional outputs

- `cohort_data/cohort_variants.tsv`  
  Cohort-level variant summary generated when `create_cohort_data` is enabled.

- `mutation_loads.xlsx`  
  Mutation load summary file generated for supported organisms and formats.

## Repository structure
```
.
├── Snakefile                     
├── workflow.config.json            
├── rules/                          
│   ├── annotate.smk         
│   └── variant_postprocessing.smk
├── wrappers/                       
│   ├── custom_annotation/
│   │   ├── custom_annotation.R
│   │   ├── custom_annotation_WGS.R         
│   │   ├── scrip.py
│   │   └── env.yaml
│   ├── merge_variants_in_samples/
│   │   ├── merge_variants_in_samples.R         
│   │   ├── scrip.py
│   │   └── env.yaml
│   ├── process_and_format_annot_variants/
│   │   ├── process_and_format_annot_variants.R         
│   │   ├── process_and_format_annot_variants_WGS.R         
│   │   ├── scrip.py
│   │   └── env.yaml
│   └── variant_annotation/
│   │   ├── formats/        
│   │   ├── dbNSFP_plugin_names.txt.R         
│   │   ├── scrip.py
│   │   └── env.yaml
├── resources/                      
└── README.md
```
