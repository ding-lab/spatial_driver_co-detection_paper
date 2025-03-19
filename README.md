# spatial_driver_co-detection_paper
This is the repository for the spatial driver detection paper: "Spatial Co-mapping of Cancer Mutations and Gene Expression at Subcellular Resolution".

Probes were designed using the code found at `Variant_site_flanking_sequence_and_probe_compatibility_scripts/Script_for_pulling_SNP_flanking_sites_v6.R`
* This code runs using the packages installed in the r4.2 environment yaml file: `spatial_driver_co-detection_paper/r4.2_env_used_to_generate_probe_flanks_last_update_20240201.yml`
    * it relies on the following packages specifically:
    ```
    r-base
    biomaRt
    BSgenome.Hsapiens.UCSC.hg38
    GenomicFeatures
    GenomicRanges
    stringi
    ```
* Input targeted sequences are provided at `spatial_driver_co-detection_paper/Variant_site_flanking_sequence_and_probe_compatibility_scripts/input_variant_targets.tsv`
* Output targeted sequences are provided at `spatial_driver_co-detection_paper/Variant_site_flanking_sequence_and_probe_compatibility_scripts/input_variant_targets_flanks.tsv`
* Flank sequences require manual review and curration as outlined in the methods. that is present here: `spatial_driver_co-detection_paper/Variant_site_flanking_sequence_and_probe_compatibility_scripts/manual_indel_flank_site_check.xlsx`

Analysis for this paper relied upon the packages installed in the r4.3.2 environment yaml file: `spatial_driver_co-detection_paper/r4.3.2_seurat5.0.1_for_analysis_20250116.yml`
* Installation can be done via `conda env create -n scrublet --file r4.3.2_seurat5.0.1_for_analysis_20250116.yml`
    * analysis relied upon the following packages:
    ```
    Seurat_5.0.1
    SeuratObject_5.0.1
    arrow_14.0.1
    future_1.33.1
    optparse_1.7.4
    magrittr_2.0.3
    tidyverse_2.0.0
    epitools_0.5-10.1
    effectsize_0.8.9
    ggrepel_0.9.5
    ggrastr_1.0.2
    ggpubr_0.6.0
    ComplexHeatmap_2.18.0
    circlize_0.4.16
    RColorBrewer_1.1-3
    grid_4.3.2
    Matrix_1.6-5
    viridisLite_0.4.2
    patchwork_1.2.0
    scales_1.3.0
    scCustomize_3.0.0
    enrichR_3.2
    reshape2_1.4.4
    ggpmisc_0.5.5
    ggh4x_0.2.8
    ```
    * the libaries ComplexHeatmap and Circlize were installed directly in R with BiocManager rather than conda.
* Code related to analysis from the paper can be found in the `spatial_driver_co-detection_paper/plotting_snippets/` folder. The worklogs in each folder describe how the code/tables were used and how code was run.

## Other tools used in the analysis of this paper include:
Copy number calling was run with [GATK4SCNA](https://github.com/Aust1nS2/GATK4SCNA).
* This relies on two docker image: `austins2/gatk4scna:v1.1` and `austins2/ggplot_gatk4scna:v.2024.08.19`.
* It is built for running on a HPC cluster running the IBM LSF job scheduler.
* The worklog for how it was run on the data for this publication can be found at `spatial_driver_co-detection_paper/bulk_WES_cnv/GATK4_somatic_cnv_detection_worklog.sh`

Scrublet for new snRNA-seq/multiome samples was run with [automated_scrublet](https://github.com/Aust1nS2/automated_scrublet).
* the environment for scrublet automation can be found on the scrublet repository. A docker image is also provided.
* Example for the stand alone and docker image can be found there.
* The worklog for running scrublet on the new snRNA-seq data present in this cohort can be found at `spatial_driver_co-detection_paper/snRNA_scRNA_preprocessing/object_generation_worklog.sh`
* Some additional python scripts used the environment found on the scrublet repository. Where applicable this is noted in the worklog.

H&E images were aligned with [HEX-SIFT](https://github.com/Aust1nS2/HEX-SIFT) as outlined in the worklog file.
* The envrionment and installation instructions for HEX-SIFT image alignment can be found in the above HEX-SIFT repository.
* The worklog for running HEX-SIFT on each sample can be found at `spatial_driver_co-detection_paper/HE_image_alignment/HEX-SIFT_alignment_worklog.sh`
* Depending on the input image file this can take up to 400GB of memory.

Mutation mapping in matching snRNA-seq data relied on [scVarScan](https://github.com/ding-lab/10Xmapping) and custom scripts for post-processing as outlined in the worklog file.
* The version of perl installed on the server that scVarScan was run on was v5.26.2. The post-processing scripts run with the conda python 3.9.6 environment from [automated_scrublet](https://github.com/Aust1nS2/automated_scrublet) active.
* The worklog for running scVarScan can be found at `spatial_driver_co-detection_paper/snRNA_mutation_mapping/matching_snRNA_mutation_mapping_worklog.sh`
* post processing scripts are located in `spatial_driver_co-detection_paper/snRNA_mutation_mapping`

Variant calling (snvs/indels) was run using [Somaticwrapper](https://github.com/ding-lab/somaticwrapper) and it's implementation in the [PECGS pipeline](https://github.com/ding-lab/pecgs-pipeline).

Layer calculation was run using [Morph](https://github.com/ding-lab/morph)
* environment for this can be found at `spatial_driver_co-detection_paper/morph_layers/`. Alternatively a docker image has been made available at: [austins2/spatial_driver_layer](https://hub.docker.com/repository/docker/austins2/spatial_driver_layer/general).
* scripts and relevant input tables are available at `spatial_driver_co-detection_paper/morph_layer_calculations/spatial_driver_layer.yml`

## Data Acccess
