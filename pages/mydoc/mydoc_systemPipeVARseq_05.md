---
title: 5. Variant calling
last_updated: Sat Apr 18 12:30:50 2020
sidebar: mydoc_sidebar
permalink: mydoc_systemPipeVARseq_05.html
---

The following performs variant calling with `GATK` and `BCFtools` on a single 
machine by `runCommandline` function can be used to run the variant 
calling with `GATK` and `BCFtools` for each sample sequentially. If a cluster 
is available, running in parallel mode on a compute cluster can be performed by 
`clusterRun` (McKenna et al., 2010; Li , 2011). Typically, the user would choose here 
only one variant caller rather than running several ones.

Not all users have a cluster system, so here to demonstrate an example variant calling 
workflow, only single-machine commands are shown. For cluster jobs, please refer 
to previous steps like code for `BWA` as a template to run on the cluster. 

## Variant calling with `GATK`

The following steps are based on `GATK 4.1.1.0` [Best Practice](https://software.broadinstitute.org/gatk/best-practices/). 
A `targets` file is needed to load samples to a `SYSargs2` intance. There are 10 
individual steps where the user can choose where to jump in and where to skip. 
All scripts are located at `param/cwl/gatk`. `BQSR` (Base Quality Score Recalibration) 
and `VQSR` (Variant Quality Score Recalibration) are very specific 
to a limited species like human, so this workflow does not support these steps. 

### Step1: `fastq` to `ubam`

Convert `fastq` files to `bam` files to prepare for the following step. It is very 
important to specific your sequencing platform, default is `illumina`. User need 
to change `gatk_fastq2ubam.cwl` if the platform is different. Platform information 
is needed for the variant caller in later steps to correct calling parameters.


```r
dir_path <- system.file("extdata/cwl/gatk", package = "systemPipeR")
targets.gatk <- "./results/targetsPE.txt"  ## targets generated from BWA
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_fastq2ubam.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FASTQ_PATH1_", 
    FileName2 = "_FASTQ_PATH2_", SampleName = "_SampleName_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args[1:2], make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "GATK_UBAM", new_col_output_index = 1, 
    overwrite = TRUE)
```

### Step2: Merge `bam` and `ubam`

This step merges a `bam` and `ubam` and creates a third `bam` file that contains 
alignment information and remaining information that was removed by the aligner like `BWA`. 
The removed information is essential for variant statistics calculation. Previous steps are 
recommended, but variant calling can still be performed without these steps.


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_mergebams.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(BWA_SAM = "_bwasam_", GATK_UBAM = "_ubam_", 
    SampleName = "_SampleName_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "GATK_MERGE", new_col_output_index = 1, 
    overwrite = TRUE)
```

### Step3: Sort `bam` files by genomic coordinates

Sort `bam` files by genomic coordinates.


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_sort.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", 
    GATK_MERGE = "_mergebam_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "GATK_SORT", new_col_output_index = 1, 
    overwrite = TRUE)
```

### Step4: Mark duplicates

Mark PCR artifacts in sequencing. A `duplicate_metrics` file will also be produced 
by this step, but will not be used for the next step. This file is just for the user 
to check duplicates status summary.


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_markduplicates.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", 
    GATK_SORT = "_sort_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = c("GATK_MARK", "GATK_MARK_METRICS"), 
    new_col_output_index = c(1, 2), overwrite = TRUE)
```

### Step5: Fixing tags

Takes the `bam` from the last step and calculates the NM, MD, and UQ tags. 
These tags are important for variant calling and filtering. 
This step is recommended but can be skipped.  


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_fixtag.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", 
    GATK_MARK = "_mark_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "GATK_FIXED", new_col_output_index = 1, 
    overwrite = TRUE)
```
  
  Up till this step, sample preprocess is done. All analysis ready `BAM` files and 
  their index `.bai` files are created. Individual and cohort calling by 
  `HaplotypeCaller` is performed from the next step.

### Step6: HaplotypeCaller `gvcf`

The `HaplotypeCaller` is running a **gvcf** mode in this step. G stands for 'genomic'. 
The file not only contains variant sites information but also non-variant sites information; 
thus, at the following step, the cohort caller can use this information to validate the true variants.


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_haplotypecaller.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", 
    GATK_FIXED = "_fixed_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "GVCF", new_col_output_index = 1, overwrite = TRUE)
```

### Step7: Import all `gvcfs`

It is recommended to import all **gvcfs** to a 
[TileDB](https://github.com/Intel-HLS/GenomicsDB/wiki) database for fast cohort 
variant calling at the following step. Note: if you are working with non-diploid data, 
use `CombineGVCFs` function from `GATK` and change the `gvcf_db_folder` parameter 
in `param/cwl/gatk/gatk.yaml` to be your combined **gvcf** file.


```r
# drop all *.g.vcf.gz files into results folder, make sure
# the tbi index is also there.
args <- loadWorkflow(targets = NULL, wf_file = "gatk_genomicsDBImport.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args)
cmdlist(args)
output(args)
args <- runCommandline(args = args, make_bam = FALSE)
```

### Step8: Cohort calling of `gvcf`

Assess variants by information from all gvcfs. A collective vcf called 
`samples.vcf.gz` is created by default naming.


```r
args <- loadWorkflow(targets = NULL, wf_file = "gatk_genotypeGVCFs.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args)
cmdlist(args)
output(args)
args <- runCommandline(args = args, make_bam = FALSE)
```

### Step9: Cohort hard filter variants

VQSR is not included in this workflow. Variants are hard filtered together.
See this [Post](https://gatkforums.broadinstitute.org/gatk/discussion/2806/howto-apply-hard-filters-to-a-call-set) for parameters for hard filtering. Change these settings in `param/cwl/gak/gatk_variantFiltration.sh` if needed.


```r
args <- loadWorkflow(wf_file = "gatk_variantFiltration.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args)
cmdlist(args)
output(args)
args <- runCommandline(args = args, make_bam = FALSE)
```

### Step10: Extract variant

After cohort calling, filtering, all variants for all samples are stored in one big file. 
Extract variants for each sample and save them separately (only variants that have 
passed the filters are stored).


```r
targets.gatk <- "./results/targets_gatk.txt"
args <- loadWorkflow(targets = targets.gatk, wf_file = "gatk_select_variant.cwl", 
    input_file = "gatk.yaml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_"))
cmdlist(args)[1:2]
output(args)[1:2]
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_gatk.txt", 
    step = 1, new_col = "FileName1", new_col_output_index = 1, 
    overwrite = TRUE)
```

## Variant calling with `BCFtools`

The following runs the variant calling with `BCFtools`. This tool takes `BWA` 
aligned `BAM` files, sort, mark duplicates by `samtools` and finally call variants 
by `BCFtools`. 


```r
dir_path <- system.file("extdata/cwl/workflow-bcftools", package = "systemPipeR")
targetsPE <- "./results/targetsPE.txt"
args <- loadWorkflow(targets = targetsPE, wf_file = "workflow_bcftools.cwl", 
    input_file = "bcftools.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(SampleName = "_SampleName_", 
    BWA_SAM = "_SAM_"))
cmdlist(args[1])
output(args[1])
args <- runCommandline(args = args, make_bam = FALSE)
writeTargetsout(x = args, file = "./results/targets_bcf.txt", 
    step = 5, new_col = "FileName1", new_col_output_index = 1, 
    overwrite = TRUE)
```

Variant calling ends here. Downstream analysis starts from the next section.

## Inspect VCF file 

Scripts of downstream analysis are stored in `param/cwl/varseq_downstream`

```r_vcf
dir_path <- system.file("extdata/cwl/varseq", package="systemPipeR")   
```

VCF files can be imported into R with the `readVcf` function. 
Both `VCF` and `VRanges` objects provide convenient data structure for
working with variant data (_e.g._ SNP quality filtering). 


```r
library(VariantAnnotation)
args <- loadWorkflow(targets = "./results/targets_gatk.txt", 
    wf_file = "filter.cwl", input_file = "varseq.yml", dir_path = dir_path)
args <- renderWF(args, inputvars = c(FileName1 = "_FILE1_"))
vcf <- readVcf(infile1(args)[1], "A. thaliana")
vcf
vr <- as(vcf, "VRanges")
vr
```

<br><br><center><a href="mydoc_systemPipeVARseq_04.html"><img src="images/left_arrow.png" alt="Previous page."></a>Previous Page &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; &nbsp; Next Page
<a href="mydoc_systemPipeVARseq_06.html"><img src="images/right_arrow.png" alt="Next page."></a></center>
