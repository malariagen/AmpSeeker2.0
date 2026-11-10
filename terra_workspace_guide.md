Ampseeker2.0 is a WDL implementation pipeline from the original Ampseeeker workflow in snakemake (https://github.com/sanjaynagi/AmpSeeker). 
It was built in the terra platform to facilitate the use and access to the pipeline as no bioinformatics skills are required nor access to a HPC environment.

It is designed to process sequencing reads from *Anopheles spp.* samples, sequenced by either Illumina or Nanopore technologies. The pipeline also includes a demultiplexing option although this is limited for illumina sequencing (bcl reads). 

The pipeline aligns the sequencing reads to the reference genome, performs variant calling, calculation of allele frequencies and calculates quality control stats as well as coverage metrics. 

When using the Ag-vampIR panel, it also performs species identification and Kdr origin analysis. This panel is based on  the *An. gambiae* AgamP4 PEST reference genome and targets 90 SNP markers in total. It includes multiple locus associated with insecticide resistance, such as Kdr, Ace1 or Rdl, markers for metabolic resistance, or SNPs that tag selective sweeps prevalent in the genomes of wild specimens from the Anopheles 1000 genomes project (Ag1000G) (https://www.biorxiv.org/content/10.1101/2025.02.14.637727v2). 

The Ampseeker pipeline is versatile so it can also be run with a custom panel of SNPs/markers for other *Anopheles* species. 

The Terra workspace contains 2 Ampseeker workflows and they should be used according to the sequencing technology that originated the reads **Ampseeker2.0_Illumina** or **Ampseeker2.0_Nanopore** .

The required inputs to run either workflow are already populated but please note that the values should be used as a guide only as changes are needed for different analysis such as the input reads, dataset, etc. The required inputs to run either workflow as well as their format are described below:

**Required inputs**

dataset 
(name of the data set to be used, it should be named after the sample_set ID provided by malariaGEN - string)

fastq_r1
(path to the fastq1 files - array)

fastq_r2 - only required for illumina
(path to the fastq2 files - array. If using nanopore reads add an empty array [])

fastq_sample_ids
(sample_ids are described in the samplesheet and metadata files - array)

from_bcl
(set to true if demultiplexing of illumina reads is required - boolean)

illumina_bcl_tars
()

illumina_sample_sheets

ploidy

qc_coverage

qc_multiqc

qc_stats

reference_fasta

sample_min_mean_target_coverage

targets_bed

ampseeker_docker

metadata_file

panel

custom_snpeffdb

reference_gff3

analysis_allele_frequencies

analysis_genetic_diversity

analysis_population_structure

analysis_sample_map

build_jupyter_book

**How to load data into Terra**

In the top right corner of the terra workspace there is an icon with a folder shape. This folder represents the allocated google cloud storage for this particular workspace. Files can be uploaded by clicking in the icon and in the upload button (Terra has access to the local computer files). The results of submitted jobs will also be available in this folder for download. 


**How to submit a job**

1. A Data table should be created (in tsv format) containing all the inputs and respective values (example in here https://app.terra.bio/#workspaces/malariagen-vector/vobs-eau-claire/data). Each column represents one of the workflow inputs. This table should be uploaded into the Data workspace (can be accessed via the link).

2.  Once all required data is uploaded into Terra as well as the data table, a copy of the required workflow should be created (Ampseeker2.0_Illumina or Ampseeker2.0_Nanopore). This can be done by clicking the 3 dots at the bottom of the workflow and selecting duplicate. The workflow should be renamed by adding the sample_set id after the original name. The required inputs should be checked and changed, in the inputs tab, as appropriate (new fastqs, sample_ids, reference, etc.). The option *Run workflow(s) with inputs defined by data table* should be selected and the data table from point 1 can be accessed in the "select data" button. In the output tab select "use defaults" and please make sure all values have the work "output" preceding the input value (e.g., this.output_allele_frequencies_nb). 

3. Once everything is complete click SAVE and Launch. Name your workflow and submit the analysis. The submited analysis will be available in the submission history tab. By clicking the submitted analysis, all steps of the workflow will be described live (click list icon in links column). 

