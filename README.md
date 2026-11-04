# Ampseeker2.0
Ampseeker code modified to be integrated into terra. The original workflow guidelines and documentation ca be accessed in https://sanjaynagi.github.io/AmpSeeker/. 

<div align="center">

[<img src="https://raw.githubusercontent.com/sanjaynagi/AmpSeeker/main/docs/ampseeker-docs/logo.png" width="400"/>](https://raw.githubusercontent.com/sanjaynagi/AmpSeeker/main/docs/ampseeker-docs/logo.png)   


</div>

AmpSeeker2.0 is a WDL workflow for Amplicon Sequencing data analysis for both Illumina and Nanopore amplicon sequencing. The pipeline is generic and can work on any data, but has extra modules tailored towards insecticide resistance monitoring. It implements:

- BCL to Fastq conversion
- Genome alignment
- Variant calling
- Quality control
- Coverage
- Visualisation of reads in IGV
- VCF to DataFrame/.xlsx 
- Allele frequency calculation
- Population structure
- Geographic sample maps
- Genetic diversity

- Kdr origin analysis (Ag-vampIR panel)
- Species assignment (Ag-vampIR panel)

The workflow uses a combination of papermill and jupyter book, so that users can visually explore the results in a local webpage for convenience.

## Usage

Ampseeker2.0 is designed to be ran in the terra platform. Please see the [documentation](https://app.terra.bio/#workspaces/malariagen-vector/vobs-eau-claire) for more information on running the workflow.


## Contributing to AmpSeeker

If you would like to contribute to Ampseeker2.0 please follo the instructions in the CONTRIBUTING page.
