# Contributing to AmpSeeker2.0

Welcome to the AmpSeeker2.0 contributing guide! This document explains how the workflow works, how the WDL/Terra implementation and Docker images fit together, and how to extend the pipeline for new panels and use cases.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Repository Layout](#repository-layout)
3. [Docker Images](#docker-images)
4. [WDL Workflow](#wdl-workflow)
5. [Running on Terra](#running-on-terra)
6. [Analysis Notebooks](#analysis-notebooks)
7. [Adding Support for New Panels](#adding-support-for-new-panels)
8. [Adding a New Analysis Module](#adding-a-new-analysis-module)
9. [Configuration System](#configuration-system)
10. [Snakemake Reference Implementation](#snakemake-reference-implementation)
11. [Development Guidelines](#development-guidelines)
12. [Testing](#testing)
13. [Troubleshooting](#troubleshooting)
14. [Contributing Process](#contributing-process)

---

## Architecture Overview

AmpSeeker2.0 is run as a **WDL workflow on Terra (Google Cloud)**, with every tool and notebook packaged in **Docker images**. The original Snakemake workflow is kept in the repository as the reference implementation (see [Snakemake Reference Implementation](#snakemake-reference-implementation)).

The pipeline is split into:

- **Core functionality:** generic amplicon sequencing analysis
- **Platform-specific processing:** Illumina vs. Nanopore
- **Panel-specific analyses:** custom analyses for specific panels (e.g. Ag-vampIR)

### Key Design Principles

- **Reproducibility:** all analyses are parameterised, version-controlled, and run in pinned container images
- **Modularity:** components can be used independently and new modules easily added
- **Interactive exploration:** Jupyter notebooks provide interactive result exploration
- **Standardised outputs:** consistent directory structure and file formats

### How the pieces fit together

```
Terra data table row  ──►  WDL workflow  ──►  one container per task  ──►  outputs in the workspace bucket
 (inputs, one run)         (task graph)        (ampseeker / bcl2fastq /     (linked back to the data table)
                                                clair3 images)
```

- The **WDL** defines the task graph (what runs, in what order, with how much disk/CPU/memory).
- Each task runs inside a **Docker image**. The main `ampseeker` image contains the command-line tools, the Python stack, and a copy of this whole repository at `/opt/AmpSeeker`, so the notebooks run from inside the image.
- **Terra** supplies the inputs (from a data table row or a JSON file), runs the tasks on Google Cloud VMs, and stores outputs in the workspace bucket.

---

## Repository Layout

```
├── config/                # Configuration files and example metadata
│   ├── config.yaml        # Main configuration file
│   ├── [panel].bed        # Target regions for specific panel
│   └── [panel]-metadata.tsv # Example Sample metadata file
├── docker/                 # Dockerfile(s) for the AmpSeeker images
├── docs/                   # Jupyter Book source (baked into the image)
├── resources/              # Reference files and panel-specific resources
│   └── [panel]/            # Panel-specific resources (e.g. ag-vampir)
├── terra/                  # Terra WDL script that runs ampseeker
├── tests/                  # Test configurations and data
├── workflow/               # Core pipeline components
│   ├── lib/                # Reusable analysis modules
│   ├── envs/               # Conda environment specifications (Snakemake)
│   ├── notebooks/          # Analysis notebooks (used by both Snakemake and WDL)
│   ├── rules/              # Snakemake rule files
│   └── Snakefile           # Snakemake entry point
└── workflows/              # WDL workflow file(s)
```

---

## Docker Images

The WDL uses three images, each passed in as a workflow input so they can be changed without editing the WDL.

| WDL input | Purpose | Notes |
|---|---|---|
| `ampseeker_docker` | Everything except BCL conversion and Clair3 | Built from this repo's Dockerfile |
| `bcl2fastq_docker` | Illumina BCL to FASTQ conversion | Illumina's bcl2fastq is proprietary and cannot be redistributed, so this is a **private** image you build yourself |
| `clair3_docker` | Nanopore variant calling | Upstream image `hkubal/clair3`; pin a specific tag |

### What is in the `ampseeker` image

The image is based on `python:3.12-slim` and contains:

- **apt packages:** samtools, bcftools, tabix, bwa, minimap2, fastqc, fastp, multiqc, OpenJDK
- **Vendored binaries:** `mosdepth`, `snpEff` (from `snpEff_latest_core.zip`), and `fastplong` (`COPY fastplong /usr/local/bin/fastplong`, so the binary must be in the build context)
- **Python stack** (`pip install -r requirements.txt`): papermill, pandas, numpy, plotly, scipy, scikit-learn, xgboost, scikit-allel, pysam, anjl, kaleido, ipyleaflet, openpyxl, jupyter-book, ipykernel and others
- **A Jupyter kernel** named `AmpSeq_python`, registered at build time
- **The whole repository** copied to `/opt/AmpSeeker` (`COPY . /opt/AmpSeeker`)

Because of that last point, the notebooks, `workflow/lib`, `config/config.yaml`, `docs/ampseeker-results` and `resources/ag-vampir` all live **inside the image**. WDL tasks refer to them by their `/opt/AmpSeeker/...` paths rather than taking them as inputs.

> **Important:** changing a notebook, a library module, or a resource file has no effect on Terra until you **rebuild the image, push it under a new tag, and point `ampseeker_docker` at that tag.**

Two things currently happen at runtime inside tasks rather than at build time: `kaleido`/Chrome setup (`kaleido_get_chrome`) for Plotly image export, and the Node.js install that Jupyter Book 2 needs in `build_report`. Tasks therefore need internet access. Moving these into the Dockerfile would make runs faster and more reproducible.

### Where images are stored

Images live in the MalariaGEN registry under the `malariagen-pipelines` path, for example:

```
gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:terrav1.0
```

Docker Hub images (e.g. `catianascimento1/ampseeker:<tag>`) also work in Terra and were used during development.

### Building the image

Build from the **repository root** (the Dockerfile copies the whole repo) and always build for `linux/amd64`, because Terra's VMs are amd64. An image built on an Apple Silicon Mac defaults to `arm64` and will fail on Terra with `exec format error`.

```bash
docker buildx build --platform linux/amd64 \
  -f docker/Dockerfile \
  -t gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:<new-tag> \
  .
```

Keep a `.dockerignore` that excludes large test data, reference indexes and `.git`, otherwise they are copied into the image.

Check the architecture of an existing image with:

```bash
docker image inspect <image> --format '{{.Architecture}}'   # should print amd64
```

### Pushing the image

One-time setup (Google Cloud CLI installed and logged in):

```bash
gcloud auth login
gcloud config set project malariagen-jupyterhub
gcloud auth configure-docker gcr.io
```

Then push:

```bash
docker push gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:<new-tag>
```

To publish an **existing local image under a new name** without rebuilding, retag it first:

```bash
docker tag <image_name> \
  gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:<new-tag>
docker push gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:<new-tag>
```

### Changing the image used by a workflow

1. **Make and test your change** (Dockerfile, notebook, library, resource).
2. **Choose a new tag.** Use an explicit version such as `terrav1.1`, never `:latest`, and never reuse a tag for different content. This keeps runs reproducible and keeps Terra's call caching reliable.
3. **Build and push** as above.
4. **Point the workflow at it**, in whichever place you configure inputs:
   - **Terra data table:** edit the `ampseeker_docker` column of the run's row (or re-upload the TSV).
   - **Terra inputs JSON/UI:** set `AmpSeeker.ampseeker_docker` to the full image path.
   - **WDL default:** update the default value of `ampseeker_docker` in the WDL if you want it to be the new default.
5. **Check permissions.** If the registry is private, the service account that runs your Terra workflows needs read access (for Artifact Registry, the *Artifact Registry Reader* role on the repository).
6. **Resubmit.** If you change the image but keep call caching on, tasks that already completed with the old image tag may be reused. Turn call caching off for the first run with a new image if you need everything to rerun.

The same applies to `bcl2fastq_docker` and `clair3_docker`: change the input, not the WDL.

### Building the bcl2fastq image

bcl2fastq (v2.20) is licensed by Illumina. Build the image yourself from Illumina's installer and keep it in a private registry. Do not publish it publicly. Pass its full path as `bcl2fastq_docker`.

---

## WDL Workflow

### Entry point

The workflow is defined in the WDL file under `workflows/` (workflow name `AmpSeeker`). It uses WDL 1.0 and avoids WDL-1.1-only functions, because Terra's Cromwell does not implement them. For example, the `sample_id` to FASTQ lookup is done in a small Python task (`resolve_fastq_array`) instead of native `Map`/`zip`.

### Main inputs

| Input | Description |
|---|---|
| `dataset` | Dataset name, used in output file names |
| `platform` | `illumina` or `nanopore` |
| `panel` | `ag-vampir` turns on species ID and Kdr-origin analysis; use `custom` otherwise |
| `from_bcl` | `true` to demultiplex from BCL, `false` to start from FASTQs |
| `illumina_bcl_tars`, `illumina_sample_sheets` | BCL run tarball(s) and matching SampleSheet(s); `[]` if `from_bcl` is false |
| `fastq_sample_ids`, `fastq_r1`, `fastq_r2` | Parallel arrays of sample IDs and FASTQs (same order); ignored if `from_bcl` is true. Use `[]` for `fastq_r2` with Nanopore |
| `metadata_file` | `.xlsx`, `.tsv` or `.csv` sample metadata |
| `reference_fasta`, `reference_gff3` | Reference genome and annotation |
| `custom_snpeffdb`, `reference_snpeffdb` | Build a custom snpEff database from the GFF3, or download a named one |
| `targets_bed` | Target regions for the panel |
| `ploidy` | Ploidy used for variant calling |
| `cohort_columns` | Metadata columns used for grouping in analyses (see below) |
| `config_yaml` | Optional config override; defaults to the one baked into the image |
| `qc_coverage`, `qc_multiqc`, `qc_stats` | QC toggles |
| `analysis_population_structure`, `analysis_genetic_diversity`, `analysis_allele_frequencies`, `analysis_sample_map` | Analysis toggles (all default to off) |
| `build_jupyter_book` | Build the final HTML report |
| `ampseeker_docker`, `bcl2fastq_docker`, `clair3_docker` | Container images (see [Docker Images](#docker-images)) |

### Flow

1. **Metadata:** `prepare_metadata` builds `metadata.tsv` and the sample list.
2. **BCL to FASTQ** (if `from_bcl`): `bcl_convert_run` (one per run) then `bcl_convert_gather`, plus `index_read_fastqc`.
3. **FASTQ resolution:** `resolve_fastq_array` aligns FASTQs to the metadata sample order.
4. **Reference prep:** `reference_index`, `bwa_index` or `minimap2_index`, `create_ploidy_file`.
5. **Per sample (scattered):** trimming, alignment, `bam_index`, coverage and stats, then variant calling for targets and amplicons.
6. **Merge and stats:** `merge_vcfs`, `vcf_stats`, `multi_qc`.
7. **Annotation:** `snp_eff_db_download` or `create_custom_snp_eff_db`, then `snp_eff` for targets and amplicons.
8. **Notebooks:** QC, `snp_dataframe`, panel-specific (`species_id`, `kdr_analysis`), then optional analysis notebooks.
9. **Report:** `build_report` assembles every notebook's results into a Jupyter Book.

### Snakemake rule to WDL task mapping

| Snakemake rule file | WDL tasks |
|---|---|
| `metadata.smk` | `prepare_metadata`, `metadata_colours`, `overwrite_metadata_colours` |
| `from-bcl.smk` | `bcl_convert_run`, `bcl_convert_gather`, `index_read_fastqc` |
| `map-call-illumina.smk` | `fastp_illumina`, `bwa_align`, `mpileup_call_amplicons_illumina`, `mpileup_call_targets` |
| `map-call-nanopore.smk` | `fastp_nanopore`, `minimap2_align`, `clair3_call_amplicons_nanopore` |
| `qc.smk` | `mosdepth_coverage`, `bam_stats`, `vcf_stats`, `multi_qc` |
| `snpeff.smk` | `snp_eff_db_download`, `create_custom_snp_eff_db`, `snp_eff` |
| `qc-notebooks.smk` / `analysis.smk` | `run_information`, `run_statistics`, `reads_per_well`, `read_qc`, `coverage`, `sample_quality_control`, `snp_dataframe`, `population_structure`, `sample_map`, `allele_frequencies`, `genetic_diversity` |
| `ag-vampir.smk` | `species_id`, `kdr_analysis` |
| `jupyter-book.smk` | `build_report` (process-toc, process-notebooks and jupyter-book merged) |

### Resources per task

Every task has a `disk_size_gb` input and a `disks: "local-disk ~{disk_size_gb} HDD"` runtime attribute. Defaults range from 20 GB (small tasks) to 500 GB (`bcl_convert_gather`). Tasks running on Terra's default disk run out of space on real data, so **new tasks must declare a disk**. Override a size without editing the WDL by setting `AmpSeeker.<task_name>.disk_size_gb` in the inputs, for example `AmpSeeker.bcl_convert_run.disk_size_gb`.

Tasks do not set `preemptible`, so they run on standard VMs.

### Metadata handling

- **`from_bcl = false`:** metadata comes from `metadata_file`. If it is empty, metadata is built from the sample sheet(s).
- **`from_bcl = true`:** metadata is built from the `[Data]` (or `[BCLConvert_Data]`) section of the sample sheet. The well column (`well` or `Sample_Well`) is split into `well_letter` and `well_number`, and `Sample_Plate`/`plate_name` becomes `plate`. If `metadata_file` is also provided, its extra columns (for example `location`, `taxon`, `latitude`, `longitude`) are merged in by `sample_id`. Columns already present from the sample sheet are not overwritten.
- `plate_info` is turned on automatically when `plate`, `well_letter` and `well_number` are all present.
- `cohort_columns` must exist in the final metadata. Column names are **case-sensitive**.

For Nanopore, the metadata must contain an `fq1` column.

Optional metadata columns: `plate`, `well_letter`, `well_number` (plate layout plots); `latitude`, `longitude` (sample map); any number of grouping columns.

### Outputs

The workflow outputs include the metadata and QC-passed metadata, BAMs/BAIs, annotated VCFs (targets and amplicons), VCF stats, SNP tables (`.xlsx`), the coverage Excel file, the MultiQC report and archive, each analysis notebook, the species-ID and Kdr outputs (Ag-vampIR), and the Jupyter Book (`jupyterbook_index`, `jupyterbook_site` as a `.tar.gz`, and `jupyterbook_standalone` as a single self-contained HTML file). Optional outputs only exist when the corresponding step ran.

---

## Running on Terra

### Option 1: inputs from file paths

Use **"Run workflow with inputs defined by file paths"** and upload an inputs JSON. This is quick for a one-off run, but **outputs are not written back to the data model** in this mode.

### Option 2: inputs from a data table (recommended)

Each row of an `ampseeker_run` table is one run.

1. **Create the table.** Upload a TSV whose first column is `entity:ampseeker_run_id`. The remaining columns are the workflow inputs (`dataset`, `platform`, `panel`, `from_bcl`, `fastq_r1`, `ampseeker_docker`, and so on). Array inputs are written as JSON arrays in a single cell, for example `["gs://bucket/a_1.fastq.gz","gs://bucket/b_1.fastq.gz"]`. Leave optional inputs empty.
2. **Map inputs once.** In the workflow's **Inputs** tab, set each input to `this.<column_name>`, then save. After that, selecting a different row uses that row's values.
3. **Map outputs once.** In the **Outputs** tab, type a **new** column name for each output you want as a link, such as `this.jupyterbook_standalone` or `this.targets_annot_vcf`, and save. Do this **before** launching: mappings only apply to runs submitted after they are saved.
4. **Launch** on the row you want to run.

Do not map an output to a column that already holds an input (for example `this.build_jupyter_book` or `this.clair3_docker`). Terra can pre-fill these incorrectly, and the output would overwrite the input.

### Call caching

Call caching reuses results from earlier runs. Leave it on while iterating on late-stage tasks. Turn it **off** when you change an image, the BCL tarball, the sample sheet, or anything upstream of a task you need to rerun.

### Preparing BCL input

- Tar the run folder from its **parent**, so the archive contains a single top-level folder with `RunInfo.xml` inside it:
  ```bash
  cd /path/to/parent
  tar -czf run.tar.gz run_folder/
  ```
- On macOS add `--no-xattrs` to avoid thousands of extended-attribute warnings.
- Verify before uploading: `tar -tzf run.tar.gz | head`, and check that `Data/Intensities/BaseCalls/L001/C1.1/` exists (watch for typos in folder names).
- `bcl_convert_run` extracts the archive and uses the first folder as the run directory, so extra nesting breaks it.

### Preparing FASTQ input

- `fastq_sample_ids`, `fastq_r1` and `fastq_r2` must be the same length and in the same order.
- Never include index reads (`_I1`/`_I2`) in `fastq_r1`/`fastq_r2`.
- Use plain straight quotes `"` in the Terra input fields. Smart quotes pasted from documents cause errors.

---

## Analysis Notebooks

Jupyter notebooks implement all analyses and are executed with papermill, which injects parameters, so the same notebooks run interactively during development and automatically in the workflow. They are shared by the Snakemake and WDL implementations.

### Quality Control Notebooks
- `read-quality.ipynb`: sequencing quality metrics
- `coverage.ipynb`: coverage analysis across targets
- `sample-quality-control.ipynb`: sample filtering and QC metrics
- `run-statistics.ipynb`: Illumina run statistics (BCL processing)
- `reads-per-well.ipynb`: per-well read distribution (plate-based data)
- `run-information.ipynb`: run and panel overview

### Population Genetics Notebooks
- `population-structure.ipynb`: PCA, neighbor-joining trees
- `genetic-diversity.ipynb`: nucleotide diversity, Tajima's D
- `allele-frequencies.ipynb`: allele frequency calculations
- `sample-map.ipynb`: geographic visualisation (needs `latitude` and `longitude`)

### Utility Notebooks
- `snp-dataframe.ipynb`: convert VCF to DataFrame/Excel
- `process-toc.ipynb`, `process-notebooks.ipynb`: assemble the Jupyter Book
- `IGV-explore.ipynb`: IGV session file generation (not part of the workflow)

### Shared library (`workflow/lib/shared.py`)
- `load_vcf()`: VCF loading with quality filtering
- `create_color_mapping()`: consistent colour schemes for visualisations
- `get_fastqs()`: FASTQ path resolution (Snakemake)
- `AmpSeekerOutputs()`: dynamic output requirement definition (Snakemake)

---

## Adding Support for New Panels

Use the Ag-vampIR panel as a reference implementation. A panel needs resources, a targets BED, notebooks, a WDL task, and report integration.

### Step 1: Panel-specific resources

```
resources/
└── your-panel/
    ├── marker_SNPs.csv      # Panel-specific variant definitions
    ├── reference_data.tsv   # Any reference datasets
    └── ...
```

Resources are copied into the image, so they are available at `/opt/AmpSeeker/resources/your-panel/`.

### Step 2: Target regions file

```
# config/your-panel.bed
chr1    1000    1001    amplicon_1
chr2    5000    5001    amplicon_2
```

The BED file is passed to the workflow as the `targets_bed` input.

### Step 3: Panel-specific notebooks

Create notebooks in `workflow/notebooks/[panel-name]/`, copying how the Ag-vampIR notebooks are structured.

```python
import pandas as pd
import numpy as np
import allel

# Parameters (injected by papermill)
dataset = "your-dataset"
metadata_path = "path/to/metadata"
vcf_path = "path/to/vcf"
panel_resources = "resources/your-panel/"

metadata = pd.read_csv(metadata_path, sep='\t')
vcf = allel.read_vcf(vcf_path)

def your_panel_analysis():
    ...  # implement your analysis
```

The second cell holds the parameters injected by papermill and **must carry the `parameters` tag**. To hide code in the Jupyter Book, add the `remove-input` tag to those cells.

### Step 4: WDL task

Add a task to the WDL, modelled on `species_id` / `kdr_analysis`. Each notebook task follows the same pattern:

```wdl
task your_panel_analysis {
  input {
    Int disk_size_gb = 50
    File targets_vcf
    File metadata_qcpass
    Array[String] cohort_columns
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    # Stage whatever the notebook reads from results/... here

    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/your-panel/your-analysis.ipynb your-analysis.ipynb \
      -k AmpSeq_python \
      -p wkdir "." -p metadata_path ~{metadata_qcpass} -p vcf_path ~{targets_vcf} \
      -p cohort_cols "~{sep="," cohort_columns}"

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf your_panel_results.tar.gz results
  >>>
  output {
    File nb = "your-analysis.ipynb"
    File docs_nb = "your-analysis.ipynb"
    File results_tar = "your_panel_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}
```

Then:

1. Call it from the workflow inside an `if (panel == "your-panel") { ... }` block.
2. Add its `docs_nb` and `results_tar` as optional inputs to `build_report`, copy the notebook into `docs/ampseeker-results/notebooks/your-panel/`, and extract the results tarball into the shared `results/` tree.
3. Add the outputs to the workflow `output {}` block.

### Step 5: Report integration

Update `workflow/notebooks/process-notebooks.ipynb` so the injected papermill parameter cells are hidden for the new notebooks:

```python
import glob

input_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/*ipynb")
if panel == 'ag-vampir':
    other_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/ag-vampir/*ipynb")
    input_nbs.extend(other_nbs)

if panel == 'your-panel':
    your_panel_nbs = glob.glob(f"{wkdir}/docs/ampseeker-results/notebooks/your-panel/*ipynb")
    input_nbs.extend(your_panel_nbs)
```

Also update `workflow/notebooks/process-toc.ipynb` and the table of contents skeleton so the new pages appear in the book.

### Step 6: Rebuild and release

Notebooks and resources live inside the image, so [rebuild and push a new image tag](#changing-the-image-used-by-a-workflow) and update `ampseeker_docker`.

### Snakemake integration (optional)

To keep the Snakemake workflow in step, add `workflow/rules/[panel-name].smk`, include it in `workflow/Snakefile` under `if panel == "your-panel":`, and extend `AmpSeekerOutputs()` in `workflow/rules/common.smk`. See the [Snakemake Reference Implementation](#snakemake-reference-implementation).

---

## Adding a New Analysis Module

1. Create the notebook in `workflow/notebooks/` with a tagged `parameters` cell.
2. Add a WDL task following the pattern above, including `disk_size_gb`.
3. Add a toggle input (for example `analysis_my_module`, default `false`) and call the task inside `if (analysis_my_module) { ... }`.
4. Feed its notebook and results tarball into `build_report`, and add its outputs to the workflow outputs.
5. Rebuild the image, push a new tag, and test on a small dataset.
6. Update the documentation.

---

## Configuration System

- `config/config.yaml` ships inside the image and is used by default. To override it for a run, pass a file as the `config_yaml` input.
- Do not add panel-specific options to `config/config.yaml`. If a panel needs extra options, create a second file (for example `your-panel.config.yaml`) and read it from the panel-specific analyses.
- Most settings that used to live in the config are now workflow inputs in the WDL (toggles, paths, thresholds).

Example Snakemake configuration:

```yaml
# config/your-panel-example-config.yaml
dataset: your-panel-study
panel: your-panel
targets: config/your-panel.bed
metadata: config/your-panel-metadata.tsv

reference-fasta: resources/reference/genome.fasta
reference-gff3: resources/reference/annotations.gff
```

---

## Snakemake Reference Implementation

The original Snakemake workflow is retained in `workflow/` (`Snakefile`, `rules/`, `envs/`). It shares the notebooks and library code with the WDL.

- **Entry point:** `workflow/Snakefile` validates configuration, loads metadata, includes platform- and panel-specific rules, and defines final outputs through `AmpSeekerOutputs()` in `rules/common.smk`.
- **Rule files:** `common.smk`, `metadata.smk`, `qc.smk`, `qc-notebooks.smk`, `from-bcl.smk`, `map-call-illumina.smk`, `map-call-nanopore.smk`, `snpeff.smk`, `analysis.smk`, `jupyter-book.smk`, and `ag-vampir.smk`.
- **Illumina:** raw reads, FastP, BWA-MEM, BCFtools, SnpEff.
- **Nanopore:** raw reads, FastP long, Minimap2, Clair3, SnpEff.
- **Environments:** `AmpSeeker-python.yaml`, `AmpSeeker-cli.yaml`, `AmpSeeker-qc.yaml`, `AmpSeeker-jupyterbook.yaml`, `AmpSeeker-snpeff.yaml` in `workflow/envs/`.

Run it with:

```bash
snakemake --configfile tests/snakemake/config/config_agvampir.yaml --use-conda -j4
```

---

## Development Guidelines

### Code Style
- **Python:** follow PEP 8
- **Plotting:** use plotly with consistent colour schemes
- **WDL:** one task per step, explicit `runtime` (docker, cpu, memory, disks), `set -euo pipefail` in commands, and no WDL-1.1-only functions
- **Notebooks:** include markdown documentation for each analysis step
- **Comments:** document complex logic and parameter choices

### Repository hygiene
- Do not commit large binaries. GitHub rejects files over 100 MB and warns above 50 MB (for example `snpEff_latest_core.zip` or bwa index files such as `*.bwt`, `*.sa`, `*.pac`). Download them in the Dockerfile or test setup, or use Git LFS.
- Keep the Dockerfile and `requirements.txt` pinned enough that rebuilding gives the same tools.

### Version control
- Use semantic versioning and tag releases.
- Keep image tags in step with releases (for example image `terrav1.0` for release 1.0), and record the tag used in the WDL default.
- Test thoroughly before merging to main.

---

## Testing

### Snakemake tests

Located in `tests/snakemake/`:

```
tests/snakemake/
├── config/
│   ├── config_agvampir.yaml      # Ag-vampIR panel test
│   ├── config_fastqauto.yaml     # Auto FASTQ detection test
│   └── config_fastqmetadata.yaml # Manual FASTQ paths test
└── resources/                    # Test data
```

```bash
snakemake --configfile tests/snakemake/config/config_agvampir.yaml --use-conda -j4
```

These are run automatically on pull requests through GitHub Actions.

### WDL / Terra testing

Before merging a change to the WDL or image, run these on a small dataset in Terra:

- **Illumina from FASTQ** with the Ag-vampIR panel
- **Illumina from BCL** (a short run, or a subset)
- **Nanopore** from FASTQ

Check that the workflow succeeds, the Jupyter Book builds, and the outputs appear as links in the data table.

### Adding new tests
1. Create a test configuration in `tests/snakemake/config/`.
2. Download test data from an external resource in the `.github/workflow` setup.
3. Ensure the test completes successfully.
4. Update the GitHub Actions workflow if needed.

---

## Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| Task fails at start, or in `Initializing`, with an image pull error | Wrong image path or tag (check for typos such as a double dot in the tag), or the Terra service account cannot read a private registry |
| `exec format error` | Image was built for `arm64`. Rebuild with `--platform linux/amd64` |
| `no space left on device` | Task disk too small. Raise `AmpSeeker.<task>.disk_size_gb` |
| `RunInfo.xml` not found in `bcl_convert_run` | BCL tarball has extra nesting. Tar from the parent folder so there is a single top-level run folder |
| bcl2fastq cannot find a BCL file | Check folder names in the archive (for example `Intensities` misspelled) |
| fastp "input files don't contain identical amount of reads" | R1/R2 arrays misaligned or contain index reads (`_I1`/`_I2`). Provide only the paired `_1`/`_2` FASTQs, in matching order |
| `AssertionError` about cohort columns | A name in `cohort_columns` is not a column in the metadata. Names are case-sensitive |
| `KeyError: ''` or `IndexError: list index out of range` on `cohort_cols` | `cohort_columns` was empty, or contained only `sample_id` where a real grouping column is needed. Use a real column such as `phenotype` or `location` |
| `sample-map` fails with `KeyError` on `latitude`/`longitude` | Metadata has no coordinates. Provide them through `metadata_file`, or set `analysis_sample_map` to `false` |
| `well_letter` missing in `run-statistics` | Sample sheet uses `Sample_Well`; make sure you are using a WDL version that maps it to `well_letter`/`well_number`, and check sample IDs match between the sample sheet and metadata |
| Terra input values are rejected | Smart quotes (`“ ”`) in the field. Use straight quotes |
| No output links in the data table | Outputs were mapped after launching, mapped to existing input columns, or the run used "file paths" mode. Map outputs to new `this.<name>` columns **before** launching and use the data-table mode |
| pip warning "multiqc requires kaleido==0.2.1" | Harmless. Newer kaleido is installed on purpose for Plotly image export |
| GitHub push rejected for large files | Remove the files from the commit or use Git LFS (see [Repository hygiene](#repository-hygiene)) |

---

## Contributing Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-analysis`)
3. Implement changes following these guidelines
4. If you changed notebooks, `lib/`, resources or the Dockerfile, build and push a new image tag and test the WDL against it
5. Add tests for new functionality
6. Update documentation as needed
7. Submit a pull request with a clear description, including the image tag used for testing
8. Address review comments
9. Merge after approval and testing

Thank you for contributing to AmpSeeker!
