version 1.0

## =========================================================================
## AmpSeeker (https://github.com/sanjaynagi/AmpSeeker) — WDL translation
## Version: 2.0
## =========================================================================

workflow AmpSeeker {
  input {
    # ---- containers --------------------------------------
    String ampseeker_docker = "gcr.io/malariagen-jupyterhub/malariagen-pipelines/ampseeker:terrav1.0"
    String bcl2fastq_docker = "gcr.io/malariagen-jupyterhub/malariagen-pipelines/bcl2fastq:v2.20"
    String clair3_docker    = "hkubal/clair3:latest"
    # ---- run identity -----------------------------------------------------
    String dataset
    String panel = "custom"          # "ag-vampir" turns on species-id/kdr-origin modules
    String platform                  # "illumina" | "nanopore"
    Array[String] cohort_columns = ["location"]
    # ---- sample input -------------------------------------------------------
    # No defaults on any of these — every run requires the user to pick them
    # explicitly in Terra's Inputs UI, pointing at that run's own GCS paths.
    Boolean from_bcl
    Array[File] illumina_bcl_tars       # pass [] if from_bcl = false
    Array[File] illumina_sample_sheets  # pass [] if from_bcl = false
    File? metadata_file
    # Data-table friendly FASTQ input: three parallel arrays, same row order
    # as a Terra sample_set (this.samples.sample_id / .fq1 / .fq2). sample_ids
    # here does NOT need to be the full cohort — it only has to align 1:1
    # with fastq_r1/fastq_r2. Ignored when from_bcl = true (pass [] then).
    Array[String] fastq_sample_ids
    Array[File] fastq_r1
    Array[File] fastq_r2   # illumina only — pass [] for nanopore
    # ---- reference ----------------------------------------------------------
    File reference_fasta
    File? reference_gff3
    String reference_snpeffdb = ""
    Boolean custom_snpeffdb = false
    File targets_bed                    # config["targets"], e.g. config/ag-vampir.bed
    Int ploidy
    # Optional override for the contig list. Left empty ([]) by default, in
    # which case run-information.ipynb and coverage.ipynb keep deriving it
    # themselves via `panel_metadata.contig.unique()` from targets_bed, same
    # as before. Set explicitly to pin/restrict the contig list instead of
    # relying on auto-derivation from the bed file.
    # NOTE: threading this value through requires the corresponding
    # papermill parameter cell in run-information.ipynb / coverage.ipynb to
    # actually read and use a `contigs` parameter (falling back to the old
    # auto-derive behavior when it's empty) — that notebook-side change is
    # not included here, since the notebook source isn't visible from the
    # WDL. Without it, this input is threaded through but has no effect.
    Array[String] contigs = []
    File? config_yaml                   # config/config.yaml, defaults to the one baked in the image
    # ---- nanopore-only ------------------------------------------------------
    File? clair3_model_tar
    # ---- QC toggles -----------------------------------------------------------
    Boolean qc_coverage
    Boolean qc_multiqc
    Boolean qc_stats
    Float sample_min_mean_target_coverage
    # ---- analysis toggles -------------------------------------------------------
    # config/config.yaml ships with all four of these OFF by default (only
    # `panel == "ag-vampir"` deployments turn most of them on) — defaulted
    # to match, override per-run as needed.
    Boolean analysis_population_structure = false
    Boolean analysis_genetic_diversity = false
    Boolean analysis_allele_frequencies = false
    Boolean analysis_sample_map = false
    Boolean build_jupyter_book = true   # matches config.yaml's build-jupyter-book: true
    File? kdr_marker_snps_csv           # override for baked-in resources/ag-vampir/Kdr_marker_SNPs.csv
    File? user_metadata_colours_json    # override for a computed config/metadata_colours.json
  }

  # ==========================================================================
  # 1. METADATA
  # ==========================================================================
  call prepare_metadata {
    input: from_bcl = from_bcl, sample_sheets = illumina_sample_sheets, metadata_file = metadata_file,
           cohort_columns = cohort_columns, platform = platform, docker = ampseeker_docker,
  }
  Array[String] sample_ids = read_lines(prepare_metadata.sample_ids_txt)
  Boolean plate_info        = read_boolean(prepare_metadata.plate_info_txt)

  # ==========================================================================
  # 2. BCL -> FASTQ
  # ==========================================================================
  if (from_bcl && platform == "illumina") {
    scatter (i in range(length(illumina_bcl_tars))) {
      call bcl_convert_run { input: bcl_tar = illumina_bcl_tars[i], sample_sheet = illumina_sample_sheets[i], docker = bcl2fastq_docker }
    }
    call bcl_convert_gather {
      input: run_output_tars = bcl_convert_run.run_output_tar,
             demultiplex_stats = bcl_convert_run.demultiplex_stats, sample_ids = sample_ids, docker = ampseeker_docker,
    }
    call index_read_fastqc {
      input: i1 = bcl_convert_gather.index1_reads, i2 = bcl_convert_gather.index2_reads, docker = ampseeker_docker
    }
  }

  # No native WDL Map/zip/as_map here (Terra's Cromwell doesn't implement the
  # WDL-1.1 collection functions) — the sample_id -> FASTQ lookup is done in
  # a plain Python task instead, against whichever source is active
  # (BCL-demultiplexed output, or the user-supplied fastq_sample_ids/
  # fastq_r1/fastq_r2 arrays).
  Array[String] fq_source_sample_ids = select_first([bcl_convert_gather.sample_order, fastq_sample_ids])
  Array[File]   fq_source_r1         = select_first([bcl_convert_gather.fastq_r1_sorted, fastq_r1])

  call resolve_fastq_array as resolve_r1 {
    input: sample_ids = sample_ids, avail_sample_ids = fq_source_sample_ids, avail_files = fq_source_r1, docker = ampseeker_docker
  }
  if (platform == "illumina") {
    Array[File] fq_source_r2 = select_first([bcl_convert_gather.fastq_r2_sorted, fastq_r2])
    call resolve_fastq_array as resolve_r2 {
      input: sample_ids = sample_ids, avail_sample_ids = fq_source_sample_ids, avail_files = fq_source_r2, docker = ampseeker_docker
    }
  }

  # ==========================================================================
  # 3. REFERENCE PREP
  # ==========================================================================
  call reference_index { input: reference_fasta = reference_fasta, docker = ampseeker_docker }
  if (platform == "illumina") { call bwa_index { input: reference_fasta = reference_fasta, docker = ampseeker_docker } }
  if (platform == "nanopore") { call minimap2_index { input: reference_fasta = reference_fasta, docker = ampseeker_docker } }
  call create_ploidy_file { input: ploidy = ploidy, docker = ampseeker_docker }

  # ==========================================================================
  # 4. PER-SAMPLE
  # ==========================================================================
  scatter (i in range(length(sample_ids))) {
    String sid = sample_ids[i]

    if (platform == "illumina") {
      call fastp_illumina {
        input: sample_id = sid, fq1 = resolve_r1.ordered[i], fq2 = select_first([resolve_r2.ordered])[i], docker = ampseeker_docker
      }
      call bwa_align {
        input: sample_id = sid, fq1 = select_first([fastp_illumina.trimmed_r1]),
               fq2 = select_first([fastp_illumina.trimmed_r2]), reference_fasta = reference_fasta,
               bwa_index_tar = select_first([bwa_index.index_tar]), docker = ampseeker_docker,
      }
    }
    if (platform == "nanopore") {
      call fastp_nanopore { input: sample_id = sid, fq1 = resolve_r1.ordered[i], docker = ampseeker_docker }
      call minimap2_align {
        input: sample_id = sid, reads = select_first([fastp_nanopore.trimmed]),
               reference_fasta = reference_fasta, reference_mmi = select_first([minimap2_index.mmi]), docker = ampseeker_docker,
      }
    }
    File sample_bam = select_first([bwa_align.bam, minimap2_align.bam])
    File fastp_json = select_first([fastp_illumina.json, fastp_nanopore.json])
    File fastp_html = select_first([fastp_illumina.html, fastp_nanopore.html])

    call bam_index { input: sample_id = sid, bam = sample_bam, docker = ampseeker_docker }

    if (qc_coverage) {
      call mosdepth_coverage {
        input: sample_id = sid, bam = sample_bam, bam_index = bam_index.bai, panel_bed = targets_bed, docker = ampseeker_docker
      }
    }
    if (qc_stats) {
      call bam_stats { input: sample_id = sid, bam = sample_bam, bam_index = bam_index.bai, docker = ampseeker_docker }
    }

    call mpileup_call_targets {
      input: sample_id = sid, bam = sample_bam, bam_index = bam_index.bai, reference_fasta = reference_fasta,
             reference_fai = reference_index.fai, ploidy_file = create_ploidy_file.ploidy_file,
             targets_bed = targets_bed, platform = platform, docker = ampseeker_docker,
    }
    call bgzip_tabix as bgzip_tabix_targets { input: vcf = mpileup_call_targets.calls_vcf, docker = ampseeker_docker }

    if (platform == "illumina") {
      call mpileup_call_amplicons_illumina {
        input: sample_id = sid, bam = sample_bam, bam_index = bam_index.bai,
               reference_fasta = reference_fasta, ploidy_file = create_ploidy_file.ploidy_file, docker = ampseeker_docker,
      }
    }
    if (platform == "nanopore") {
      call clair3_call_amplicons_nanopore {
        input: sample_id = sid, bam = sample_bam, bam_index = bam_index.bai, reference_fasta = reference_fasta,
               reference_fai = reference_index.fai, model_tar = clair3_model_tar, ploidy = ploidy,
               docker = clair3_docker,
      }
    }
    File amplicons_calls_vcf = select_first([
      mpileup_call_amplicons_illumina.calls_vcf, clair3_call_amplicons_nanopore.calls_vcf
    ])
    call bgzip_tabix as bgzip_tabix_amplicons { input: vcf = amplicons_calls_vcf, docker = ampseeker_docker }
  }

  # ==========================================================================
  # 5. MERGE + STATS
  # ==========================================================================
  call merge_vcfs as merge_vcfs_targets {
    input: call_type = "targets", dataset = dataset, vcfs = bgzip_tabix_targets.vcfgz, tbis = bgzip_tabix_targets.tbi,
           docker = ampseeker_docker,
  }
  call merge_vcfs as merge_vcfs_amplicons {
    input: call_type = "amplicons", dataset = dataset, vcfs = bgzip_tabix_amplicons.vcfgz, tbis = bgzip_tabix_amplicons.tbi,
           docker = ampseeker_docker,
  }
  if (qc_stats) {
    call vcf_stats as vcf_stats_targets {
      input: call_type = "targets", dataset = dataset, merged_vcf = merge_vcfs_targets.merged_vcf, docker = ampseeker_docker
    }
    call vcf_stats as vcf_stats_amplicons {
      input: call_type = "amplicons", dataset = dataset, merged_vcf = merge_vcfs_amplicons.merged_vcf, docker = ampseeker_docker
    }
  }
  if (qc_multiqc) {
    call multi_qc {
      input: fastp_jsons = fastp_json, flagstats = select_all(bam_stats.flagstat),
             vcf_stats_files = select_all([vcf_stats_targets.stats, vcf_stats_amplicons.stats]),
             per_base_beds = select_all(mosdepth_coverage.per_base_bed),
             summary_files = select_all(mosdepth_coverage.summary), dist_files = select_all(mosdepth_coverage.global_dist),
             docker = ampseeker_docker,
    }
  }

  # ==========================================================================
  # 6. snpEff
  # ==========================================================================
  if (!custom_snpeffdb) { call snp_eff_db_download { input: reference_snpeffdb = reference_snpeffdb, docker = ampseeker_docker } }
  if (custom_snpeffdb) {
    call create_custom_snp_eff_db {
      input: reference_fasta = reference_fasta, reference_gff3 = select_first([reference_gff3]), docker = ampseeker_docker
    }
  }
  Array[File] snpeff_db_ready = select_all([snp_eff_db_download.done_marker, create_custom_snp_eff_db.db_tar])

  call snp_eff as snp_eff_targets {
    input: call_type = "targets", dataset = dataset, merged_vcf = merge_vcfs_targets.merged_vcf,
           custom_snpeffdb = custom_snpeffdb, reference_snpeffdb = reference_snpeffdb, db_ready = snpeff_db_ready[0],
           docker = ampseeker_docker,
  }
  call snp_eff as snp_eff_amplicons {
    input: call_type = "amplicons", dataset = dataset, merged_vcf = merge_vcfs_amplicons.merged_vcf,
           custom_snpeffdb = custom_snpeffdb, reference_snpeffdb = reference_snpeffdb, db_ready = snpeff_db_ready[0],
           docker = ampseeker_docker,
  }

  # ==========================================================================
  # 7. TOP-LEVEL METADATA COLOURS
  # ==========================================================================
  call metadata_colours {
    input: metadata_tsv = prepare_metadata.metadata_tsv, cohort_columns = cohort_columns,
           user_colours_json = user_metadata_colours_json, docker = ampseeker_docker,
  }

  # ==========================================================================
  # 8. QC NOTEBOOKS
  # ==========================================================================
  call run_information {
    input: metadata = prepare_metadata.metadata_tsv, targets_bed = targets_bed, panel = panel, dataset = dataset,
           platform = platform, cohort_columns = cohort_columns, contigs = contigs, config_yaml = config_yaml,
           docker = ampseeker_docker,
  }

  if (from_bcl && platform == "illumina") {
    call run_statistics {
      input: metadata = prepare_metadata.metadata_tsv, stats_json = select_first([bcl_convert_gather.stats_json]),
             colours_json = metadata_colours.colours_json, cohort_columns = cohort_columns, plate_info = plate_info,
             docker = ampseeker_docker,
    }
  }

  if (plate_info) {
    call reads_per_well {
      input: sample_ids = sample_ids, bams = sample_bam, flagstats = select_all(bam_stats.flagstat),
             metadata = prepare_metadata.metadata_tsv, docker = ampseeker_docker,
    }
  }

  call read_qc {
    input: metadata = prepare_metadata.metadata_tsv, index_read_qc = (from_bcl && platform == "illumina"),
           docker = ampseeker_docker,
  }

  if (qc_coverage) {
    call coverage {
      input: sample_ids = sample_ids, per_base_beds = select_all(mosdepth_coverage.per_base_bed),
             regions_beds = select_all(mosdepth_coverage.regions_bed), metadata = prepare_metadata.metadata_tsv,
             targets_bed = targets_bed, contigs = contigs, docker = ampseeker_docker,
    }
  }

  call sample_quality_control {
    input: sample_ids = sample_ids, regions_beds = select_all(mosdepth_coverage.regions_bed),
           targets_vcf = snp_eff_targets.annot_vcf, amplicons_vcf = snp_eff_amplicons.annot_vcf,
           metadata = prepare_metadata.metadata_tsv, targets_bed = targets_bed, cohort_columns = cohort_columns,
           sample_min_mean_target_coverage = sample_min_mean_target_coverage, panel = panel, dataset = dataset,
           platform = platform, colours_json = metadata_colours.colours_json, docker = ampseeker_docker,
  }

  # ==========================================================================
  # 9. AG-VAMPIR: species_id -> overwrite_metadata_colours -> kdr_analysis
  # ==========================================================================
  call snp_dataframe {
    input: targets_vcf = snp_eff_targets.annot_vcf, amplicons_vcf = snp_eff_amplicons.annot_vcf,
           dataset = dataset, platform = platform, docker = ampseeker_docker,
  }

  if (panel == "ag-vampir") {
    call species_id {
      input: targets_vcf = snp_eff_targets.annot_vcf, targets_bed = targets_bed,
             metadata_qcpass = sample_quality_control.qcpass_metadata, cohort_columns = cohort_columns,
             colours_json = metadata_colours.colours_json, docker = ampseeker_docker,
    }

    call overwrite_metadata_colours {
      input: metadata_qcpass = species_id.updated_metadata, cohort_columns = cohort_columns,
             user_colours_json = user_metadata_colours_json, docker = ampseeker_docker,
    }

    call kdr_analysis {
      input: targets_vcf = snp_eff_targets.annot_vcf, amplicons_vcf = snp_eff_amplicons.annot_vcf,
             metadata_qcpass = species_id.updated_metadata, kdr_marker_snps_csv = kdr_marker_snps_csv,
             targets_bed = targets_bed, snp_dataframe_targets_xlsx = snp_dataframe.targets_xlsx,
             colours_json = overwrite_metadata_colours.colours_json, dataset = dataset,
             cohort_columns = cohort_columns, docker = ampseeker_docker,
    }
  }

  File final_metadata = select_first([species_id.updated_metadata, sample_quality_control.qcpass_metadata])
  File final_colours   = select_first([overwrite_metadata_colours.colours_json, metadata_colours.colours_json])

  # ==========================================================================
  # 10. ANALYSIS NOTEBOOKS
  # ==========================================================================
  if (analysis_population_structure) {
    call population_structure {
      input: amplicons_vcf = snp_eff_amplicons.annot_vcf, metadata_qcpass = final_metadata, dataset = dataset,
             cohort_columns = cohort_columns, platform = platform, colours_json = final_colours, docker = ampseeker_docker,
    }
  }

  if (analysis_sample_map) {
    call sample_map { input: metadata_qcpass = final_metadata, dataset = dataset, docker = ampseeker_docker }
  }

  if (analysis_allele_frequencies) {
    call allele_frequencies {
      input: metadata_qcpass = final_metadata, targets_bed = targets_bed, targets_vcf = snp_eff_targets.annot_vcf,
             amplicons_vcf = snp_eff_amplicons.annot_vcf, dataset = dataset, cohort_columns = cohort_columns,
             platform = platform, docker = ampseeker_docker,
    }
  }

  if (analysis_genetic_diversity) {
    call genetic_diversity {
      input: amplicons_vcf = snp_eff_amplicons.annot_vcf, metadata_qcpass = final_metadata, dataset = dataset,
             cohort_columns = cohort_columns, platform = platform, docker = ampseeker_docker,
    }
  }

  # ==========================================================================
  # 11. REPORT (process-toc + process-notebooks + jupyter-book, merged)
  # ==========================================================================
  if (build_jupyter_book) {
    call build_report {
      input:
        from_bcl = from_bcl, panel = panel, plate_info = plate_info, dataset = dataset,
        docs_run_info = run_information.docs_nb, docs_run_statistics = run_statistics.docs_nb,
        docs_snp_df = snp_dataframe.docs_nb, docs_sample_quality_control = sample_quality_control.docs_nb,
        docs_coverage = coverage.docs_nb, docs_pop_structure = population_structure.docs_nb,
        docs_af = allele_frequencies.docs_nb, docs_sample_map = sample_map.docs_nb,
        docs_genetic_diversity = genetic_diversity.docs_nb, docs_read_quality = read_qc.docs_nb,
        docs_reads_per_well = reads_per_well.docs_nb,
        docs_species_id = species_id.docs_nb, docs_kdr_analysis = kdr_analysis.docs_nb,
        run_information_results_tar = run_information.results_tar,
        run_statistics_results_tar = run_statistics.results_tar,
        reads_per_well_results_tar = reads_per_well.results_tar,
        read_qc_results_tar = read_qc.results_tar,
        coverage_results_tar = coverage.results_tar,
        sample_quality_control_results_tar = sample_quality_control.results_tar,
        population_structure_results_tar = population_structure.results_tar,
        sample_map_results_tar = sample_map.results_tar,
        allele_frequencies_results_tar = allele_frequencies.results_tar,
        snp_dataframe_results_tar = snp_dataframe.results_tar,
        genetic_diversity_results_tar = genetic_diversity.results_tar,
        species_id_results_tar = species_id.results_tar,
        kdr_analysis_results_tar = kdr_analysis.results_tar,
        sample_ids = sample_ids, fastp_htmls = fastp_html,
        docker = ampseeker_docker,
    }
  }

  output {
    File  metadata_tsv          = prepare_metadata.metadata_tsv
    File  metadata_qcpass_tsv   = final_metadata
    File  metadata_colours_json = final_colours
    Array[File] sample_bams     = sample_bam
    Array[File] sample_bais     = bam_index.bai
    File targets_annot_vcf      = snp_eff_targets.annot_vcf
    File amplicons_annot_vcf    = snp_eff_amplicons.annot_vcf
    File? targets_vcf_stats     = vcf_stats_targets.stats
    File? amplicons_vcf_stats   = vcf_stats_amplicons.stats
    File? multiqc_report        = multi_qc.report_html
    File? multiqc_report_tar    = multi_qc.report_tar
    File? coverage_excel        = coverage.excel
    File sample_quality_control_nb = sample_quality_control.nb
    File? population_structure_nb  = population_structure.nb
    File? sample_map_nb            = sample_map.nb
    File? allele_frequencies_nb    = allele_frequencies.nb
    File? genetic_diversity_nb     = genetic_diversity.nb
    File  snp_dataframe_nb         = snp_dataframe.nb
    File  targets_snps_xlsx        = snp_dataframe.targets_xlsx
    File  amplicons_snps_xlsx      = snp_dataframe.amplicons_xlsx
    File? species_id_nb          = species_id.nb
    File? species_aims_tsv       = species_id.aims_tsv
    File? kdr_analysis_nb        = kdr_analysis.nb
    File? kdr_origins_tsv        = kdr_analysis.kdr_origins_tsv
    File? kdr_origin_counts_xlsx = kdr_analysis.kdr_origin_counts_xlsx
    File? kdr_origin_freqs_xlsx  = kdr_analysis.kdr_origin_freqs_xlsx
    File? jupyterbook_index      = build_report.home_page
    File? jupyterbook_site       = build_report.site_tar
    File? jupyterbook_standalone = build_report.standalone_html
  }
}

# ============================================================================
#                                  TASKS
# ============================================================================

# ---------------------------------------------------------------------------
# rules/metadata.smk
# ---------------------------------------------------------------------------
task prepare_metadata {
  input {
    Int disk_size_gb = 20
    Boolean from_bcl
    Array[File] sample_sheets
    File? metadata_file
    Array[String] cohort_columns
    String platform
    String docker
  }

  command <<<
    set -euo pipefail

    SHEETS_FILE="~{write_lines(sample_sheets)}"
    COHORT_COLS_FILE="~{write_json(cohort_columns)}"

    python3 <<PYEOF
import pandas as pd, re, os, json, numpy as np

with open("$COHORT_COLS_FILE") as f:
    cohort_cols = json.load(f)

from_bcl = "~{from_bcl}" == "true"
platform = "~{platform}"

def metadata_from_sample_sheet(p):
    with open(p) as f: content = f.read()
    m = re.search(r'\[(BCLConvert_)?Data\],+\n(.*?)($|\[)', content, re.DOTALL)
    if not m: raise ValueError("Could not find [Data] or [BCLConvert_Data] section in the sample sheet")
    lines = m.group(2).strip().split('\n')
    df = pd.DataFrame([l.split(',') for l in lines[1:]], columns=lines[0].split(','))
    df = df.rename(columns={'sample_ID': 'sample_id'})
    well_col = next((c for c in df.columns if c.lower() in ('well', 'sample_well')), None)
    if well_col:
        df['well_letter'] = df[well_col].str[0]
        df['well_number'] = df[well_col].str[1:].astype(str)
        df = df.drop(columns=[well_col])
    if 'plate_name' in df.columns: df = df.rename(columns={'plate_name': 'plate'})
    plate_col = next((c for c in df.columns if c.lower() in ('plate', 'sample_plate')), None)
    if plate_col and plate_col != 'plate':
        df = df.rename(columns={plate_col: 'plate'})
    return df

def metadata_from_sample_sheets(paths):
    metadata = []
    for p in paths:
        df = metadata_from_sample_sheet(p).copy()
        df["illumina_run_dir"] = os.path.dirname(p)
        metadata.append(df)
    metadata = pd.concat(metadata, ignore_index=True)
    dup = metadata["sample_id"].duplicated(keep=False)
    if dup.any():
        cols = [c for c in metadata.columns if c != "illumina_run_dir"]
        conflicts = [sid for sid, g in metadata.loc[dup, cols].groupby("sample_id")
                     if any(g[c].dropna().astype(str).nunique() > 1 for c in g.columns if c != "sample_id")]
        if conflicts: raise ValueError("Duplicate sample_id values with conflicting metadata: " + ", ".join(map(str, conflicts)))
        metadata = metadata.drop_duplicates(subset=["sample_id"], keep="first")
    return metadata

sheet_paths = [l.strip() for l in open("$SHEETS_FILE") if l.strip()] if os.path.exists("$SHEETS_FILE") else []
if from_bcl:
    metadata = metadata_from_sample_sheets(sheet_paths)
    # If a metadata file is also provided, merge its extra columns in by sample_id
    p = "~{if defined(metadata_file) then metadata_file else ''}"
    if p:
        if p.endswith(".xlsx"): extra = pd.read_excel(p, engine="openpyxl")
        elif p.endswith(".tsv"): extra = pd.read_csv(p, sep="\t")
        elif p.endswith(".csv"): extra = pd.read_csv(p, sep=",")
        else: extra = None
        if extra is not None:
            if 'sample_ID' in extra.columns: extra = extra.rename(columns={'sample_ID': 'sample_id'})
            # Only merge columns not already in metadata (avoid duplicates)
            new_cols = [c for c in extra.columns if c != 'sample_id' and c not in metadata.columns]
            if new_cols:
                metadata = metadata.merge(extra[['sample_id'] + new_cols], on='sample_id', how='left')
else:
    p = "~{if defined(metadata_file) then metadata_file else ''}"
    if p.endswith(".xlsx"): metadata = pd.read_excel(p, engine="openpyxl")
    elif p.endswith(".tsv"): metadata = pd.read_csv(p, sep="\t")
    elif p.endswith(".csv"): metadata = pd.read_csv(p, sep=",")
    elif not p and len(sheet_paths) > 0:
        metadata = metadata_from_sample_sheets(sheet_paths)
    elif not p:
        raise ValueError("Either metadata_file or illumina_sample_sheets must be provided")
    else: raise ValueError("Metadata file must be .xlsx, .tsv or .csv")

if not cohort_cols:
    cohort_cols = ["sample_id"]
if cohort_cols:
    assert all(c in metadata.columns for c in cohort_cols), f"Not all provided cohort columns are present in the metadata file. Available columns: {list(metadata.columns)}"
if 'sample_ID' in metadata.columns: metadata = metadata.rename(columns={'sample_ID': 'sample_id'})

os.makedirs("results/config", exist_ok=True)
metadata.to_csv("results/config/metadata.tsv", sep="\t", index=False)
with open("sample_ids.txt", "w") as f:
    for sid in metadata["sample_id"]: f.write(f"{sid}\n")
plate_info = np.isin(["plate", "well_letter", "well_number"], metadata.columns).all()
open("plate_info.txt", "w").write("true" if plate_info else "false")
if platform == "nanopore":
    assert 'fq1' in metadata.columns, "fq1 column required in metadata for nanopore data"
PYEOF
    cp results/config/metadata.tsv metadata.tsv
  >>>
  output {
    File metadata_tsv = "metadata.tsv"
    File sample_ids_txt = "sample_ids.txt"
    File plate_info_txt = "plate_info.txt"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/bcl-convert.smk (bcl2fastq needs its own licensed image)
# ---------------------------------------------------------------------------
task bcl_convert_run {
  input {
    Int disk_size_gb = 200
    File bcl_tar
    File sample_sheet
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p run_in run_out
    tar -xf ~{bcl_tar} -C run_in
    RUNDIR=$(find run_in -mindepth 1 -maxdepth 1 -type d | head -n1); RUNDIR=${RUNDIR:-run_in}
    bcl2fastq --runfolder-dir "$RUNDIR" --output-dir run_out --sample-sheet ~{sample_sheet} \
      --ignore-missing-bcls --ignore-missing-filter --ignore-missing-positions --ignore-missing-controls
    tar -cf run_out.tar -C run_out .
  >>>
  output {
    File run_output_tar = "run_out.tar"
    File demultiplex_stats = "run_out/Stats/DemultiplexingStats.xml"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 8
    memory: "16 GB"
  }
}

task bcl_convert_gather {
  input {
    Array[File] run_output_tars
    Array[File] demultiplex_stats
    Array[String] sample_ids
    String docker
    Int disk_size_gb = 500
  }
  command <<<
    set -euo pipefail
    mkdir -p all_runs out
    i=0
    for t in ~{sep=" " run_output_tars}; do i=$((i+1)); mkdir -p all_runs/run_$i; tar -xf "$t" -C all_runs/run_$i; done
    SAMPLE_IDS_FILE="~{write_lines(sample_ids)}"
    python3 <<PYEOF
import glob, json, shutil, os
samples = [l.strip() for l in open("$SAMPLE_IDS_FILE") if l.strip()]
os.makedirs("out/r1", exist_ok=True)
os.makedirs("out/r2", exist_ok=True)
sorted_samples = sorted(set(samples))

for idx, sample in enumerate(sorted_samples):
    # Find R1 and R2 — search recursively under all_runs
    r1 = sorted(glob.glob(f"all_runs/**/{sample}_*_R1_001.fastq.gz", recursive=True))
    r2 = sorted(glob.glob(f"all_runs/**/{sample}_*_R2_001.fastq.gz", recursive=True))
    if not r1: raise ValueError(f"No R1 FASTQ found for sample: {sample}. Available: {glob.glob('all_runs/**/*.fastq.gz', recursive=True)[:5]}")
    if not r2: raise ValueError(f"No R2 FASTQ found for sample: {sample}")
    if len(r1) != len(r2): raise ValueError(f"R1/R2 count mismatch for {sample}: {len(r1)} R1 vs {len(r2)} R2")
    # Use zero-padded index to guarantee glob returns in same order for r1 and r2
    tag = f"{idx:04d}_{sample}"
    if len(r1) == 1:
        # Single run — symlink directly, no concatenation
        os.symlink(os.path.abspath(r1[0]), f"out/r1/{tag}_R1.fastq.gz")
        os.symlink(os.path.abspath(r2[0]), f"out/r2/{tag}_R2.fastq.gz")
    else:
        # Multiple runs — concatenate keeping R1 with R1, R2 with R2
        with open(f"out/r1/{tag}_R1.fastq.gz", "wb") as o:
            for f in r1:
                with open(f, "rb") as fh: shutil.copyfileobj(fh, o)
        with open(f"out/r2/{tag}_R2.fastq.gz", "wb") as o:
            for f in r2:
                with open(f, "rb") as fh: shutil.copyfileobj(fh, o)

open("out/sample_order.txt", "w").write("\n".join(sorted_samples) + "\n")

i1 = sorted(glob.glob("all_runs/**/*_I1_001.fastq.gz", recursive=True))
i2 = sorted(glob.glob("all_runs/**/*_I2_001.fastq.gz", recursive=True))
if i1:
    with open("out/I1_combined.fastq.gz", "wb") as o:
        for f in i1:
            with open(f, "rb") as fh: shutil.copyfileobj(fh, o)
if i2:
    with open("out/I2_combined.fastq.gz", "wb") as o:
        for f in i2:
            with open(f, "rb") as fh: shutil.copyfileobj(fh, o)

merged = {"ConversionResults": []}
for sp in glob.glob("all_runs/**/Stats.json", recursive=True):
    merged["ConversionResults"].extend(json.load(open(sp)).get("ConversionResults", []))
json.dump(merged, open("out/Stats.json", "w"))
PYEOF
    cp ~{demultiplex_stats[0]} out/DemultiplexingStats.xml
  >>>
  output {
    Array[File] fastq_r1_sorted = glob("out/r1/*_R1.fastq.gz")
    Array[File] fastq_r2_sorted = glob("out/r2/*_R2.fastq.gz")
    Array[String] sample_order  = read_lines("out/sample_order.txt")
    File index1_reads = "out/I1_combined.fastq.gz"
    File index2_reads = "out/I2_combined.fastq.gz"
    File stats_json = "out/Stats.json"
    File demultiplex_stats = "out/DemultiplexingStats.xml"
  }
  runtime {
    docker: docker
    cpu: 2
    memory: "4 GB"
    disks: "local-disk ~{disk_size_gb} HDD"
  }
}

# ---------------------------------------------------------------------------
# sample_id -> FASTQ lookup, done in Python rather than WDL Map/zip/as_map
# (Terra's Cromwell doesn't implement the WDL-1.1 collection functions).
# Given the canonical sample order (from metadata) and whatever
# (id, file) pairs are actually available (from BCL demux or user input),
# emits one File per sample_id, in canonical order. Errors loudly if a
# sample_id has no matching file rather than silently misaligning.
#
# NOTE (fixed): this task used to write raw absolute filesystem paths into
# a text file and declare `Array[File] ordered = read_lines("ordered.txt")`.
# That pattern is unreliable for delocalization on GCP Batch when the
# referenced paths sit outside the task's own working directory — one
# shard's fastq silently failed to make it back to GCS, causing a
# downstream "URLs matched no objects" failure in fastp_illumina. Fixed by
# physically copying each matched file into a local `ordered/` subdirectory
# and outputting via glob(), which Cromwell always delocalizes correctly.
# ---------------------------------------------------------------------------
task resolve_fastq_array {
  input {
    Int disk_size_gb = 50
    Array[String] sample_ids
    Array[String] avail_sample_ids
    Array[File] avail_files
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p ordered
    python3 <<PYEOF
import shutil, os
sample_ids = open("~{write_lines(sample_ids)}").read().splitlines()
avail_ids  = open("~{write_lines(avail_sample_ids)}").read().splitlines()
avail_files = open("~{write_lines(avail_files)}").read().splitlines()
if len(avail_ids) != len(avail_files):
    raise SystemExit(f"avail_sample_ids ({len(avail_ids)}) and avail_files ({len(avail_files)}) length mismatch")
lookup = dict(zip(avail_ids, avail_files))
missing = [s for s in sample_ids if s not in lookup]
if missing:
    raise SystemExit("No FASTQ provided for sample_id(s): " + ", ".join(missing))
for i, s in enumerate(sample_ids):
    src = lookup[s]
    ext = ".fastq.gz" if src.endswith(".gz") else os.path.splitext(src)[1]
    dst = f"ordered/{i:04d}{ext}"
    shutil.copy(src, dst)
PYEOF
  >>>
  output {
    Array[File] ordered = glob("ordered/*")
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "1 GB"
  }
}

task index_read_fastqc {
  input {
    Int disk_size_gb = 20
    File i1
    File i2
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p results/qc/index-read-qc
    zcat ~{i1} | fastqc stdin --outdir results/qc/index-read-qc/
    mv results/qc/index-read-qc/stdin_fastqc.html results/qc/index-read-qc/I1.html
    mv results/qc/index-read-qc/stdin_fastqc.zip results/qc/index-read-qc/I1.zip
    zcat ~{i2} | fastqc stdin --outdir results/qc/index-read-qc/
    mv results/qc/index-read-qc/stdin_fastqc.html results/qc/index-read-qc/I2.html
    mv results/qc/index-read-qc/stdin_fastqc.zip results/qc/index-read-qc/I2.zip
  >>>
  output {
    File index1_qc = "results/qc/index-read-qc/I1.html"
    File index2_qc = "results/qc/index-read-qc/I2.html"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/utilities.smk
# ---------------------------------------------------------------------------
task reference_index {
  input {
    Int disk_size_gb = 20
    File reference_fasta
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{reference_fasta} reference.fasta
    samtools faidx reference.fasta
  >>>
  output {
    File fai = "reference.fasta.fai"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

task bwa_index {
  input {
    Int disk_size_gb = 20
    File reference_fasta
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{reference_fasta} reference.fasta
    bwa index reference.fasta
    tar -cf bwa_index.tar reference.fasta.amb reference.fasta.ann reference.fasta.bwt reference.fasta.pac reference.fasta.sa
  >>>
  output {
    File index_tar = "bwa_index.tar"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task minimap2_index {
  input {
    Int disk_size_gb = 20
    File reference_fasta
    String docker
  }
  command <<< set -euo pipefail; minimap2 -t 4 -d reference.mmi ~{reference_fasta} >>>
  output {
    File mmi = "reference.mmi"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "4 GB"
  }
}

task bam_index {
  input {
    Int disk_size_gb = 50
    String sample_id
    File bam
    String docker
  }
  command <<< set -euo pipefail; samtools index ~{bam} ~{sample_id}.bam.bai >>>
  output {
    File bai = "~{sample_id}.bam.bai"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

task bgzip_tabix {
  input {
    Int disk_size_gb = 50
    File vcf
    String docker
  }
  String base = basename(vcf)
  command <<< set -euo pipefail; bgzip -c ~{vcf} > ~{base}.gz; tabix ~{base}.gz >>>
  output {
    File vcfgz = "~{base}.gz"
    File tbi = "~{base}.gz.tbi"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

task create_ploidy_file {
  input {
    Int disk_size_gb = 20
    Int ploidy
    String docker
  }
  command <<<
    set -euo pipefail
    printf '* * * M %d\n' ~{ploidy} > ploidy.txt
    printf '* * * F %d\n' ~{ploidy} >> ploidy.txt
  >>>
  output {
    File ploidy_file = "ploidy.txt"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "1 GB"
  }
}

task merge_vcfs {
  input {
    Int disk_size_gb = 50
    String call_type
    String dataset
    Array[File] vcfs
    Array[File] tbis
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p in
    i=0; for f in ~{sep=" " vcfs}; do i=$((i+1)); ln -s "$f" in/$i.vcf.gz; done
    j=0; for f in ~{sep=" " tbis}; do j=$((j+1)); ln -s "$f" in/$j.vcf.gz.tbi; done
    n=$(ls in/*.vcf.gz | wc -l)
    if [ "$n" -gt 1000 ]; then
      half=$((n/2))
      ls in/*.vcf.gz | head -n "$half" > half1.txt
      ls in/*.vcf.gz | tail -n +"$((half+1))" > half2.txt
      bcftools merge --threads 12 -o merged.1.vcf -Ov -l half1.txt
      bcftools merge --threads 12 -o merged.2.vcf -Ov -l half2.txt
      bgzip -c merged.1.vcf > merged.1.vcf.gz && tabix merged.1.vcf.gz
      bgzip -c merged.2.vcf > merged.2.vcf.gz && tabix merged.2.vcf.gz
      bcftools merge -o ~{dataset}.merged.vcf -Ov merged.1.vcf.gz merged.2.vcf.gz
    else
      bcftools merge --threads 12 -o ~{dataset}.merged.vcf -Ov in/*.vcf.gz --force-samples
    fi
  >>>
  output {
    File merged_vcf = "~{dataset}.merged.vcf"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 12
    memory: "16 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/map-call-illumina.smk
# ---------------------------------------------------------------------------
task fastp_illumina {
  input {
    Int disk_size_gb = 50
    String sample_id
    File fq1
    File fq2
    String docker
  }
  command <<<
    set -euo pipefail
    fastp -i ~{fq1} -I ~{fq2} -o ~{sample_id}_1.fastq.gz -O ~{sample_id}_2.fastq.gz --html ~{sample_id}.html --json ~{sample_id}.json
  >>>
  output {
    File trimmed_r1 = "~{sample_id}_1.fastq.gz"
    File trimmed_r2 = "~{sample_id}_2.fastq.gz"
    File html = "~{sample_id}.html"
    File json = "~{sample_id}.json"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "4 GB"
  }
}

task bwa_align {
  input {
    Int disk_size_gb = 100
    String sample_id
    File fq1
    File fq2
    File reference_fasta
    File bwa_index_tar
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{reference_fasta} reference.fasta
    tar -xf ~{bwa_index_tar}
    bwa mem -t 4 -R '@RG\tID:~{sample_id}\tSM:~{sample_id}\tPU:nye\tPL:nye\tLB:~{sample_id}_lb~{sample_id}' \
      reference.fasta ~{fq1} ~{fq2} | samtools sort -@4 -o ~{sample_id}.bam
  >>>
  output {
    File bam = "~{sample_id}.bam"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "8 GB"
  }
}

task mpileup_call_amplicons_illumina {
  input {
    Int disk_size_gb = 50
    String sample_id
    File bam
    File bam_index
    File reference_fasta
    File ploidy_file
    String docker
  }
  command <<<
    set -euo pipefail
    bcftools mpileup -Ov -I -f ~{reference_fasta} -a AD --max-depth 2000 ~{bam} |
      bcftools call -m --ploidy-file ~{ploidy_file} -Ov | bcftools sort -Ov -o ~{sample_id}.calls.vcf
  >>>
  output {
    File calls_vcf = "~{sample_id}.calls.vcf"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "4 GB"
  }
}

task mpileup_call_targets {
  input {
    Int disk_size_gb = 50
    String sample_id
    File bam
    File bam_index
    File reference_fasta
    File reference_fai
    File ploidy_file
    File targets_bed
    String platform
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{bam} ~{sample_id}.bam
    ln -s ~{bam_index} ~{sample_id}.bam.bai
    ln -s ~{reference_fasta} reference.fasta
    ln -s ~{reference_fai} reference.fasta.fai
    ONT_FLAG=""
    if [ "~{platform}" = "nanopore" ]; then ONT_FLAG="-X ont-sup"; fi
    bcftools mpileup $ONT_FLAG -Ov -f reference.fasta -R ~{targets_bed} -a AD --max-depth 2000 ~{sample_id}.bam |
      bcftools call -f GQ,GP -m --ploidy-file ~{ploidy_file} -Ov | bcftools sort -Ov -o ~{sample_id}.calls.vcf
  >>>
  output {
    File calls_vcf = "~{sample_id}.calls.vcf"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "4 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/map-call-nanopore.smk (fastplong needs adding to the image, see DEPLOY.md)
# ---------------------------------------------------------------------------
task fastp_nanopore {
  input {
    Int disk_size_gb = 50
    String sample_id
    File fq1
    String docker
  }
  command <<<
    set -euo pipefail
    fastplong -i ~{fq1} -o ~{sample_id}.fastq.gz --html ~{sample_id}.html --json ~{sample_id}.json
  >>>
  output {
    File trimmed = "~{sample_id}.fastq.gz"
    File html = "~{sample_id}.html"
    File json = "~{sample_id}.json"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "4 GB"
  }
}

task minimap2_align {
  input {
    Int disk_size_gb = 100
    String sample_id
    File reads
    File reference_fasta
    File reference_mmi
    String docker
  }
  command <<<
    set -euo pipefail
    minimap2 -t 8 -ax map-ont -R '@RG\tID:~{sample_id}\tSM:~{sample_id}\tPL:ONT\tLB:~{sample_id}_lib' \
      ~{reference_mmi} ~{reads} | samtools sort -@8 -o ~{sample_id}.bam
  >>>
  output {
    File bam = "~{sample_id}.bam"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 8
    memory: "8 GB"
  }
}

task clair3_call_amplicons_nanopore {
  input {
    Int disk_size_gb = 100
    String sample_id
    File bam
    File bam_index
    File reference_fasta
    File reference_fai
    File? model_tar
    Int ploidy
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{bam} ~{sample_id}.bam
    ln -s ~{bam_index} ~{sample_id}.bam.bai
    ln -s ~{reference_fasta} reference.fasta
    ln -s ~{reference_fai} reference.fasta.fai
    mkdir -p outdir
    if [ -n "~{if defined(model_tar) then model_tar else ''}" ]; then
      mkdir -p model
      tar -xf ~{if defined(model_tar) then model_tar else ''} -C model
      MODEL_PATH=$(find model -name "pileup.pt" | head -1 | xargs dirname 2>/dev/null || find model -maxdepth 2 -mindepth 1 -type d | head -1 || echo "model")
    else
      # Use model bundled in the clair3 image at /opt/models/
      MODEL_PATH="/opt/models/r1041_e82_400bps_sup_v500"
    fi
    echo "Using model path: $MODEL_PATH"
    PLOIDY_FLAGS=""
    if [ "~{ploidy}" = "1" ]; then PLOIDY_FLAGS="--no_phasing_for_fa --haploid_sensitive"; fi
    run_clair3.sh --bam_fn=~{sample_id}.bam --ref_fn=reference.fasta --threads=1 --platform=ont \
      --model_path=$MODEL_PATH --output=outdir --include_all_ctgs $PLOIDY_FLAGS --sample_name=~{sample_id}
    gunzip -c outdir/merge_output.vcf.gz > ~{sample_id}.calls.vcf
  >>>
  output {
    File calls_vcf = "~{sample_id}.calls.vcf"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "16 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/qc.smk
# ---------------------------------------------------------------------------
task mosdepth_coverage {
  input {
    Int disk_size_gb = 50
    String sample_id
    File bam
    File bam_index
    File panel_bed
    String docker
  }
  command <<<
    set -euo pipefail
    ln -s ~{bam} ~{sample_id}.bam
    ln -s ~{bam_index} ~{sample_id}.bam.bai
    mosdepth ~{sample_id} ~{sample_id}.bam --fast-mode --by ~{panel_bed} --threads 4
  >>>
  output {
    File per_base_bed = "~{sample_id}.per-base.bed.gz"
    File regions_bed = "~{sample_id}.regions.bed.gz"
    File summary = "~{sample_id}.mosdepth.summary.txt"
    File global_dist = "~{sample_id}.mosdepth.global.dist.txt"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "4 GB"
  }
}

task bam_stats {
  input {
    Int disk_size_gb = 50
    String sample_id
    File bam
    File bam_index
    String docker
  }
  command <<< set -euo pipefail; samtools flagstat ~{bam} > ~{sample_id}.flagstat >>>
  output {
    File flagstat = "~{sample_id}.flagstat"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

task vcf_stats {
  input {
    Int disk_size_gb = 50
    String call_type
    String dataset
    File merged_vcf
    String docker
  }
  command <<< set -euo pipefail; bcftools stats ~{merged_vcf} > ~{dataset}-~{call_type}.merged.vcf.txt >>>
  output {
    File stats = "~{dataset}-~{call_type}.merged.vcf.txt"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "4 GB"
  }
}

task multi_qc {
  input {
    Int disk_size_gb = 100
    Array[File] fastp_jsons
    Array[File] flagstats
    Array[File] vcf_stats_files
    Array[File] per_base_beds
    Array[File] summary_files
    Array[File] dist_files
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p qc_in
    cp ~{sep=" " fastp_jsons} qc_in/ 2>/dev/null || true
    cp ~{sep=" " flagstats} qc_in/ 2>/dev/null || true
    cp ~{sep=" " vcf_stats_files} qc_in/ 2>/dev/null || true
    # Note: per_base_beds and global_dist files are not passed to MultiQC for
    # amplicon data — global coverage appears as ~0 since reads only align to
    # tiny target regions, which is misleading. Summary files still show
    # per-region depth correctly.
    cp ~{sep=" " summary_files} qc_in/ 2>/dev/null || true
    mkdir -p results/qc/multiqc
    multiqc qc_in -n multiqc_report -o results/qc/multiqc

    # multiqc_report.html can reference sibling assets (plot data / flat-mode
    # plot images) depending on sample count and config. Delocalizing only
    # the HTML file drops those, same issue as jupyter-book's index.html —
    # archive the whole output directory too.
    tar -czf multiqc_report_full.tar.gz -C results/qc/multiqc .
  >>>
  output {
    File report_html = "results/qc/multiqc/multiqc_report.html"
    File report_tar  = "multiqc_report_full.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "16 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/snpeff.smk
# ---------------------------------------------------------------------------
task snp_eff_db_download {
  input {
    Int disk_size_gb = 50
    String reference_snpeffdb
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p mysnpeffdb
    snpEff download ~{reference_snpeffdb} -dataDir mysnpeffdb
    touch mysnpeffdb/.db.dl
    tar -cf db.tar -C mysnpeffdb .
  >>>
  output {
    File done_marker = "db.tar"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "8 GB"
  }
}

task create_custom_snp_eff_db {
  input {
    Int disk_size_gb = 50
    File reference_fasta
    File reference_gff3
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p mysnpeffdb
    ln -s ~{reference_fasta} mysnpeffdb/sequences.fa
    ln -s ~{reference_gff3} mysnpeffdb/genes.gff
    cp /opt/snpEff/snpEff/snpEff.config snpEff.config
    echo 'mysnpeffdb.genome : mysnpeffdb' >> snpEff.config
    snpEff build -gff3 -v -c snpEff.config -dataDir "$(pwd)" mysnpeffdb -noCheckCds -noCheckProtein
    tar -cf db.tar mysnpeffdb snpEff.config
  >>>
  output {
    File db_tar = "db.tar"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}

task snp_eff {
  input {
    Int disk_size_gb = 50
    String call_type
    String dataset
    File merged_vcf
    Boolean custom_snpeffdb
    String reference_snpeffdb
    File db_ready
    String docker
  }
  command <<<
    set -euo pipefail

    if [ "~{custom_snpeffdb}" = "true" ]; then
      tar -xf ~{db_ready}

      test -f mysnpeffdb/snpEffectPredictor.bin

      snpEff eff \
        -c snpEff.config \
        -dataDir "$(pwd)" \
        -stats ~{dataset}.summary.html \
        -csvStats ~{dataset}.summary.csv \
        -ud 0 \
        mysnpeffdb \
        ~{merged_vcf} \
        > ~{dataset}.annot.vcf
    else
      mkdir -p dataDir
      tar -xf ~{db_ready} -C dataDir

      snpEff eff         -dataDir dataDir \
        -stats ~{dataset}.summary.html \
        -csvStats ~{dataset}.summary.csv \
        -ud 0 \
        ~{reference_snpeffdb} \
        ~{merged_vcf} \
        > ~{dataset}.annot.vcf
    fi
  >>>
  output {
    File annot_vcf = "~{dataset}.annot.vcf"
    File stats = "~{dataset}.summary.html"
    File csv_stats = "~{dataset}.summary.csv"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "8 GB"
  }
}

# ---------------------------------------------------------------------------
# Snakefile top-level (create_color_mapping)
# ---------------------------------------------------------------------------
task metadata_colours {
  input {
    Int disk_size_gb = 20
    File metadata_tsv
    Array[String] cohort_columns
    File? user_colours_json
    String docker
  }
  command <<<
    set -euo pipefail
    if [ -n "~{if defined(user_colours_json) then user_colours_json else ''}" ]; then
      cp ~{user_colours_json} metadata_colours.json
    else
      COHORT_COLS_FILE="~{write_json(cohort_columns)}"
      python3 <<PYEOF
import pandas as pd, json
PALETTES = {
    'Plotly': ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52'],
    'Safe': ['#88CCEE','#CC6677','#DDCC77','#117733','#332288','#AA4499','#44AA99','#999933','#882255','#661100','#6699CC','#888888'],
    'Vivid': ['#E58606','#5D69B1','#52BCA3','#99C945','#CC61B0','#24796C','#DAA51B','#2F8AC4','#764E9F','#ED645A','#CC3A8E','#A5AA99'],
    'Pastel': ['#B6B6F2','#FFB6B6','#BFFFBF','#FFE2BF','#BFE2FF','#FFB6FF','#E2FFB6','#FFD9B6','#B6FFE2','#D9B6FF'],
    'Dark2': ['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666'],
}
def create_color_mapping(df, columns, repeated_colors=True, unassigned_value='unassigned'):
    m = {}
    for i, col in enumerate(columns):
        colors = PALETTES[list(PALETTES.keys())[i % len(PALETTES)]]
        cc = {}
        for j, v in enumerate(df[col].unique()):
            cc[v] = "#808080" if v == unassigned_value else colors[j % len(colors) if repeated_colors else j]
        m[col] = cc
    return m
df = pd.read_csv("~{metadata_tsv}", sep="\t")
cols = json.load(open("$COHORT_COLS_FILE"))
json.dump(create_color_mapping(df, columns=cols), open("metadata_colours.json", "w"), indent=2)
PYEOF
    fi
  >>>
  output {
    File colours_json = "metadata_colours.json"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

# ---------------------------------------------------------------------------
# Notebook tasks. All notebooks, workflow/lib, docs/ampseeker-results, and
# resources/ag-vampir/* live baked into the image at /opt/AmpSeeker (from
# the Dockerfile's `COPY . /opt/AmpSeeker`); tasks reference those baked-in
# paths directly. `wkdir="."` is the task's own ephemeral working directory,
# which is why each task still stages the specific results/*, config/*
# files a given notebook hardcodes reads to — that part hasn't changed,
# only where the notebook/library code itself comes from.
# ---------------------------------------------------------------------------
task run_information {
  input {
    Int disk_size_gb = 20
    File metadata
    File targets_bed
    String panel
    String dataset
    String platform
    Array[String] cohort_columns
    Array[String] contigs
    File? config_yaml
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow && ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    CONFIG_PATH="~{if defined(config_yaml) then config_yaml else ''}"
    if [ -z "$CONFIG_PATH" ]; then CONFIG_PATH="/opt/AmpSeeker/config/config.yaml"; fi

    papermill /opt/AmpSeeker/workflow/notebooks/run-information.ipynb run-information.ipynb -k AmpSeq_python \
      -p config_path "$CONFIG_PATH" -p metadata_path ~{metadata} -p bed_targets_path ~{targets_bed} \
      -p panel "~{panel}" -p dataset "~{dataset}" -p platform "~{platform}" \
      -p cohort_cols "~{sep="," cohort_columns}" -p contigs "~{sep="," contigs}" -p wkdir "."

    # Collect whatever this notebook wrote under results/ (plots, tables)
    # so build_report can reconstruct a shared results/ tree and resolve
    # any relative link/embed the rendered page makes to them.
    mkdir -p results
    tar -czf run_information_results.tar.gz results
  >>>
  output {
    File nb_out = "run-information.ipynb"
    File docs_nb = "run-information.ipynb"
    File results_tar = "run_information_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task run_statistics {
  input {
    Int disk_size_gb = 50
    File metadata
    File stats_json
    File colours_json
    Array[String] cohort_columns
    Boolean plate_info
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/config resources/reads/Stats
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{colours_json} results/config/metadata_colours.json
    cp ~{stats_json} resources/reads/Stats/Stats.json

    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome
    papermill /opt/AmpSeeker/workflow/notebooks/run-statistics.ipynb run-statistics.ipynb -k AmpSeq_python \
      -p plate_info ~{plate_info} -p metadata_path ~{metadata} -p cohort_cols "~{sep="," cohort_columns}" -p wkdir "."

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf run_statistics_results.tar.gz results
  >>>
  output {
    File nb_out = "run-statistics.ipynb"
    File docs_nb = "run-statistics.ipynb"
    File results_tar = "run_statistics_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task reads_per_well {
  input {
    Int disk_size_gb = 20
    Array[String] sample_ids
    Array[File] bams
    Array[File] flagstats
    File metadata
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/alignments/bamStats
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    SIDS="~{write_lines(sample_ids)}"; BAMS="~{write_lines(bams)}"; FLAGSTATS="~{write_lines(flagstats)}"
    paste "$SIDS" "$BAMS" | while IFS=$'\t' read -r sid bam; do cp "$bam" "results/alignments/${sid}.bam"; done
    paste "$SIDS" "$FLAGSTATS" | while IFS=$'\t' read -r sid fs; do cp "$fs" "results/alignments/bamStats/${sid}.flagstat"; done

    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome
    papermill /opt/AmpSeeker/workflow/notebooks/reads-per-well.ipynb reads-per-well.ipynb -k AmpSeq_python \
      -p wkdir "." -p metadata_path ~{metadata}

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf reads_per_well_results.tar.gz results
  >>>
  output {
    File nb_out = "reads-per-well.ipynb"
    File docs_nb = "reads-per-well.ipynb"
    File results_tar = "reads_per_well_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task read_qc {
  input {
    Int disk_size_gb = 50
    File metadata
    Boolean index_read_qc
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow && ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    papermill /opt/AmpSeeker/workflow/notebooks/read-quality.ipynb read-quality.ipynb -k AmpSeq_python \
      -p metadata_path ~{metadata} -p index_read_qc ~{index_read_qc} -p wkdir "."

    mkdir -p results
    tar -czf read_qc_results.tar.gz results
  >>>
  output {
    File nb_out = "read-quality.ipynb"
    File docs_nb = "read-quality.ipynb"
    File results_tar = "read_qc_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task coverage {
  input {
    Int disk_size_gb = 50
    Array[String] sample_ids
    Array[File] per_base_beds
    Array[File] regions_beds
    File metadata
    File targets_bed
    Array[String] contigs
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/coverage
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    SIDS="~{write_lines(sample_ids)}"; PB="~{write_lines(per_base_beds)}"; RB="~{write_lines(regions_beds)}"
    paste "$SIDS" "$PB" | while IFS=$'\t' read -r sid f; do cp "$f" "results/coverage/${sid}.per-base.bed.gz"; done
    paste "$SIDS" "$RB" | while IFS=$'\t' read -r sid f; do cp "$f" "results/coverage/${sid}.regions.bed.gz"; done

    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome
    papermill /opt/AmpSeeker/workflow/notebooks/coverage.ipynb coverage.ipynb -k AmpSeq_python \
      -p metadata_path ~{metadata} -p bed_targets_path ~{targets_bed} -p contigs "~{sep="," contigs}" -p wkdir "."

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf coverage_results.tar.gz results
  >>>
  output {
    File nb_out = "coverage.ipynb"
    File docs_nb = "coverage.ipynb"
    File excel = "results/coverage/amplicon_by_sample_depth.xlsx"
    File results_tar = "coverage_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}

task sample_quality_control {
  input {
    Int disk_size_gb = 50
    Array[String] sample_ids
    Array[File] regions_beds
    File targets_vcf
    File amplicons_vcf
    File metadata
    File targets_bed
    Array[String] cohort_columns
    Float sample_min_mean_target_coverage
    String panel
    String dataset
    String platform
    File colours_json
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/coverage results/config results/vcfs/amplicons
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    SIDS="~{write_lines(sample_ids)}"; RB="~{write_lines(regions_beds)}"
    paste "$SIDS" "$RB" | while IFS=$'\t' read -r sid f; do cp "$f" "results/coverage/${sid}.regions.bed.gz"; done
    cp ~{colours_json} results/config/metadata_colours.json
    cp ~{amplicons_vcf} results/vcfs/amplicons/~{dataset}.annot.vcf

    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome
    papermill /opt/AmpSeeker/workflow/notebooks/sample-quality-control.ipynb sample-quality-control.ipynb -k AmpSeq_python \
      -p platform "~{platform}" -p dataset "~{dataset}" -p panel "~{panel}" -p metadata_path ~{metadata} \
      -p cohort_cols "~{sep="," cohort_columns}" -p bed_targets_path ~{targets_bed} -p vcf_path ~{targets_vcf} \
      -p wkdir "." -p sample_min_mean_target_coverage ~{sample_min_mean_target_coverage}

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf sample_quality_control_results.tar.gz results
  >>>
  output {
    File nb = "sample-quality-control.ipynb"
    File docs_nb = "sample-quality-control.ipynb"
    File qcpass_metadata = "results/config/metadata.qcpass.tsv"
    File results_tar = "sample_quality_control_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}

task population_structure {
  input {
    Int disk_size_gb = 50
    File amplicons_vcf
    File metadata_qcpass
    String dataset
    Array[String] cohort_columns
    String platform
    File colours_json
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/config
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{colours_json} results/config/metadata_colours.json

    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/population-structure.ipynb population-structure.ipynb -k AmpSeq_python \
      -p wkdir "." -p platform "~{platform}" -p metadata_path ~{metadata_qcpass} \
      -p dataset "~{dataset}" -p vcf_path ~{amplicons_vcf} -p cohort_cols "~{sep="," cohort_columns}"

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf population_structure_results.tar.gz results
  >>>
  output {
    File nb = "population-structure.ipynb"
    File docs_nb = "population-structure.ipynb"
    File results_tar = "population_structure_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 4
    memory: "8 GB"
  }
}

task sample_map {
  input {
    Int disk_size_gb = 20
    File metadata_qcpass
    String dataset
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow && ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/sample-map.ipynb sample-map.ipynb -k AmpSeq_python \
      -p metadata_path ~{metadata_qcpass} -p dataset "~{dataset}"

    mkdir -p results
    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf sample_map_results.tar.gz results
  >>>
  output {
    File nb = "sample-map.ipynb"
    File docs_nb = "sample-map.ipynb"
    File results_tar = "sample_map_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "4 GB"
  }
}

task allele_frequencies {
  input {
    Int disk_size_gb = 50
    File metadata_qcpass
    File targets_bed
    File targets_vcf
    File amplicons_vcf
    String dataset
    Array[String] cohort_columns
    String platform
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/vcfs/amplicons
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{amplicons_vcf} results/vcfs/amplicons/~{dataset}.annot.vcf

    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/allele-frequencies.ipynb allele-frequencies.ipynb -k AmpSeq_python \
      -p metadata_path ~{metadata_qcpass} -p platform "~{platform}" -p dataset "~{dataset}" -p wkdir "." \
      -p bed_path ~{targets_bed} -p cohort_cols "~{sep="," cohort_columns}" -p vcf_path ~{targets_vcf}

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf allele_frequencies_results.tar.gz results
  >>>
  output {
    File nb = "allele-frequencies.ipynb"
    File docs_nb = "allele-frequencies.ipynb"
    File results_tar = "allele_frequencies_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}

task snp_dataframe {
  input {
    Int disk_size_gb = 50
    File targets_vcf
    File amplicons_vcf
    String dataset
    String platform
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/vcfs/targets results/vcfs/amplicons
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{targets_vcf} results/vcfs/targets/~{dataset}.annot.vcf
    cp ~{amplicons_vcf} results/vcfs/amplicons/~{dataset}.annot.vcf

    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/snp-dataframe.ipynb snp-dataframe.ipynb -k AmpSeq_python \
      -p dataset "~{dataset}" -p platform "~{platform}" -p wkdir "."

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf snp_dataframe_results.tar.gz results
  >>>
  output {
    File nb = "snp-dataframe.ipynb"
    File docs_nb = "snp-dataframe.ipynb"
    File targets_xlsx = "results/vcfs/targets/~{dataset}-snps.xlsx"
    File amplicons_xlsx = "results/vcfs/amplicons/~{dataset}-snps.xlsx"
    File results_tar = "snp_dataframe_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task genetic_diversity {
  input {
    Int disk_size_gb = 50
    File amplicons_vcf
    File metadata_qcpass
    String dataset
    Array[String] cohort_columns
    String platform
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/genetic-diversity
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/genetic-diversity.ipynb genetic-diversity.ipynb -k AmpSeq_python \
      -p dataset "~{dataset}" -p platform "~{platform}" -p vcf_path ~{amplicons_vcf} \
      -p wkdir "." -p cohort_cols "~{sep="," cohort_columns}" -p metadata_path ~{metadata_qcpass}

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf genetic_diversity_results.tar.gz results
  >>>
  output {
    File nb = "genetic-diversity.ipynb"
    File docs_nb = "genetic-diversity.ipynb"
    Array[File] pi_tables = glob("results/genetic-diversity/*.pi.tsv")
    File samples_pi_table = "results/genetic-diversity/samples.pi.tsv"
    File results_tar = "genetic_diversity_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/ag-vampir.smk : species_id — resources/ag-vampir/* reference files
# are baked into the image; no longer WDL inputs.
# ---------------------------------------------------------------------------
task species_id {
  input {
    Int disk_size_gb = 50
    File targets_vcf
    File targets_bed
    File metadata_qcpass
    Array[String] cohort_columns
    File colours_json
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/config results/ag-vampir/aims

    # The species-id notebook expects the AG-VAMPIR reference resources
    # relative to its working directory:
    #   ./resources/ag-vampir/*.npy
    # These resources are baked into the AmpSeeker Docker image under
    # /opt/AmpSeeker/resources/ag-vampir.
    if [ ! -d /opt/AmpSeeker/resources/ag-vampir ]; then
      echo "ERROR: /opt/AmpSeeker/resources/ag-vampir is missing from the Docker image." >&2
      exit 1
    fi
    cp -a /opt/AmpSeeker/resources ./resources

    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{colours_json} results/config/metadata_colours.json
    cp ~{metadata_qcpass} results/config/metadata.qcpass.tsv

    # Plotly image export in species-id requires Kaleido and Chrome.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    # Filter cohort_cols to only columns that exist in metadata and are non-empty
    COHORT_COLS=$(python3 -c "
import pandas as pd, sys
cols = '~{sep="," cohort_columns}'.split(',')
meta = pd.read_csv('results/config/metadata.qcpass.tsv', sep='\t')
valid = [c.strip() for c in cols if c.strip() and c.strip() in meta.columns and c.strip() != 'sample_id']
print(','.join(valid))
")
    papermill /opt/AmpSeeker/workflow/notebooks/ag-vampir/species-id.ipynb species-id.ipynb -k AmpSeq_python \
      -p platform "illumina" -p metadata_path results/config/metadata.qcpass.tsv -p wkdir "." \
      -p vcf_path ~{targets_vcf} -p bed_targets_path ~{targets_bed} -p cohort_cols "$COHORT_COLS"

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf species_id_results.tar.gz results
  >>>
  output {
    File nb = "species-id.ipynb"
    File docs_nb = "species-id.ipynb"
    File aims_tsv = "results/ag-vampir/aims/taxon_aims.tsv"
    File updated_metadata = "results/config/metadata.qcpass.tsv"
    File results_tar = "species_id_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

task overwrite_metadata_colours {
  input {
    Int disk_size_gb = 20
    File metadata_qcpass
    Array[String] cohort_columns
    File? user_colours_json
    String docker
  }
  command <<<
    set -euo pipefail
    if [ -n "~{if defined(user_colours_json) then user_colours_json else ''}" ]; then
      cp ~{user_colours_json} metadata_colours.json
    else
      COHORT_COLS_FILE="~{write_json(cohort_columns)}"
      python3 <<PYEOF
import pandas as pd, json
PALETTES = {
    'Plotly': ['#636EFA','#EF553B','#00CC96','#AB63FA','#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52'],
    'Safe': ['#88CCEE','#CC6677','#DDCC77','#117733','#332288','#AA4499','#44AA99','#999933','#882255','#661100','#6699CC','#888888'],
    'Vivid': ['#E58606','#5D69B1','#52BCA3','#99C945','#CC61B0','#24796C','#DAA51B','#2F8AC4','#764E9F','#ED645A','#CC3A8E','#A5AA99'],
    'Pastel': ['#B6B6F2','#FFB6B6','#BFFFBF','#FFE2BF','#BFE2FF','#FFB6FF','#E2FFB6','#FFD9B6','#B6FFE2','#D9B6FF'],
    'Dark2': ['#1B9E77','#D95F02','#7570B3','#E7298A','#66A61E','#E6AB02','#A6761D','#666666'],
}
def create_color_mapping(df, columns, repeated_colors=True, unassigned_value='unassigned'):
    m = {}
    for i, col in enumerate(columns):
        colors = PALETTES[list(PALETTES.keys())[i % len(PALETTES)]]
        cc = {}
        for j, v in enumerate(df[col].unique()):
            cc[v] = "#808080" if v == unassigned_value else colors[j % len(colors) if repeated_colors else j]
        m[col] = cc
    return m
df = pd.read_csv("~{metadata_qcpass}", sep="\t")
cols = json.load(open("$COHORT_COLS_FILE"))
json.dump(create_color_mapping(df, columns=cols), open("metadata_colours.json", "w"), indent=2)
PYEOF
    fi
  >>>
  output {
    File colours_json = "metadata_colours.json"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 1
    memory: "2 GB"
  }
}

task kdr_analysis {
  input {
    Int disk_size_gb = 50
    File targets_vcf
    File amplicons_vcf
    File metadata_qcpass
    File? kdr_marker_snps_csv
    File targets_bed
    File snp_dataframe_targets_xlsx
    File colours_json
    String dataset
    Array[String] cohort_columns
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow results/vcfs/amplicons results/config config results/vcfs/targets results/ag-vampir/kdr-origins
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp ~{amplicons_vcf} results/vcfs/amplicons/~{dataset}.annot.vcf
    cp ~{colours_json} results/config/metadata_colours.json
    cp ~{targets_bed} config/ag-vampir.bed
    cp ~{snp_dataframe_targets_xlsx} results/vcfs/targets/~{dataset}-snps.xlsx

    KDR_CSV="~{if defined(kdr_marker_snps_csv) then kdr_marker_snps_csv else ''}"
    if [ -z "$KDR_CSV" ]; then KDR_CSV="/opt/AmpSeeker/resources/ag-vampir/Kdr_marker_SNPs.csv"; fi

    # Plotly image export dependency for this notebook.
    python3 -m pip install --no-cache-dir --upgrade kaleido
    kaleido_get_chrome

    papermill /opt/AmpSeeker/workflow/notebooks/ag-vampir/kdr-analysis.ipynb kdr-analysis.ipynb -k AmpSeq_python \
      -p platform "illumina" -p dataset "~{dataset}" -p metadata_path ~{metadata_qcpass} \
      -p vcf_path ~{targets_vcf} -p cohort_cols "~{sep="," cohort_columns}" -p wkdir "." \
      -p kdr_marker_snps_path "$KDR_CSV"

    python3 /opt/AmpSeeker/workflow/notebooks/save_pdfs.py
    tar -czf kdr_analysis_results.tar.gz results
  >>>
  output {
    File nb = "kdr-analysis.ipynb"
    File docs_nb = "kdr-analysis.ipynb"
    File kdr_origins_tsv = "results/ag-vampir/kdr-origins/kdr_origins.tsv"
    File kdr_genhap_origins_tsv = "results/ag-vampir/kdr-origins/kdr_genhap_origins.tsv"
    File kdr_origin_counts_xlsx = "results/ag-vampir/kdr-origins/kdr_origin_counts.xlsx"
    File kdr_origin_freqs_xlsx = "results/ag-vampir/kdr-origins/kdr_origin_freqs.xlsx"
    File results_tar = "kdr_analysis_results.tar.gz"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "4 GB"
  }
}

# ---------------------------------------------------------------------------
# rules/jupyter-book.smk, merged: process_toc + process_notebooks + build.
# docs/ampseeker-results and its _toc_skeleton.yml are baked into the image.
# ---------------------------------------------------------------------------
task build_report {
  input {
    Int disk_size_gb = 50
    Boolean from_bcl
    String panel
    Boolean plate_info
    String dataset
    File docs_run_info
    File? docs_run_statistics
    File docs_snp_df
    File docs_sample_quality_control
    File? docs_coverage
    File? docs_pop_structure
    File? docs_af
    File? docs_sample_map
    File? docs_genetic_diversity
    File docs_read_quality
    File? docs_reads_per_well
    File? docs_species_id
    File? docs_kdr_analysis
    # Each notebook task's full local results/ tree, tarred. Needed because
    # every notebook here writes plots/tables to results/... in its own
    # isolated task, and links/embeds them by relative path — but only its
    # .ipynb file was ever passed to build_report before, so those links
    # and images were dead. Extracting all of these into one shared
    # results/ tree here (mirroring the single shared working directory
    # the original Snakemake pipeline used) makes those references resolve.
    File run_information_results_tar
    File? run_statistics_results_tar
    File? reads_per_well_results_tar
    File read_qc_results_tar
    File? coverage_results_tar
    File sample_quality_control_results_tar
    File? population_structure_results_tar
    File? sample_map_results_tar
    File? allele_frequencies_results_tar
    File snp_dataframe_results_tar
    File? genetic_diversity_results_tar
    File? species_id_results_tar
    File? kdr_analysis_results_tar
    Array[String] sample_ids
    Array[File] fastp_htmls
    String docker
  }
  command <<<
    set -euo pipefail
    mkdir -p workflow docs
    ln -s /opt/AmpSeeker/workflow/lib workflow/lib
    cp -r /opt/AmpSeeker/docs/ampseeker-results docs/ampseeker-results
    mkdir -p docs/ampseeker-results/notebooks/ag-vampir

    # Re-stage each sample's fastp HTML report at the results/fastp/ path
    # read-quality.ipynb links to. Previously these were generated by
    # fastp_illumina/fastp_nanopore but never collected anywhere, so every
    # per-sample "fastp report" link on the Sample read QC page was dead.
    mkdir -p results/fastp
    SIDS="~{write_lines(sample_ids)}"; HTMLS="~{write_lines(fastp_htmls)}"
    paste "$SIDS" "$HTMLS" | while IFS=$'\t' read -r sid f; do cp "$f" "results/fastp/${sid}.html"; done

    # Reconstruct one shared results/ tree from every notebook task's own
    # tarred output, so relative links/embedded images in the rendered
    # pages below resolve. Later tars intentionally never overwrite an
    # earlier task's files, since each task only ever wrote into its own
    # results/<module>/ subdirectory.
    for t in "~{run_information_results_tar}" \
             "~{if defined(run_statistics_results_tar) then run_statistics_results_tar else ''}" \
             "~{if defined(reads_per_well_results_tar) then reads_per_well_results_tar else ''}" \
             "~{read_qc_results_tar}" \
             "~{if defined(coverage_results_tar) then coverage_results_tar else ''}" \
             "~{sample_quality_control_results_tar}" \
             "~{if defined(population_structure_results_tar) then population_structure_results_tar else ''}" \
             "~{if defined(sample_map_results_tar) then sample_map_results_tar else ''}" \
             "~{if defined(allele_frequencies_results_tar) then allele_frequencies_results_tar else ''}" \
             "~{snp_dataframe_results_tar}" \
             "~{if defined(genetic_diversity_results_tar) then genetic_diversity_results_tar else ''}" \
             "~{if defined(species_id_results_tar) then species_id_results_tar else ''}" \
             "~{if defined(kdr_analysis_results_tar) then kdr_analysis_results_tar else ''}"; do
      [ -n "$t" ] && tar -xzf "$t" -C .
    done

    cp ~{docs_run_info} docs/ampseeker-results/notebooks/run-information.ipynb
    cp ~{docs_snp_df} docs/ampseeker-results/notebooks/snp-dataframe.ipynb
    cp ~{docs_sample_quality_control} docs/ampseeker-results/notebooks/sample-quality-control.ipynb
    cp ~{docs_read_quality} docs/ampseeker-results/notebooks/read-quality.ipynb
    [ -n "~{if defined(docs_run_statistics) then docs_run_statistics else ''}" ] && cp ~{docs_run_statistics} docs/ampseeker-results/notebooks/run-statistics.ipynb || true
    [ -n "~{if defined(docs_reads_per_well) then docs_reads_per_well else ''}" ] && cp ~{docs_reads_per_well} docs/ampseeker-results/notebooks/reads-per-well.ipynb || true
    [ -n "~{if defined(docs_coverage) then docs_coverage else ''}" ] && cp ~{docs_coverage} docs/ampseeker-results/notebooks/coverage.ipynb || true
    [ -n "~{if defined(docs_pop_structure) then docs_pop_structure else ''}" ] && cp ~{docs_pop_structure} docs/ampseeker-results/notebooks/population-structure.ipynb || true
    [ -n "~{if defined(docs_af) then docs_af else ''}" ] && cp ~{docs_af} docs/ampseeker-results/notebooks/allele-frequencies.ipynb || true
    [ -n "~{if defined(docs_sample_map) then docs_sample_map else ''}" ] && cp ~{docs_sample_map} docs/ampseeker-results/notebooks/sample-map.ipynb || true
    [ -n "~{if defined(docs_genetic_diversity) then docs_genetic_diversity else ''}" ] && cp ~{docs_genetic_diversity} docs/ampseeker-results/notebooks/genetic-diversity.ipynb || true
    if [ "~{panel}" = "ag-vampir" ]; then
      [ -n "~{if defined(docs_species_id) then docs_species_id else ''}" ] && cp ~{docs_species_id} docs/ampseeker-results/notebooks/ag-vampir/species-id.ipynb || true
      [ -n "~{if defined(docs_kdr_analysis) then docs_kdr_analysis else ''}" ] && cp ~{docs_kdr_analysis} docs/ampseeker-results/notebooks/ag-vampir/kdr-analysis.ipynb || true
    fi

    papermill /opt/AmpSeeker/workflow/notebooks/process-toc.ipynb process-toc.ipynb -k AmpSeq_python \
      -p wkdir "." -p bcl_convert ~{from_bcl} -p panel "~{panel}" -p plate_info ~{plate_info}

    papermill /opt/AmpSeeker/workflow/notebooks/process-notebooks.ipynb process-notebooks.ipynb -k AmpSeq_python \
      -p wkdir "." -p panel "~{panel}"

    # Symlink results/ inside book source so relative paths from notebooks resolve at build time.
    ln -sfn "$(pwd)/results" docs/ampseeker-results/results

    # Jupyter Book 2 requires Node.js. Do not rely on its interactive installer
    # because Terra/Cromwell has no interactive stdin.
    apt-get update
    DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends nodejs npm
    node --version
    npm --version

    jupyter-book build --all docs/ampseeker-results --path-output results/ampseeker-results

    # index.html is only a redirect stub pointing at intro.html; the real
    # report is the whole _build/html tree (intro.html, per-notebook pages,
    # _static/ assets, search index). Delocalizing a single File would drop
    # all of that, so archive the entire site and return it too.
    # Archive the built site; also include any saved PDFs alongside it.
    mkdir -p ampseeker-results-export
    cp -r results/ampseeker-results/_build/html/. ampseeker-results-export/
    [ -d results/plots ] && cp -r results/plots ampseeker-results-export/plots || true
    tar -czf ampseeker-results-site.tar.gz -C ampseeker-results-export .

    # Produce a single self-contained HTML file by inlining all CSS, JS and
    # images into intro.html — this is what Terra shows as a clickable output
    # that opens directly in the browser without needing to unzip anything.
    pip install --quiet htmlark
    htmlark results/ampseeker-results/_build/html/intro.html \
        -o ampseeker-results-standalone.html || \
        cp results/ampseeker-results/_build/html/intro.html ampseeker-results-standalone.html
  >>>
  output {
    File home_page       = "results/ampseeker-results/_build/html/index.html"
    File site_tar        = "ampseeker-results-site.tar.gz"
    File standalone_html = "ampseeker-results-standalone.html"
  }
  runtime {
    disks: "local-disk ~{disk_size_gb} HDD"
    docker: docker
    cpu: 2
    memory: "8 GB"
  }
}