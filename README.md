# TREC-Biobank / IMATREC central data processing workflows

This repository contains code to process all data modalities of the TREC-Biobank / IMATREC project by the Schwab team, CBB unit EMBL Heidelberg.

The repo is sub-structured according to the various modalities and downstream analyses:

| Directory | Content |
|---|---|
| `tem/` | Transmission electron microscopy (TEM) |
| `sem/` | Scanning electron microscopy (SEM) |
| `hitt/` | High-Throughput Tomography; X-ray scans |
| `cryo/` | Cryo confocal images of plunge-frozen (PF) and waffle-method (WM) grids, and HPF carrier screens |
| `plastic/` | Confocal images of freeze-substituted resin blocks (Leica `.lif`) |
| `exm/` | Expansion microscopy (Nikon `.nd2`, one file per channel) |
| `metab/` | Metabarcoding (16S/18S amplicons): nf-core/ampliseq plus SILVA/PR2 taxonomy |
| `overview/` | Biobank overview (samples, sampling events, imaging per modality) and collection-table overview |
| `annotations/` | TEM annotation table setup |
| `analysis/` | TEM annotation analysis (annotator QC, site features, feature associations, organelle clustering) |
| `microsam/` | micro-SAM segmentation, DINOv2 embedding and clustering of TEM images (separate pipeline, see `microsam/README.txt`) |
| `container/` | Apptainer/Singularity definition files |

The imaging workflows follow the same strategy: starting from raw data as it arrives from the microscopes / input sources and ending in conversion to OME-Zarr (zarr format 3 by default), storing everything in S3 storage and creating collection tables to browse the data through MoBIE.
Each workflow is implemented in Nextflow and has its dependencies handled via Apptainer/Singularity containers (see definition files in `container/`).

## Usage instructions

The current implementation is optimized for the EMBL Slurm HPC environment. Every imaging workflow has a launcher in the repo root (`launch_wf<MODALITY>.sh`) that holds the complete run definition: input sheet, collection table, paths, S3 target and workflow stage. Edit the launcher, then submit it:

```bash
sbatch /g/schwab/marco/repos/tem_classification/launch_wfTEM.sh
```

The launcher runs as one Slurm job (2 CPUs, 16 GB, 24 h by default, logs in `/scratch/rheinnec/<modality>_screen/`). The first argument of the `*_main.sh` call inside it selects where the Nextflow tasks run:

- `cluster`: each Nextflow process is submitted as its own Slurm job (except light-weight guardrail or file-collection steps).
- `interactive`: all processes run inside the launcher job / current node.
- `local`: lightweight discovery/debugging on local paths.

Recommended settings:

| Launcher | Mode | Notes |
|---|---|---|
| `launch_wfTEM.sh` | `cluster` | Heavy (IMOD blending, gradient correction, conversion) |
| `launch_wfHITT.sh` | `cluster` | Extremely heavy. Export the Hamburg cerberus password first: `export HITT_SSHPASS='PASSWORD'`; allow at least 48 h |
| `launch_wfSEM.sh` | `interactive` | Lightweight |
| `launch_wfCRYO.sh` | `interactive` | Lightweight; input from the `cryo_input_table` sheet |
| `launch_wfPLASTIC.sh` | `interactive` | Lightweight |
| `launch_wfEXM.sh` | `interactive` | EuBI memory is sized per image from the input size |
| `launch_collection_overview.sh` | - | Collection-table overview (see [Overview](#overview-and-annotations)) |
| `metab/launch_ampliseq_cluster.sh` | - | Run with `bash` on the login node; it submits the ampliseq job itself (see [Metabarcoding](#metabarcoding-metab)) |

Common launcher options:

- `--workflow_stage`: `discover` (only find images and write logs), `process` (process and convert, no S3/collection table), `all` (full run including S3 and collection table); CRYO also has `collection` (only rebuild the collection table).
- `--dryrun TRUE --dryrun_n N`: process only N candidates (EXM: `--dryrun_images` picks fixed images).
- `--resume`: reuse the Nextflow cache.
- `--logdir`: run logs, Nextflow trace/report/timeline/DAG and per-dataset metadata go to `/g/schwab/biobank_logs/<modality>/wf<MODALITY>_<timestamp>`.

### S3 storage

All datasets go to `imatrec/central_data_processing/<modality>` on the EMBL S3 (`schwab.buckets.embl.de`).

- TEM, SEM and HITT upload with the MinIO client (`mc`). The root `zarr.json` is uploaded last, so a dataset only counts as present (and is skipped in later runs) once its upload is complete.
- CRYO, PLASTIC and EXM copy into the mounted S3 share (`--publish_dir`, `/g/schwab/s3-public/imatrec/central_data_processing/<modality>`). The copy goes to a hidden `.<name>.ome.zarr.partial` folder that is renamed when complete; incomplete entries are listed in `incomplete_s3_entries.txt` and processed again. Non-dry `all`/`collection` runs stop if the share is not mounted.

## Repository Layout

```text
.
+-- launch_wf{TEM,SEM,HITT,CRYO,PLASTIC,EXM}.sh   # run definitions, submit with sbatch
+-- launch_collection_overview.sh
+-- container/        # *.def: eubibridge, imod(_procps), mc, py_mrcfile, segmentation,
|                     #        czi_to_tiff, nd2_to_tiff, lif_metadata, sshpass, microsam
+-- tem/              # wfTEM.nf, main.sh, nextflow.config, process_images.py, correct_gradient.py,
|                     # extract_image_stats.py, make_collection_table.R, make_split_collection_tables.R,
|                     # imaging_ov.R, annotation_ov.R, figures/
+-- sem/              # wfSEM.nf, sem_main.sh, nextflow.config, select_images.R, extract_metadata.py,
|                     # make_collection_table.R
+-- hitt/             # wfHITT.nf, hitt_main.sh, nextflow.config, prepare_slices.py,
|                     # extract_stack_stats.py, make_collection_table.R
+-- cryo/             # wfCRYO.nf, cryo_main.sh, nextflow.config, extract_metadata.py, prepare_input.py,
|                     # patch_omezarr_metadata.py, make_collection_table.R, cryo_ov_figure.R
+-- plastic/          # wfPLASTIC.nf, plastic_main.sh, nextflow.config, extract_metadata.py,
|                     # patch_omezarr_metadata.py, make_collection_table.R
+-- exm/              # wfEXM.nf, exm_main.sh, nextflow.config, prep_input_table.R, extract_metadata.py,
|                     # prepare_input.py, patch_omezarr_metadata.py, make_collection_table.R
+-- metab/            # ampliseq launch scripts, samplesheet builder, datahub download, AMPLISEQ_SETUP.md
|   +-- taxonomy/     # wfTAXONOMY.nf + R scripts (SILVA/PR2 classification)
+-- overview/         # biobank_overview.R, biobank_overview_figure.R, collection_table_overview.R
+-- annotations/      # setup_main_table.R, launcher.sh
+-- analysis/         # annotation analysis: common.R, 01-05 scripts, run_analysis.sh
+-- microsam/         # segmentation / embedding / clustering pipeline
```

Each imaging module has its own README (`tem/`, `hitt/`, `cryo/`, `plastic/`, `exm/`) with the detailed input-table format and processing steps.

## TEM Workflow

Defined in `tem/wfTEM.nf`, configured by `tem/nextflow.config`, launched through `tem/main.sh` (via `launch_wfTEM.sh`).

1. Scans the raw TEM directory for `.mrc` files and matching `.mdoc` files.
2. Maintains an image log in a local TSV file or a Google Sheet.
3. Selects images that are not yet completely uploaded to S3 (root `zarr.json` present).
4. Runs IMOD tools (`justblend`, `blendmont`) for montages.
5. Detects and optionally corrects broad low-frequency grayscale gradients (metrics/QC to the log dir).
6. Generates overview PNGs and intensity statistics for contrast limits.
7. Converts processed MRC images to OME-Zarr using EuBI-Bridge.
8. In `all` mode, uploads to S3 and writes the `tem_collection_table`.

```bash
bash tem/main.sh cluster --sheet_mode google --workflow_stage all --dryrun TRUE --dryrun_n 2 [...]
```

Options can also be given as environment variables (`WORKFLOW_STAGE`, `DRYRUN`, `RESUME`, `SHEET_MODE`, `RAWDIR`, ...). See `tem/README.md`.

## HITT Workflow

Defined in `hitt/wfHITT.nf`, launched through `hitt/hitt_main.sh`.

HITT does not scan a raw directory. It reads a Google Sheet or local table with a remote `source_path` column (Hamburg cerberus server) and processes exactly those image directories:

1. Incrementally copies each remote stack into scratch (`/scratch/rheinnec/hitt_screen/...`).
2. Stages a renamed copy of the `tomo/slice_*.tif` stack so it starts at `Z0001.tif`, converted to `uint16` by default.
3. Calculates stack-wide gray-value display limits.
4. Converts to OME-Zarr, uploads to `s3embl/imatrec/central_data_processing/hitt/<image_name>/` and updates the `hitt_collection_table`.

See `hitt/README.md`.

## SEM Workflow

Defined in `sem/wfSEM.nf`, launched through `sem/sem_main.sh`.

1. Scans the raw SEM directory (`/g/schwab/Chandni/SEM/IMATREC SEM`) for `.tif`/`.tiff` files.
2. Builds a local or Google-backed image log.
3. Extracts TIFF/Zeiss SEM metadata including pixel size (kept in the log dir, no intermediate data published).
4. Converts to OME-Zarr, uploads to S3 and builds the SEM collection table.

## CRYO Workflow

Defined in `cryo/wfCRYO.nf`, launched through `cryo/cryo_main.sh`.

The data is registered in labID; the labID table is copied into the `cryo_input_table` Google Sheet (paths only), which is the workflow input. For every listed file the workflow extracts X/Y/Z pixel size and channels, prepares the input (Zeiss `.czi` is converted to an intermediate OME-TIFF in the `czi_to_tiff` container), converts with EuBI-Bridge, patches the OME-Zarr metadata (channel names/colors) and publishes into the mounted S3 share.

- PF and WM grids go to `cryo_collection_table`.
- HPF carrier screens (`_HPF_` in the file name) go to a separate `cryo_hpf_collection_table` (`--hpf_collection_table_sheet`). These are unstitched tile scans; the tiles are placed at their stage positions, without registration or blending.
- Channels are read from the acquired image dimensions; the reflection channel is shown in white.

See `cryo/README.md`.

## PLASTIC Workflow

Defined in `plastic/wfPLASTIC.nf`, launched through `plastic/plastic_main.sh`.

Same table/S3/collection-table pattern as CRYO, for Leica `.lif` files. There is no intermediate OME-TIFF: the raw `.lif` goes directly to EuBI-Bridge (extra flags via `--eubi_extra_args`). See `plastic/README.md`.

## EXM Workflow

Defined in `exm/wfEXM.nf`, launched through `exm/exm_main.sh` (started as a copy of CRYO).

```text
CHECKEXISTINGEXMS3FILES -> SELECTEXMIMAGES -> EXTRACTEXMMETADATA -> PREPAREEXMINPUT
  -> EUBIEXMCONVERSION -> PATCHEXMOMEZARRMETADATA (patch + publish) -> COLLECTEXMS3FILES -> MAKEEXMCOLLECTIONTABLE
```

- Input: `exm_input_table` (one row per `.nd2` file or folder). All channel files of one image (`Img0098_Channel*.nd2`) are merged into one multi-channel OME-Zarr, ordered DAPI, FITC, TRITC, Cy5.
- Datasets are named by `shortname` (e.g. `BAR_Img0093`); the site prefix is used as view/grid in the collection table.
- Channel names in the collection table are the stained structures (FITC Centrin, TRITC Tubulin, Cy5 NHS), with MoBIE colors per channel.
- EuBI memory is 4x the prepared OME-TIFF size + 8 GB (at least `req_mem`); time limit 2 h, rerun timed-out images with `--eubi_time 6h`.
- The collection table only contains datasets from the input table and is rebuilt also when no new images were converted.

See `exm/README.md`.

## Metabarcoding (`metab/`)

16S/18S metabarcoding of the TREC samples (primers 515F-Y / 926R, which amplify both prokaryotic/organellar 16S and eukaryotic 18S). Details, run structure and rationale are in `metab/AMPLISEQ_SETUP.md`.

**Data.** Raw FASTQs come from the EMBL datahub (`metab/datahub_links.txt`, download with `sbatch metab/download_datahub_slurm.sh`). `metab/sampleinfo.tsv` maps the files to samples. There are four sequencing runs (flowcells): three from 2023 (no controls) and one from 2024 with 6 extraction/PCR blanks.

**Step 1: nf-core/ampliseq 2.18.0.** Edit the paths at the top of `metab/launch_ampliseq_cluster.sh`, then on the login node:

```bash
bash metab/launch_ampliseq_cluster.sh
```

It builds the samplesheet (`make_ampliseq_samplesheet.sh`: skips incomplete downloads, adds the flowcell as `run` so DADA2 learns one error model per run), writes the params file and a `resourceLimits` config, and submits `run_ampliseq_slurm.sh` as the Nextflow controller job. Pipeline settings:

- `ampliseq_orientation.config` adds Cutadapt `--revcomp`, because the library contains both read orientations.
- `mergepairs_strategy: consensus`: overlapping 16S pairs are merged, the longer non-overlapping 18S pairs are concatenated.
- `skip_taxonomy: true`, taxonomy is done in step 2.
- Blanks that lose all reads are dropped (`ignore_failed_trimming/filtering`) instead of failing the run.

`metab/launcher.sh` runs a small local pilot without Slurm (subsampled reads).

**Step 2: taxonomy (`metab/taxonomy/wfTAXONOMY.nf`).** Download the references once (`wget` commands in `AMPLISEQ_SETUP.md`: SILVA 138.2 and PR2 5.1.1 in DADA2 format, plus the DADA2 container as `.sif`), then:

```bash
sbatch metab/taxonomy/run_taxonomy_slurm.sh
```

| Process | What it does |
|---|---|
| `SPLITASVS` | Splits ASVs by barrnap's best hit: bacteria/archaea/mitochondria -> 16S, eukaryotes -> 18S; ASVs without a hit go to both |
| `ASSIGNTAXONOMY` | DADA2 `assignTaxonomy` (minBoot 50): 16S against SILVA (to genus), 18S against PR2 (to species), in parallel |
| `MERGETAXONOMY` | Combines both into `<ampliseq results>/taxonomy/`: `asv_taxonomy.tsv` (lineage, lowest rank, bootstrap, category: prokaryote / plastid / mitochondrion / eukaryote_nuclear / unassigned, all SILVA and PR2 columns), `asv_table_taxonomy.tsv` (counts per sample with taxonomy), `taxonomy_overview.tsv` |

## Overview and annotations

**Biobank overview (`overview/biobank_overview.R`, `overview/biobank_overview_figure.R`).** Run locally in R. Reads the sampling events, HPF, cryo PF/WM, MetaB and chemical-fixation sheets and the collection tables of all modalities, maps every sample and imaged dataset to its sampling event, and writes event maps and coverage tables (`*_event_map.tsv`, `imaging_map.tsv`, `event_coverage.tsv`). The figure script draws the overview figures: circle, radial, HPF flow and the chain figure (samples -> blocks/PF -> imaging, coloured by source). Set `refresh_from_sheets` to re-export the sheets (needs Google auth).

**Collection-table overview (`launch_collection_overview.sh`).** Summarises all collection tables (rows, dataset metadata, counts per modality) in the tidyverse container.

**Annotations (`annotations/`).**
- `setup_main_table.R` builds the main annotation table from the TEM collection table and the existing annotations.
- `analysis/` (repo root) holds the annotation analysis, one script per analysis. All scripts take `--input <annotations.tsv> --outdir <dir>` and write plots to `<dir>/plots` and tables to `<dir>/tables`. Shared data preparation is in `common.R`: it keeps the first annotation of each cell across annotation runs and splits the `undescribed` column into one 0/1 column per number.
  - `01_annotator_qc.R`: beauty scores per annotator; agreement between annotation and validation.
  - `02_site_overview.R`: annotated cells per site and major group.
  - `03_site_features.R` (dinoflagellates): per-site frequency of every organelle and undescribed structure, with site tests (Kruskal-Wallis, chi-squared, binomial GLMM with sample as random effect).
  - `04_associations.R` (dinoflagellates): pairwise association tests between all features and metadata. `05_association_plots.R` plots the significant pairs and needs the output of `04`.
  - `06_organelle_clustering.R` (dinoflagellates): clusters the cells by organelle presence/absence (Jaccard distance, Ward linkage, k picked by silhouette width; settings at the top of the script). Outputs a heatmap, a PCoA, the organelle profile of each cluster, and cluster composition per site, size fraction, sampling time and life status.
  - `run_analysis.sh <annotations.tsv> <outdir> [steps]` runs the steps in an R container (`R_SIF`, which needs tidyverse, cowplot, patchwork and lme4).
- `annotations/launcher.sh` is temporary until all old annotations are adapted; use with care.

## External Dependencies

- Nextflow; Slurm for cluster runs; Apptainer/Singularity
- R with tidyverse and googlesheets4 (discovery, logs, collection tables, overview, annotations)
- Python image libraries (`mrcfile`, `numpy`, `scipy`, `scikit-image`, `tifffile`, `nd2`, `bioio`/`aicspylibczi`)
- IMOD (TEM), EuBI-Bridge (OME-Zarr conversion), MinIO client `mc` (TEM/SEM/HITT upload)
- nf-core/ampliseq 2.18.0 and DADA2 1.38 (metab), SILVA 138.2 and PR2 5.1.1 reference files
- Google service-account key (`--google_key`) for `sheet_mode=google`

## Notes and Caveats

- Paths, bucket names and Google Sheet URLs are project-specific and live in the launchers and `nextflow.config` files.
- `tem/crop_omezarr_by_mask.py` is experimental/stale and not wired into the TEM workflow.
- `microsam/` uses its own absolute cluster paths and environment (see `microsam/README.txt`).

## Development Notes

Before changing workflow behavior, check:

- The relevant `*.nf` process definitions.
- The matching `nextflow.config` profile and container assignments.
- Whether the script is used by the workflow or is a downstream/manual analysis helper.
- Whether the change affects local, interactive, and cluster modes differently.
