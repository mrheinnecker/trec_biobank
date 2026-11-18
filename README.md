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
sbatch /g/schwab/marco/repos/trec_biobank/launch_wfTEM.sh
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
|                     #        czi_to_tiff, nd2_to_tiff, lif_metadata, sshpass, microsam,
|                     #        r_analysis (R for annotations/ and analysis/; on WSL2 build it with
|                     #        container/build_r_analysis.sh)
+-- tem/              # wfTEM.nf, main.sh, nextflow.config, process_images.py, correct_gradient.py,
|                     # extract_image_stats.py, make_collection_table.R, imaging_ov.R, figures/
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
+-- annotations/      # launcher.sh, make_split_collection_tables.R, annotation_rounds_common.R,
|                     # migrate_annotation_rounds.R, setup_main_table.R, annotation_ov.R, bio_anno.R
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
| `SELECTORGANELLES` / `ASSIGNORGANELLES` | ASVs that SILVA calls Chloroplast or Mitochondria are classified again against PR2, which has organelle 16S labelled with the eukaryotic host (e.g. `Eukaryota:plas;...;Thalassiosira:plas`) |
| `MERGETAXONOMY` | Combines all into `<ampliseq results>/taxonomy/`: `asv_taxonomy.tsv` (lineage, lowest rank, bootstrap, category: prokaryote / plastid / mitochondrion / eukaryote_nuclear / unassigned, all SILVA and PR2 columns; plastid/mitochondrion ASVs also get `organelle_host_taxonomy`), `asv_table_taxonomy.tsv` (counts per sample with taxonomy and host lineage), `taxonomy_overview.tsv` |

## Overview and annotations

**Biobank overview (`overview/biobank_overview.R`, `overview/biobank_overview_figure.R`).** Run locally in R. Reads the sampling events, HPF, cryo PF/WM, MetaB and chemical-fixation sheets and the collection tables of all modalities, maps every sample and imaged dataset to its sampling event, and writes event maps and coverage tables (`*_event_map.tsv`, `imaging_map.tsv`, `event_coverage.tsv`). The figure script draws the overview figures: circle, radial, HPF flow and the chain figure (samples -> blocks/PF -> imaging, coloured by source). Set `refresh_from_sheets` to re-export the sheets (needs Google auth).

**Collection-table overview (`launch_collection_overview.sh`).** Summarises all collection tables (rows, dataset metadata, counts per modality) in the tidyverse container.

**Annotations (`annotations/`).**
- `annotations/make_split_collection_tables.R` (commands in `annotations/launcher.sh`) distributes the images of the `tem_collection_table` to the annotators and collects their work. Every image is annotated twice by different people: round 1 (name in `annotated_by`) and round 2, the validation (name in `validated_by`). The name counts, not the tab it was entered in. Sheets in the annotation spreadsheet:
  - `main`: one row per finished round (`round`, `annotated_by` = who did it, annotation values); a validated image has two rows.
  - `status`: one row per image with `annotated_by` (round 1) and `validated_by` (round 2).
  - one tab per assignee with the current assignment (`round` 1 or 2; for round 2 `annotated_by` shows the first annotator).

  Every run backs up all sheets to `<annotation_log_dir>/<timestamp>/` before writing anything (plus the new state as `new_*.tsv` and `issues.tsv`). Rows with values and a name move to `main`; rows with values but no name stay in the same tab; unknown names stop the run (known misspellings are in `name_aliases` in `annotations/annotation_rounds_common.R`). The tabs are then refilled: all open images are split equally over `--assignees`, images without any annotation first (contiguous blocks), validations behind them (never to the first annotator); `--max_rows_per_person` optionally caps the rows per person. When someone leaves, move the name to `--former_assignees`: the tab is still backed up and its finished rows collected, then it is emptied; unsigned rows there go back to the pool.
- `annotations/migrate_annotation_rounds.R`: one-off conversion of the old main sheet to the round format, rebuilt from all backups (`--write TRUE` to overwrite `main` and `status`).
- `annotations/setup_main_table.R`: original setup of the main table (historical).
- `annotations/annotation_ov.R`, `annotations/bio_anno.R`: early overview plots of the quality and biological annotations in the old `image_log` sheet (historical).
- `analysis/` (repo root) holds the annotation analysis, one script per analysis. Input is the `main` sheet of the annotation spreadsheet, read with googlesheets4 (`--sheet_url`, `--sheet`, `--google_key`), or a tsv snapshot of it via `--input`. All scripts take `--outdir <dir>` and write plots to `<dir>/plots` and tables to `<dir>/tables`. `--rounds` selects the annotation rounds (default `1`: only the first annotation of every cell; `1,2` also the validations). Shared data preparation is in `common.R`: it selects the rounds, counts uncertain entries (`?`) as absent (0/1 columns: `??`, `1?` -> 0; `undescribed`: numbers with `?` are dropped; categorical columns unchanged) and splits the `undescribed` column into one 0/1 column per number.
  - `00_download_annotations.R`: saves the `main` sheet as `<dir>/tables/annotations_main.tsv` (plus a timestamped copy), so all steps of a run use the same data.
  - `01_annotator_qc.R`: beauty scores per annotator; overview of the uncertain entries (`plots/uncertain_entries`, `tables/uncertain_entries.tsv`); agreement between annotation and validation.
  - `02_site_overview.R`: annotated cells per site and major group.
  - `03_site_features.R` (dinoflagellates): per-site frequency of the organelles in `use_features` (`assoc_common.R`) and of every undescribed structure, with site tests (Kruskal-Wallis, chi-squared, binomial GLMM with sample as random effect). The top bar is split into the samples (barcodes), the site labels give `samples / sampling dates`; `tables/site_samples.tsv` lists every sample with date, size fraction, time and number of cells.
  - `04_organelle_associations.R` (dinoflagellates): pairwise association tests between the selected organelles (`use_features` + undescribed, set in `assoc_common.R`; features present in fewer than `min_feature_cells` = 10 cells are dropped for all of 04-09), Bonferroni within this analysis. `05_organelle_association_plots.R` plots the significant pairs and needs the output of `04`.
  - `09_negative_pair_collections.R`: for every significant negative organelle pair of `04` a MoBIE collection table of the cells split into `<A>_only`, `<B>_only`, `both`, `none` (one view/grid each); written to `tables/collections/` and as one tab per pair (+ `overview`) to the spreadsheet in `TEM_COLLECTION_SHEET_URL`.
  - `10_collection_montages.py`: one JPEG overview per collection table of `06` (clusters, groups C1..Ck) and `09` (`montages_f<factor>/<A>_vs_<B>.jpg`): the four groups stacked under header bands, all cells true to scale (all overview PNGs have the same pixel size, 7.06 nm; all reduced by `TEM_SCALE_FACTOR`, default 4 -> 28.3 nm/px, 10 um scale bar in the title), packed with the skyline bottom-left heuristic into an about square figure, each cell labelled with its name. Images from `TEM_PNG_DIR/<site>/<source_name>.png`, reduced once to grayscale JPEGs cached in `TEM_THUMB_DIR`. For the clusters additionally `clusters_k<k>_tree.jpg`: the clustering tree on top (branches coloured by cluster), below it the metadata tracks and the organelle presence/absence per cell as in the clustering heatmap (from `tables/clust_tree_k<k>_*.tsv` of `06`), and one column of cells per cluster, aligned with its branch. Cells with a `taxo_class` get a coloured frame (one colour per class) and a `taxo_class` track below the organelles. Column width proportional to the number of cells of the cluster (equal strip per cell); the cell images are reduced (same factor for all) until the widest cell of each cluster fits its column (sizes: `TEM_TREE_PLOT_SCALE`, upper limit of the image size `TEM_TREE_IMAGE_SCALE`).
  - `11_chloro_foodvac_sites.R` (dinoflagellates): three groups from `chloropl` and `food_vac` (chloroplast + food vacuole, chloroplast only, no chloroplast) and their fraction per site (A) and per sample (B, date order within the site). Site tests: Fisher site x group (cells independent, Cramer's V) and per group a binomial GLMM with sample as random effect (LRT, BH). Tables `tables/chloro_foodvac_*.tsv`.
  - `12_chloro_foodvac_map.R`: the groups of `11` as one stacked bar per site on a map of Europe (orthographic projection, Natural Earth countries, as in `metab/analysis/ov_map_plot.R`), each bar standing on the site's coordinates (`site_coords` in the script). Needs `tables/chloro_foodvac_site_counts.tsv` of `11`.
  - `13_trophic_map_combined.R`: imaging and sequencing on one map: per site the TEM bar of `11` (black outline) next to one metabarcoding bar per size fraction (`metab_fractions`, default 5-40, 10-40, 10-200, 20-200) with the Tara W4 trophic modes of the dinoflagellate reads (table of `metab/analysis/06_trophic_modes_map.R`, `TEM_METAB_TROPHIC_TABLE`). Sites without TEM cells only get the metabarcoding bars. Second figure `trophic_map_combined_by_fraction`: TEM also split by size fraction (per site and size fraction a pair of bars TEM | metaB), the sampling month under each site name, and on the right the sampling time span per site, separately for TEM (dates of the annotated cells) and metaB (dates of the samples of the shown size fractions, `05_trophic_mode_per_sample.tsv`, `TEM_METAB_SAMPLE_TABLE`); dates also in `tables/trophic_map_sampling_dates.tsv`.
  - `14_tem_feature_maps.R` (TEM only): the same map layout as the by-fraction figure of `13` (one bar per site and size fraction, sampling month, sampling dates on the right) for two other groupings of the dinoflagellate cells: crystalline reticulate network (`cryst_rich_RN`) present / absent (`plots/tem_map_cryst_rich_RN`) and starch + lipid droplets / starch only / lipid droplets only / none (`plots/tem_map_starch_lipid`). Groupings are defined in the `groupings` list at the top. The map code shared by `13` and `14` (site coordinates and bar shifts, base map, bars, legends, timeline) is in `map_common.R`.
  - `07_metadata_associations.R` (dinoflagellates): every selected organelle against the metadata (site, size fraction, sampling time, life status, beauty), separate from 04 with its own Bonferroni correction. `08_metadata_association_plots.R`: heatmap of all organelle x metadata pairs, composition plots and contingency tables of the significant ones, and `assoc_metadata_overview`: everything in one figure (rows = organelles, columns = every category of every metadata variable except site, which has its own figure in `03`; number = % of the cells of the category with the organelle, colour = difference to the overall prevalence, categories with < 5 cells not coloured; per variable a `test` column with Cramer's V / Spearman rho and stars and a `p_adj` column; values in `tables/assoc_metadata_prevalence.tsv`); needs the output of `07`.
  - `06_organelle_clustering.R` (dinoflagellates): clusters the cells by organelle presence/absence (Jaccard distance, Ward linkage, k picked by silhouette width or fixed with `TEM_CLUSTER_K`; settings at the top of the script). Also writes a collection table with one view/grid per cluster (`tables/collections/clusters_k<k>.tsv`, tab `clusters_k<k>` in the collection spreadsheet). Outputs a heatmap, a PCoA, the organelle profile of each cluster, cluster composition per site, size fraction, sampling time and life status, and the taxonomic class (`taxo_class`, where filled) per cluster (`plots/clust_taxo_class`, `tables/clust_taxo_class_counts.tsv`, `tables/clust_taxo_class_cells.tsv`).
  - `launcher.sh`: all settings for a local run (container, key, rounds, steps, output folder) in one place; edit and run with `bash analysis/launcher.sh`.
  - `run_analysis.sh <outdir> [steps]` creates a new folder `<outdir>/<timestamp>_round<ROUNDS>/` per run, downloads the snapshot (00) into it and runs the steps on it in an R container (`R_SIF`, which needs tidyverse, googlesheets4, cowplot, patchwork, lme4, cluster and ComplexHeatmap). Env: `ROUNDS` (default 1), `INPUT` (existing snapshot instead of downloading), `TEM_GOOGLE_KEY`.
  - `01_annotator_qc.R`: the annotation vs validation comparison only runs with `--rounds 1,2`.

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
