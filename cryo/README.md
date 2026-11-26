# CRYO table-driven workflow

This workflow reads listed light-microscopy Z-stack files, extracts pixel-size metadata, prepares each raw file for conversion, converts it to OME-Zarr with EuBI-Bridge, copies completed datasets into the mounted S3 share (`--publish_dir`), and writes a collection table.

Before selecting datasets, the workflow lists the top level of `--publish_dir` (the S3 bucket mounted on `/g`, default `/g/schwab/s3-public/imatrec/central_data_processing/cryo`). Any dataset with a complete `<dataset_name>.ome.zarr/` (root `zarr.json` or `.zattrs` present) is excluded from `images_to_process.csv`, independently of the Nextflow `-resume` cache; incomplete ones are listed in `incomplete_s3_entries.txt` and processed again.

For a non-dry-run `all` or `collection` run, the workflow stops if `--publish_dir` is not accessible. This avoids accidentally reprocessing every dataset or writing an incomplete collection table when S3 is unavailable.

## Input table

In contrast to the other workflows, the data for cryo and plastic is already registred in 
labID, which means we rely on thsi as out main input source. This also measn that we need to query labID first... this needs password authentification which is currently not supported in the labid command line tool
So whenever new data is there we need to query via this command the table in labID, copy it to a google sheet (https://docs.google.com/spreadsheets/d/1ePRpa56mmMvCeRTLXmwOywOLy5_I3AFrxJepSUYGR1s/edit?gid=1442254503#gid=1442254503)
and then rerun the workflow.

The input table can be a Google Sheet or local TSV/CSV. It should contain one raw-file column. Accepted names are:

```text
raw_path
file_path
filepath
file
source_path
path
```

Example:

```text
raw_path
/g/example/light_microscopy/sample_001.ome.tif
```

If a `convert` column is present, values such as `0`, `FALSE`, or `no` keep the row in `all_datasets.tsv` but skip processing. If `convert` is absent, rows are processed by default unless already present in S3.

Dataset names are derived from the raw filename unless a `filename` column is present.

## Pixel-size metadata

The workflow needs X pixel size, Y pixel size, and Z spacing in nanometers for EuBI-Bridge conversion.

For CZI inputs, `extract_metadata.py` reads Zeiss scaling metadata from the CZI XML and converts meter-scale values to nanometers.

For TIFF/OME-TIFF inputs, `extract_metadata.py` tries to read:

- OME `PhysicalSizeX`, `PhysicalSizeY`, and `PhysicalSizeZ`
- ImageJ `spacing`
- TIFF resolution tags for X/Y when available

You can override or provide missing values in the table with:

```text
x_scale
y_scale
z_scale
scale_unit
```

Aliases such as `x_pixel_size_nm`, `y_pixel_size_nm`, and `z_spacing_nm` are also accepted. `scale_unit` supports `nm`, `um`, and `mm`; the workflow converts values to nanometers before calling EuBI-Bridge.

Launcher defaults are also available:

```bash
bash cryo_main.sh interactive \
  --default_x_scale 250 \
  --default_y_scale 250 \
  --default_z_scale 1000 \
  --scale_unit nm
```

## CZI preparation

Zeiss `.czi` files are not sent directly to EuBI-Bridge. The workflow first runs `prepare_input.py`, which converts each CZI into an intermediate OME-TIFF inside the Nextflow work directory. EuBI-Bridge then converts that OME-TIFF to OME-Zarr.

Non-CZI inputs are linked or copied into the same prepared-input directory and then passed to EuBI-Bridge unchanged. The original raw files are not modified.

The intermediate OME-TIFF stores the image axes, physical X/Y/Z pixel sizes, and channel names where available. After EuBI-Bridge conversion, `patch_omezarr_metadata.py` updates the final OME-Zarr root `.zattrs` with `omero.channels` display metadata and a compact `cryo_metadata` block derived from the original CZI metadata JSON. For CZI inputs, the full raw Zeiss metadata XML is also copied into the OME-Zarr under `czi_metadata/czi_metadata.xml`.

## Run examples

Discover selected rows without conversion:

```bash
cd cryo
bash cryo_main.sh local \
  --sheet_mode local \
  --input_table /path/to/cryo_images.tsv \
  --workflow_stage discover
```

Run conversion without upload:

```bash
cd cryo
bash cryo_main.sh interactive \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/..." \
  --sheet_name "cryo_input_table" \
  --workflow_stage process
```

Run conversion, upload, and write the collection table:

```bash
cd cryo
bash cryo_main.sh cluster \
  --sheet_mode google \
  --sheet_url "https://docs.google.com/spreadsheets/d/..." \
  --sheet_name "cryo_input_table" \
  --collection_table_url "https://docs.google.com/spreadsheets/d/..." \
  --collection_table_sheet "cryo_collection_table" \
  --s3_bucket "s3embl/imatrec/central_data_processing/cryo" \
  --workflow_stage all
```

Rebuild only the collection table from S3 and persisted metadata:

```bash
bash cryo_main.sh cluster --workflow_stage collection
```

## Workflow stages

- `discover`: parse the table and write `images_to_process.csv` / `all_datasets.tsv`.
- `process`: extract metadata, prepare inputs, and convert listed raw files to OME-Zarr.
- `all`: extract metadata, prepare inputs, convert, copy into the S3 share, and write the collection table.
- `collection`: list the S3 share and rebuild the collection table from the input table plus persisted metadata.

Metadata JSON and pixel-size TSV files are persisted under `--persistent_metadata_dir`, defaulting to `<main_dir>/metadata` from the launcher or the profile-specific central-data path on the cluster.

With `--workflow_stage all`, datasets are copied from the work directory into the mounted S3 share, first into a hidden `.<dataset_name>.ome.zarr.partial` directory that is renamed when the copy is complete:

```text
<publish_dir>/<dataset_name>.ome.zarr/
```

`--s3_bucket` only sets the public URLs (`https://schwab.buckets.embl.de/<bucket path without the mc alias>/...`). The collection table contains only datasets that are in the input table, and it is rebuilt also when no new image is processed. Grid positions are counted per site.

No image data is written outside the Nextflow work directory except the final OME-Zarr in `--publish_dir`; the log directory and `--persistent_metadata_dir` only receive metadata and logs.

HPF carrier screens (`_HPF_` in the file name) are written to a separate collection table, `--hpf_collection_table_sheet` (default `cryo_hpf_collection_table`, locally `cryo_hpf_collection_table.tsv`); PF and WM stay in `--collection_table_sheet`. Both are written by the same run from the same input table. The HPF CZIs are unstitched tile scans; the tiles are placed at their stage positions when the CZI is read in `prepare_input.py` (bioio-czi), without registration or blending.

The input table only needs a path column (`path`, `FilePath`, `raw_path`, ...). Pixel sizes, channels and the site come from the CZI and the file name.

The collection table keeps `name` as `<dataset_name>` but points `uri` to the `.ome.zarr` prefix. Existing legacy `.zarr` datasets are still recognized.

For multi-channel CZI data, the collection table is expanded to one row per channel. The original dataset name is kept in `source_name`, while `name` becomes channel-specific, for example:

```text
sample_c0_DAPI
sample_c1_GFP
```

The table also writes `channel`, `channel_label`, `display`, `color`, `blend=sum`, `format=OmeZarr`, `grid`, and `grid_position`. Channels from the same image share the same `grid_position`, and `display` is channel-specific so MoBIE keeps independent display controls for DAPI/GFP/etc.

During CZI preparation, the workflow computes per-channel statistics from the `TCZYX` image data: min/max (data range) and the 0.1/99.5 percentiles (display range). The percentiles are stored as `contrast_limits` in the collection table and as the `omero` window `start`/`end`, because min/max is useless when a few pixels are saturated.

All channels are kept in their original order, so `omero.channels[i]`, the OME-XML channel `i` and the collection-table `channel` column all describe image channel `i`.

Optional input-table columns `site` and `grid` set the MoBIE `view` and `grid` of a dataset (`grid` defaults to `site`). Without them the site is the leading letters of the dataset name, e.g. `Roscoff` in `Roscoff_HPF_...`.

`build_10x_input_table.py` builds an input table with the 10x overview of every grid from a `find` inventory of the 2021-Villefranche / 2022-Iceland folders (see the script header).

## Containers

The CRYO workflow uses two separate containers:

- CZI metadata/preparation: `/g/schwab/marco/container_devel/czi_to_tiff.sif`
- EuBI-Bridge conversion: `/g/schwab/marco/container_devel/eubibridge.sif`

Build the CZI preparation container from:

```text
container/czi_to_tiff.def
```

For example:

```bash
singularity build /g/schwab/marco/container_devel/czi_to_tiff.sif container/czi_to_tiff.def
```

This keeps the existing EuBI-Bridge container unchanged and only adds the Zeiss CZI dependencies where they are needed.
