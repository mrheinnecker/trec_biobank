# nf-core/ampliseq setup for the mixed-domain SSU dataset

This setup targets nf-core/ampliseq 2.18.0 with paired Illumina reads that
contain both primer orientations. It does not assume that the amplicons are
exclusively bacterial 16S; the experiment is expected to recover prokaryotic
16S and eukaryotic 18S SSU sequences.

## Required information

Before running, supply:

1. Exact forward and reverse biological primer sequences, 5' to 3'. Exclude
   adapters, indices, pads, linkers, heterogeneity spacers, and barcodes.
2. Cluster-visible absolute paths for FASTQs, output, work, and the Apptainer
   cache.
3. Any required Slurm account, partition, or QoS.
4. Experimental metadata and biological sample names, if available. The
   supplied `sampleinfo.tsv` identifies six paired negative-control samples;
   ENA run accessions are used as sample IDs until better metadata is supplied.

## Why the custom config is required

In a 100,000-pair check of ERR16149169, 45.55% of pairs had the forward primer
in read 1 and 47.72% had it in read 2. The stock ampliseq paired-Illumina
Cutadapt command accepts only the first orientation and discards the other.
`ampliseq_orientation.config` adds Cutadapt 5.2's `--revcomp` option, which
normalizes paired orientation by swapping mates when appropriate.

## Pilot procedure

There are two entry points:

- `launcher.sh`: a minimal local test. It writes its own samplesheet for a
  fixed list of samples and runs Nextflow in the foreground with
  `nextflow.local.config`. It does not use `main.sh`.
- `main.sh`: the cluster run. It validates its arguments, gzip-tests and
  selects the requested samples, creates the samplesheet and parameter file,
  prepares scratch directories, and submits `run_ampliseq_slurm.sh` as the
  Nextflow controller job. Use `all` for either sample count to include every
  available sample of that type. All paths, including the orientation config,
  are passed to the controller job via `sbatch --export`, because sbatch runs
  the script from a copy in the Slurm spool directory.

## Sequencing runs and controls

With `sampleinfo.tsv`, `make_ampliseq_samplesheet.sh` writes a `run` column
(the Illumina flowcell from the submitted filename) and a `control` column.
ampliseq uses `run` to learn one DADA2 error model per sequencing run.

| Dataset | Flowcell | Samples | Controls |
|---|---|---|---|
| 2023 (`MetaB_16SV4V5_Plankton_AML40um`) | HGH7NDRX5 | 105 | 0 |
| | HJ5N5DRX5 | 8 | 0 |
| | HTFWJDRX5 | 7 | 0 |
| 2024 (`..._2024` and `analysis-1`) | HCJF5DRX7 | 92 | 6 |

The six negative controls (one extraction blank, five PCR/tag water blanks)
were all sequenced with 2024. Index cross-talk and reagent contamination are
run- and batch-specific, so these controls are not used to remove features
from the 2023 runs. ampliseq's `decontam` is not run per run, so it stays off
(`decontam: none`, the default); decontamination is done downstream on the
2024 samples, and the 2024 blanks only serve as a check on 2023.

The pilot deliberately sets `skip_taxonomy: true`. First inspect Cutadapt read
retention, DADA2 filtering/merging, ASV lengths, FastQC, and MultiQC. Taxonomic
classification should then be configured with a reference that represents both
prokaryotic 16S and eukaryotic 18S; the pipeline's built-in DADA2 SILVA dataset
is documented as unsuitable for Eukaryotes, while PR2 targets 18S.

## Taxonomy (`taxonomy/wfTAXONOMY.nf`)

ampliseq runs with `skip_taxonomy: true`, because no single DADA2 reference
covers both 16S and 18S. With `mergepairs_strategy: consensus`, 16S pairs are
merged (~370 bp ASVs) and the longer 18S pairs are concatenated with 10 N.
The taxonomy workflow then:

1. `SPLITASVS`: assigns each ASV to 16S (barrnap bac/arc/mito) or 18S
   (barrnap euk); ASVs without a barrnap hit go to both.
2. `ASSIGNTAXONOMY`: DADA2 `assignTaxonomy` (minBoot 50), 16S against SILVA
   138.2 (to genus), 18S against PR2 5.1.1 (to species), in parallel.
3. `SELECTORGANELLES` + `ASSIGNORGANELLES`: ASVs that SILVA files as
   Chloroplast or Mitochondria (SILVA keeps them under Cyanobacteria /
   Rickettsiales) are classified again against PR2, which contains plastid and
   mitochondrial 16S labelled with the eukaryotic host (`Eukaryota:plas;...`).
4. `MERGETAXONOMY`: writes `asv_taxonomy.tsv` (chosen lineage, lowest rank,
   bootstrap, category prokaryote / plastid / mitochondrion /
   eukaryote_nuclear / unassigned, plus all SILVA and PR2 columns; for
   organelle ASVs also `organelle_host_taxonomy`, `organelle_host_lowest_taxon`,
   `organelle_host_lowest_boot` and the raw `organelle_pr2_*` columns),
   `asv_table_taxonomy.tsv` (counts with taxonomy and host lineage) and
   `taxonomy_overview.tsv`. The main `taxonomy` of organelle ASVs stays SILVA.

References (DADA2 format):

```bash
wget -c https://github.com/pr2database/pr2database/releases/download/v5.1.1/pr2_version_5.1.1_SSU_dada2.fasta.gz
wget -c https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz
```

Run (outputs go to `<ampliseq_results>/taxonomy`):

```bash
nextflow run taxonomy/wfTAXONOMY.nf -c taxonomy/nextflow.config [-profile slurm] \
    --ampliseq_results /path/to/ampliseq/results \
    --silva_db /path/to/silva_nr99_v138.2_toGenus_trainset.fa.gz \
    --pr2_db /path/to/pr2_version_5.1.1_SSU_dada2.fasta.gz \
    [--container_dada2 /path/to/dada2_1.38.sif] [--cpus 16 --memory "32 GB"]
```

On the cluster, adjust the paths at the top of `taxonomy/run_taxonomy_slurm.sh`
(or pass them via `sbatch --export=ALL,AMPLISEQ_RESULTS=...`) and submit it with
`sbatch taxonomy/run_taxonomy_slurm.sh`. It checks all inputs, then runs the
command above with `-profile slurm` as a small controller job.

On the pilot (6,573 ASVs, 8 CPUs) the 16S assignment took ~4 min and the 18S
assignment ~22 min.
