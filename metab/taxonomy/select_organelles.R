#!/usr/bin/env Rscript
# Pick the organelle 16S ASVs for a second classification against PR2.
# SILVA files plastids under Cyanobacteria (Order Chloroplast) and mitochondria
# under Rickettsiales (Family Mitochondria) without the eukaryotic host; PR2
# contains plastid/mitochondrial 16S labelled with the host lineage (":plas",
# ":mito"). The same rules as the category in merge_taxonomy.R are used.
#
# Usage: select_organelles.R tax_16S_silva.tsv asv_16S.fasta asv_domains.tsv

suppressPackageStartupMessages(library(Biostrings))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: select_organelles.R tax_16S_silva.tsv asv_16S.fasta asv_domains.tsv")

silva <- read.delim(args[1], check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")
seqs <- readDNAStringSet(args[2])
domains <- read.delim(args[3], check.names = FALSE, stringsAsFactors = FALSE, na.strings = "")

organelle <- silva$ASV_ID[silva$Order %in% "Chloroplast" | silva$Family %in% "Mitochondria"]
organelle <- union(organelle, domains$ASV_ID[domains$barrnap_kingdom %in% "mito"])
selected <- seqs[names(seqs) %in% organelle]

writeXStringSet(selected, "asv_organelle.fasta")
cat("organelle ASVs selected for PR2:", length(selected), "of", length(seqs), "16S ASVs\n")
