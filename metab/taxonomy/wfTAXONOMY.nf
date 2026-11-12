// Taxonomy for mixed 16S/18S ampliseq output (515F-Y/926R amplicons).
// ampliseq runs with skip_taxonomy; this workflow splits its ASVs by barrnap's
// best hit and classifies 16S ASVs against SILVA and 18S ASVs against PR2
// (DADA2 assignTaxonomy), then merges both into one table. ASVs that SILVA
// calls plastid or mitochondrion are classified again against PR2, which has
// organelle 16S labelled with the eukaryotic host lineage.

params.ampliseq_results = params.ampliseq_results ?: ""
params.asv_fasta = params.asv_fasta ?: "${params.ampliseq_results}/dada2/ASV_seqs.fasta"
params.asv_table = params.asv_table ?: "${params.ampliseq_results}/dada2/ASV_table.tsv"
params.barrnap_summary = params.barrnap_summary ?: "${params.ampliseq_results}/barrnap/summary.tsv"
params.outdir = params.outdir ?: "${params.ampliseq_results}/taxonomy"
params.script_dir = params.script_dir ?: projectDir.toString()
params.silva_db = params.silva_db ?: ""
params.pr2_db = params.pr2_db ?: ""
params.silva_taxlevels = params.silva_taxlevels ?: "Kingdom,Phylum,Class,Order,Family,Genus"
params.pr2_taxlevels = params.pr2_taxlevels ?: "Domain,Supergroup,Division,Subdivision,Class,Order,Family,Genus,Species"
params.min_boot = params.min_boot ?: 50
params.cpus = params.cpus ?: 8
params.memory = params.memory ?: "16 GB"

// Modules only see params that are set before the include.
include { ASSIGNTAXONOMY } from './modules/assign_taxonomy.nf'
include { ASSIGNTAXONOMY as ASSIGNORGANELLES } from './modules/assign_taxonomy.nf'


process SPLITASVS {

    cpus 1
    memory "4GB"
    time "1h"

    publishDir "${params.outdir}", mode: "copy", pattern: "asv_domains.tsv"

    input:
    path r_script
    path asv_fasta
    path barrnap_summary

    output:
    path "asv_16S.fasta", emit: fasta_16S
    path "asv_18S.fasta", emit: fasta_18S
    path "asv_domains.tsv", emit: domains

    script:
    """
    Rscript "${r_script}" "${asv_fasta}" "${barrnap_summary}"
    """
}


process SELECTORGANELLES {

    cpus 1
    memory "4GB"
    time "1h"

    input:
    path r_script
    path tax_16S
    path asv_16S
    path domains

    output:
    path "asv_organelle.fasta", emit: fasta

    script:
    """
    Rscript "${r_script}" "${tax_16S}" "${asv_16S}" "${domains}"
    """
}


process MERGETAXONOMY {

    cpus 1
    memory "8GB"
    time "1h"

    publishDir "${params.outdir}", mode: "copy"

    input:
    path r_script
    path domains
    path tax_16S
    path tax_18S
    path tax_organelle
    path asv_table
    path asv_fasta

    output:
    path "asv_taxonomy.tsv"
    path "asv_table_taxonomy.tsv"
    path "taxonomy_overview.tsv"

    script:
    """
    Rscript "${r_script}" \\
        "${domains}" "${tax_16S}" "${tax_18S}" "${tax_organelle}" "${asv_table}" "${asv_fasta}" ${params.min_boot}
    """
}


workflow {

    if (!params.silva_db || !params.pr2_db) {
        error "Set --silva_db and --pr2_db to the DADA2-formatted reference files."
    }

    asv_fasta = file(params.asv_fasta, checkIfExists: true)
    asv_table = file(params.asv_table, checkIfExists: true)
    barrnap_summary = file(params.barrnap_summary, checkIfExists: true)
    silva_db = file(params.silva_db, checkIfExists: true)
    pr2_db = file(params.pr2_db, checkIfExists: true)

    script_dir = params.script_dir
    SPLITASVS(file("${script_dir}/split_asvs.R"), asv_fasta, barrnap_summary)

    to_assign = SPLITASVS.out.fasta_16S
        .map { fasta -> tuple("16S_silva", fasta, silva_db, params.silva_taxlevels) }
        .mix(SPLITASVS.out.fasta_18S
            .map { fasta -> tuple("18S_pr2", fasta, pr2_db, params.pr2_taxlevels) })

    ASSIGNTAXONOMY(file("${script_dir}/assign_taxonomy.R"), to_assign)

    tax = ASSIGNTAXONOMY.out.tsv
        .toSortedList { a, b -> a[0] <=> b[0] }
        .map { results -> results.collectEntries { label, tsv -> [label, tsv] } }

    tax_16S = tax.map { it["16S_silva"] }

    SELECTORGANELLES(
        file("${script_dir}/select_organelles.R"),
        tax_16S,
        SPLITASVS.out.fasta_16S,
        SPLITASVS.out.domains
    )
    ASSIGNORGANELLES(
        file("${script_dir}/assign_taxonomy.R"),
        SELECTORGANELLES.out.fasta
            .map { fasta -> tuple("organelle_pr2", fasta, pr2_db, params.pr2_taxlevels) }
    )

    MERGETAXONOMY(
        file("${script_dir}/merge_taxonomy.R"),
        SPLITASVS.out.domains,
        tax_16S,
        tax.map { it["18S_pr2"] },
        ASSIGNORGANELLES.out.tsv.map { label, tsv -> tsv },
        asv_table,
        asv_fasta
    )
}
