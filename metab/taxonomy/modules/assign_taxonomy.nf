// DADA2 assignTaxonomy for one ASV set against one reference. In a module so
// that wfTAXONOMY.nf can use it twice (16S/18S, then the organelle ASVs).

process ASSIGNTAXONOMY {

    tag "${label}"
    cpus params.cpus
    memory params.memory
    time "12h"

    publishDir "${params.outdir}", mode: "copy"

    input:
    path r_script
    tuple val(label), path(asv_fasta), path(reference), val(taxlevels)

    output:
    tuple val(label), path("tax_${label}.tsv"), emit: tsv
    path "tax_${label}.args.txt"

    script:
    """
    Rscript "${r_script}" \\
        "${asv_fasta}" "${reference}" "${taxlevels}" ${params.min_boot} ${task.cpus} "tax_${label}"
    """
}
