params.catalog = "${projectDir}/catalogues/ExpansionHunter_hg38.json"  // path to ExpansionHunter catalog, use 'default' for the default catalog

process call {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "progress/eh5/", mode: "symlink"
    tag { sam }

    input:
        tuple val(sam), val(type), path(bam), path(bai)

    output:
        tuple val(sam), path("${sam}_realigned.bam"), path("${sam}.json"), path("${sam}.vcf")

    script:
    """
        ExpansionHunter \
		--reads ${bam} \
		--reference ${params.ref_fasta} \
		--variant-catalog ${params.catalog} \
		--output-prefix ${sam}
    """
}