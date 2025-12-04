process merge {
    cpus 1
	memory {'2 GB'}
	time '5 h'
    publishDir "output/ehdn/", mode: "copy"

    container 'quay.io/biocontainers/expansionhunterdenovo:0.9.0--h6ac36c1_11'

    input:
		path(manifest)

    output:
        path '*multisample_profile.json'

    script:
    """
        ExpansionHunterDenovo merge \
		--manifest ${manifest} \
		--reference ${params.illumina_ref_fasta} \
		--output-prefix ${params.id}

        # Create a ~100MB file
        dd if=/dev/zero of=step1.txt bs=1M count=100
        echo "Step 1 complete" >> step1.txt
    """
}