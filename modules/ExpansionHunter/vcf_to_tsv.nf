process vcf_to_tsv {
    
    publishDir "output/eh5/tsvs/", mode: "copy"
    tag { sam }

    input:
    tuple val(sam), path(vcf)
    
    output:
    tuple val(sam), path(out_tsv)
    
    script:
    out_tsv = "${sam}_allele_counts.tsv"
    """
    python prcoess_alleles.py "${vcf}" "${sample}"
    """
}