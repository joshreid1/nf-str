process concat_tsvs {
    
    publishDir: "output/", mode: "copy"

    input:
    tuple val(cohort), path(tsv_dir)
    
    output:
    val(out_parquet)
    
    script:
    out_tsv = "${cohort}_allele_counts.parquet"
    """
    python prcoess_alleles.py "${vcf}" "${sample}"
    """
}