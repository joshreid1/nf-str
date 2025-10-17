params.caller = 'ExpnansionHunter'
params.catalog = ''

include { call } from './ExpansionHunter/call.nf'
include { kmer_filter } from './ExpansionHunter/kmer_filter.nf'
include { vcf_to_tsv } from './ExpansionHunter/vcf_to_tsv.nf'
include { concat_tsvs } from './common/concat_tsvs.nf'

workflow run_expansion_hunter {
    take:
        ref
        sam_bam_ch
    main:
        eh5_results = sam_bam_ch |
            combine(ref) |
            call |
            kmer_filter |
            vcf_to_tsv |
            combine(vcf_to_tsv.out_tsv) |
            concat_tsvs(cohort: params.caller)
    emit:
        eh5_results  
}