params.caller = 'ExpnansionHunter'

include { call } from './ExpansionHunter/call.nf'
include { kmer_filter } from './ExpansionHunter/kmer_filter.nf'
include { sort_bamlet } from './ExpansionHunter/sort_bamlet.nf'

workflow run_expansion_hunter {
    take:
        ref
        sam_bam_ch
    main:
        eh5_results = sam_bam_ch |
            combine(ref) |
            call |
            sort_bamlet |
            kmer_filter 
    emit:
        eh5_results  
}