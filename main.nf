#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.bams = ''
params.ref_fasta = ''
params.sequencing_mode = 'srWGS'
params.assembly = 'hg38'
params.callers = ['ExpansionHunter5', 'scattr']

include { path; read_tsv; date_ymd } from './modules/functions'

include { run_expansion_hunter } from './modules/ExpansionHunter.nf'
include { run_scattr } from './modules/scattr.nf'

bams = read_tsv(path(params.bams), ['iid', 'bam'])
ref_fa = path(params.ref_fasta)
ref_fai = path(params.ref_fasta + '.fai')

workflow {

    ref_ch = Channel.value([ref_fa, ref_fai])

    sam_bam_ch =
        Channel.from(bams) |
        map { [it.iid, path(it.bam), path(it.bam + '.bai')] }
    
    if (params.callers.contains('ExpansionHunter5')) {
        run_expansion_hunter(ref, sam_bam_ch)
    }

    if (params.callers.contains('scattr')) {
        run_scattr(ref, sam_bam_ch)
    }

}