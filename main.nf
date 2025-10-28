#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.manifest = ''
params.ref_fasta = ''

include { path; read_tsv; date_ymd } from './modules/functions'

//ExpansionHunter
include { run_expansion_hunter } from './modules/ExpansionHunter/ExpansionHunter.nf'

//Scatter
include { run_scattr } from './modules/scattr.nf'

//Straglr
include { run_straglr as run_straglr_ont } from './modules/Straglr/Straglr.nf'
include { run_straglr as run_straglr_pacbio } from './modules/Straglr/Straglr.nf'

//LongTR
include { run_longtr as run_longtr_ont } from './modules/LongTR/LongTR.nf'
include { run_longtr as run_longtr_pacbio } from './modules/LongTR/LongTR.nf'

manifest = read_tsv(path(params.manifest), ['sample', 'type', 'bam'])

workflow {

    // Read manifest TSV and create channel
    Channel
        .from(manifest)
        .map { record -> 
            def sample = record.sample
            def type = record.type
            def bam_file = file(record.bam)
            def index_file
            
            // Determine index file extension
            if (bam_file.toString().endsWith('.cram')) {
                index_file = file(record.bam + '.crai')
            } else {
                index_file = file(record.bam + '.bai')
            }
            [sample, type, bam_file, index_file]
        }
        .branch { sample, type, bam, index ->
            illumina: type == 'illumina'
                return [sample, bam, index]
            ont: type == 'ont'
                return [sample, bam, index]
            pacbio: type == 'pacbio'
                return [sample, bam, index]
        }
        .set { samples }
    
     // Run platform-specific workflows
    illumina_results = run_illumina(samples.illumina)
    ont_results = run_ont(samples.ont)
    pacbio_results = run_pacbio(samples.pacbio)
    
    // Combine all results
    //all_results = illumina_results.mix(ont_results, pacbio_results)
}

workflow run_illumina {
    take:
        sample_ch
    main:
        eh_results = sample_ch  | run_expansion_hunter
            
    emit:
        eh_results
}

workflow run_ont {
    take:
        sample_ch
    main:
        straglr_results = sample_ch | run_straglr_ont
        longtr_results = sample_ch | run_longtr_ont
        
    emit:
        straglr_results
        longtr_results
}

workflow run_pacbio {
    take:
        sample_ch
    main:
        straglr_results = sample_ch | run_straglr_pacbio
        longtr_results = sample_ch | run_longtr_pacbio
        
    emit:
        straglr_results
        longtr_results
}