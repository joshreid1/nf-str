#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.id = ''
params.manifest = ''
params.ref_fasta = ''

include { path; read_tsv; date_ymd } from './modules/functions'
include { index_bam } from './modules/common/sort_bam.nf'
include { minimap2_ubam_illumina; minimap2_ubam_ont; minimap2_ubam_pacbio } from './modules/common/map_ubam.nf'

//ExpansionHunter
include { run_expansion_hunter } from './modules/ExpansionHunter/ExpansionHunter.nf'

//Scatter
include { run_scattr } from './modules/ScatTR/ScatTR.nf'

//Straglr
include { run_straglr as run_straglr } from './modules/Straglr/Straglr.nf'

//LongTR
include { run_longtr } from './modules/LongTR/LongTR.nf'

//STRkit
include { run_strkit } from './modules/STRkit/STRkit.nf'

// atarva
include { run_atarva as run_atarva } from './modules/atarva/atarva.nf'

// TRGT
include { run_trgt } from './modules/TRGT/TRGT.nf'

manifest = read_tsv(path(params.manifest), ['sample', 'type', 'bam', 'align'])

workflow {
    // Read manifest TSV and create channel
    Channel
        .from(manifest)
        .map { record -> 
            def sample = record.sample
            def type = record.type
            def bam_file = file(record.bam) 
            def align = record.align?.toLowerCase()?.trim() == 'yes'
            
            // Validate file exists
            if (!bam_file.exists()) {
                error "BAM file does not exist: ${record.bam} for sample: ${sample}"
            }
            
            // Validate extension
            if (!bam_file.name.endsWith('.bam') && !bam_file.name.endsWith('.cram')) {
                error "Invalid file extension for sample ${sample}"
            }
            
            tuple(sample, type, bam_file, align)
        }
        .branch { 
            to_align: it[3] == true
                return tuple(it[0], it[1], it[2])
            already_aligned: it[3] == false
                return tuple(it[0], it[1], it[2])
        }
        .set { alignment_check }
    
    // Route alignment jobs by sequencing type
    alignment_check.to_align
        .branch { sample, type, bam_file ->
            illumina: 
                type == 'illumina'
                return tuple(sample, type, bam_file)
            ont: 
                type == 'ont'
                return tuple(sample, type, bam_file)
            pacbio: 
                type == 'pacbio'
                return tuple(sample, type, bam_file)
        }
        .set { unaligned_by_type }

    // Align by platform (each returns tuple(sample, type, bam, bai))
    illumina_aligned = minimap2_ubam_illumina(unaligned_by_type.illumina)
    ont_aligned = minimap2_ubam_ont(unaligned_by_type.ont)
    pacbio_aligned = minimap2_ubam_pacbio(unaligned_by_type.pacbio)
    
    // Combine all newly aligned samples (already have index from alignment)
    all_aligned = illumina_aligned
        .mix(ont_aligned, pacbio_aligned)
    
    // Already aligned samples - need to check/generate index
    already_aligned_samples = alignment_check.already_aligned
    
    // Check for index files on all already-aligned samples
    already_aligned_samples
        .map { sample, type, bam_file ->
            def index_file
            def index_exists = false
            
            // Determine expected index file
            if (bam_file.toString().endsWith('.cram')) {
                index_file = file(bam_file.toString() + '.crai')
            } else {
                index_file = file(bam_file.toString() + '.bai')
            }
            
            // Check if index exists
            if (index_file.exists()) {
                index_exists = true
            }
            
            tuple(sample, type, bam_file, index_file, index_exists)
        }
        .branch { sample, type, bam, index, exists ->
            has_index: 
                exists == true
                return tuple(sample, type, bam, index)
            needs_index: 
                exists == false
                return tuple(sample, type, bam)
        }
        .set { indexed_check }
    
    // Generate missing indexes
    generated_indexes = index_bam(indexed_check.needs_index)
    
    // Combine all samples with indexes (both newly aligned and already aligned)
    all_samples_with_index = all_aligned
        .mix(indexed_check.has_index, generated_indexes)
        .branch { sample, type, bam, index ->
            illumina: 
                type == 'illumina'
                return tuple(sample, type, bam, index)
            ont: 
                type == 'ont'
                return tuple(sample, type, bam, index)
            pacbio: 
                type == 'pacbio'
                return tuple(sample, type, bam, index)
        }
        .set { samples }
    
    // Run platform-specific workflows
    illumina_results = run_illumina(samples.illumina)
    ont_results = run_ont(samples.ont)
    pacbio_results = run_pacbio(samples.pacbio)
}

workflow.onComplete {

    println ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """ : """
        Failed: ${workflow.errorReport}
        exit status : ${workflow.exitStatus}
        """
    )
}

workflow run_illumina {
    take:
        sample_ch
    main:
        eh_results = sample_ch  | run_expansion_hunter
        //scattr_results = sample_ch | run_scattr
            
    emit:
        eh_results
        //scattr_results
}

workflow run_ont {
    take:
        sample_ch
    main:
        atarva_results = sample_ch | run_atarva
        straglr_results = sample_ch | run_straglr
        longtr_results = sample_ch | run_longtr
        strkit_results = sample_ch | run_strkit
        
    emit:
        straglr_results
        longtr_results
        atarva_results
        strkit_results
}

workflow run_pacbio {
    take:
        sample_ch
    main:
        atarva_results = sample_ch | run_atarva
        trgt_results = sample_ch | run_trgt
        straglr_results = sample_ch | run_straglr
        longtr_results = sample_ch | run_longtr
        strkit_results = sample_ch | run_strkit
        
    emit:
        straglr_results
        longtr_results
        atarva_results
        strkit_results
}