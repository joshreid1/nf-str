# nf-str: Nextflow pipeline for analysis of short tandem repeats (STRs)


Nextflow cohort level STR calling pipeline for short read and long read sequencing. This pipeline is a work in progress.

## Usage (WEHI HPC)

* Clone this repository
* Create and navigate to run directory
* See [bahlolab/nextflow-config](https://github.com/bahlolab/nextflow-config) for generic Nextflow configuration for Milton/SLURM.
* Create configuration file in run directory named `nextflow.config`, e.g.:
  ```Nextflow
    params {
      // inputs
      id = 'str-run'
      manifest = 'bams.tsv'
      
      // run config
      ref_fasta = '/stornext/Bioinf/data/lab_bahlo/ref_db/human/hg38/GATK/fasta_no_alt/hg38.no_alt.fasta'
    }
    ```
* **Params**  
  * `id` - Unique name for run. Used to name output files.
  * `manifest` - Path to TSV file with first column containing sample ID, second column containing sequencing type (i.e. illumina, ont or pacbio), and third column containing path to BAM/CRAM file. No headers or row names. 
  * `callers` - List of STR callers, default is all of them for the releveant sequencing type.
* First run:  
`nextflow run /PATH/TO/nf-str -profile singularity`
* Resume run:  
`nextflow run /PATH/TO/nf-str -resume`


## Output
* Outputs are created in the folder `output` in the run directory

## Implementation
* Calling is implemented as recommended for individual callers using default parameters.

#### Not implemented yet

* For each caller + technology, a directory of standardised parquet files named `<id>.<caller>.<tech>.parquet` is output
* Visualisations
* Allow user to specify parameters for individual callers

