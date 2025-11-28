process download_s3_files {
    tag "${sample}_${type}"

    publishDir "s3_files/${sample}/${type}", mode: 'copy'

	//module 'awscli'
	container 'quay.io/biocontainers/awscli:1.8.3--0'

    cpus 4
    memory '4 GB'
    time '12h'

    input:
    tuple val(sample), val(type), val(s3_uri), val(align)

    output:
    tuple val(sample), val(type), path("${filename}"), val(align)

    script:
    filename = s3_uri.tokenize('/')[-1]
    """
    # Configure AWS CLI for optimized parallel transfers
    aws configure set default.s3.max_concurrent_requests 50
    aws configure set default.s3.max_bandwidth 500MB/s
    aws configure set default.s3.multipart_threshold 64MB
    aws configure set default.s3.multipart_chunksize 16MB

    # Download file from S3
    aws s3 cp --no-sign-request ${s3_uri} ${filename}

    # Verify download
    if [ ! -f "${filename}" ]; then
        echo "Error: Failed to download ${filename}" >&2
        exit 1
    fi

    echo "Successfully downloaded ${filename} (\$(du -h ${filename} | cut -f1))"
    """
}
