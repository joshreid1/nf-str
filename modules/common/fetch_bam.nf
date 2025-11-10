process fetch_bam {
    cpus 1
    memory '1 GB'
    time '3 hours'
    tag { sam }
    
    input:
    tuple val(sam), path(url)

    output:
    tuple val(sam), path("${sam}*")
    
    script:
    """
    ext="${url##*.}"
    output="${sam}.${ext}"
    wget -O ${url} ${output}
    """
}

process fetch_bam_s3 {
    cpus 1
    memory '1 GB'
    time '3 hours'
    tag { sam }
    
    container 'amazon/aws-cli'

    input:
    tuple val(sam), path(url)

    output:
    tuple val(sam), path("*")
    
    script:
    """
    amazon/aws-cli s3 --no-sign-request cp ${s3_url} .
    """
}