process DOWNLOAD_SRAS {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/downloaded_samples", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(list_file)
 
    output:
    path("*/*.sra")

    script:
        def args = task.ext.args ?: ''
        """
        prefetch --option-file $list_file $args
        """
}

process FASTERQ_DUMP {
    tag {"$sra_file"}
    cpus = params.cpu_num
    publishDir "${params.output_dir}/downloaded_samples", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(sra_file)
 
    output:
    tuple eval("echo \$filename"), path("*.fastq.gz")

    script:
        def args = task.ext.args ?: ''
        """
        prefix=`basename $sra_file .sra`
        fasterq-dump \$prefix --threads ${params.cpu_num} $args
        gzip \${prefix}*.fastq
        if [ -f \${prefix}_1.fastq.gz ]; then
            mv \${prefix}_1.fastq.gz \${prefix}_R1.fastq.gz
            mv \${prefix}_2.fastq.gz \${prefix}_R2.fastq.gz
            filename=\${prefix}
        else
            filename=\${prefix}.fastq.gz
        fi
        """
}