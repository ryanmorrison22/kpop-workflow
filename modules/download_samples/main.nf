process DOWNLOAD_SRAS {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/downloaded_samples"

    input:
    path(list_file)
 
    output:
    path("*/*.sra")

    script:
        """
        prefetch --option-file $list_file
        """
}

process FASTERQ_DUMP {
    tag {"$sra_file"}
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/downloaded_samples"

    input:
    path(sra_file)
 
    output:
    tuple eval("echo \$filename"), path("*.fastq.gz")

    script:
        """
        prefix=`basename $sra_file .sra`
        fasterq-dump \$prefix
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