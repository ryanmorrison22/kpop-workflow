process KPOPCOUNT_BY_CLASS {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/kmer_counts"

    input:
    path(fasta_list)
 
    output:
    path("train.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for file in $fasta_list ; do
            class_name=\$(basename \$file _modified.fasta.gz)
            KPopCount -l \$class_name -f <(gzip -c -d \$file) -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$class_name" -L "\$class_name" -N -D -t /dev/stdout 2> /dev/null
        done | KPopCountDB -k /dev/stdin -o train -v $args2
        """
}

process KPOPCOUNT {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/kmer_counts"

    input:
    tuple path(combined_fasta_file), val(prefix)
 
    output:
    path("${prefix}.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        KPopCount -L -f <(gzip -c -d $combined_fasta_file) -k ${params.kmer_len} $args | \\
            KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}