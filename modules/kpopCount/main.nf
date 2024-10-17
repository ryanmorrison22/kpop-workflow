process KPOPCOUNT_BY_CLASS {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/kmer_counts", mode: 'copy'

    input:
    tuple path(fasta_list), val(prefix)
 
    output:
    path("*.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for file in $fasta_list ; do
            class_name=\$(basename \$file _modified.fasta.gz)
            KPopCount -l \$class_name -f <(gzip -c -d \$file) -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$class_name" -L "\$class_name" -N -D -t /dev/stdout 2> /dev/null
        done | KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}

process KPOPCOUNT {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/kmer_counts", mode: 'copy'

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

process KPOPCOUNT_READS {
    cpus = params.cpu_num
    debug true
    publishDir "${params.output_dir}/kmer_counts", mode: 'copy'

    input:
    tuple val(fastq_file_list), val(prefix)
 
    output:
    path("${prefix}.KPopCounter")
    
    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for fileNames in $fastq_file_list ; do
            file=\$(echo \$fileNames | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/,//g')
            num_of_separators=\$(echo \$file | grep -o "?" | wc -l)
            if [[ \$file = *.gz ]]; then
                open_file=zcat
            else 
                open_file=cat
            fi
            if [ \$num_of_separators == 1 ]; then
                R1=\$(echo \$file | cut -d"?" -f1)
                R2=\$(echo \$file | cut -d"?" -f2)
                baseName=\$(basename \$R1 | \\
                    sed 's/\\(.*\\).fastq.*/\\1/' | \\
                    sed 's/\\(.*\\).fq.*/\\1/' | \\
                    sed 's/_R1//g' | \\
                    sed 's/.R1//g' | \\
                    sed 's/-R1//g')
                \$open_file \$R1 \$R2 | KPopCount -l \$baseName -s /dev/stdin -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$baseName" -L "\$baseName" -N -D -t /dev/stdout 2> /dev/null
            elif [ \$num_of_separators == 2 ]; then
                unmated=\$(echo \$file | cut -d"?" -f1)
                R1=\$(echo \$file | cut -d"?" -f2)
                R2=\$(echo \$file | cut -d"?" -f3)
                baseName=\$(basename \$R1 | \\
                    sed 's/\\(.*\\).fastq.*/\\1/' | \\
                    sed 's/\\(.*\\).fq.*/\\1/' | \\
                    sed 's/_R1//g' | \\
                    sed 's/.R1//g' | \\
                    sed 's/-R1//g')
                \$open_file \$R1 \$R2 \$unmated | KPopCount -l \$baseName -s /dev/stdin -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$baseName" -L "\$baseName" -N -D -t /dev/stdout 2> /dev/null
            else
                baseName=\$(basename \$file | \\
                    sed 's/\\(.*\\).fastq.*/\\1/' | \\
                    sed 's/\\(.*\\).fq.*/\\1/')
                \$open_file \$file | KPopCount -l \$baseName -s /dev/stdin -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$baseName" -L "\$baseName" -N -D -t /dev/stdout 2> /dev/null
            fi
        done | KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}

process KPOPCOUNT_READS_BY_CLASS {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/kmer_counts"

    input:
    tuple path(fastq_list), val(prefix)
 
    output:
    path("*.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for file in $fastq_list ; do
            class_name=\$(basename \$file _modified.fastq.gz)
            KPopCount -l \$class_name -s <(gzip -c -d \$file) -k ${params.kmer_len} $args | \\
                KPopCountDB -k /dev/stdin -R "~." -A "\$class_name" -L "\$class_name" -N -D -t /dev/stdout 2> /dev/null
        done | KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}