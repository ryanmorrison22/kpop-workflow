process COMBINE_FILES {
    tag {"$fileName.fileName"}
    cpus = params.cpu_num
    publishDir "${params.output_dir}/modified_fasta_files", mode: 'copy'

    input:
    tuple val(fileName), path(fasta_file)
 
    output:
    path("*modified.fasta.gz")

    script:
        def is_compressed = fasta_file.getName().endsWith(".gz") ? true : false
        """
        baseName=\$(basename $fasta_file | \\
        sed 's/\\(.*\\).fasta.*/\\1/' | \\
        sed 's/\\(.*\\).fa.*/\\1/' | \\
        sed 's/_matched//g')
        echo ">\$baseName" >> \${baseName}_modified.fasta
        if [ "$is_compressed" == "true" ]; then
            gzip -c -d $fasta_file | \\
                grep -v ">" | \\
                sed 's/^/N/g' >> \${baseName}_modified.fasta
        else
            cat $fasta_file | \\
                grep -v ">" | \\
                sed 's/^/N/g' >> \${baseName}_modified.fasta
        fi
        echo 'N' >> \${baseName}_modified.fasta
        gzip \${baseName}_modified.fasta
        """
}

process COMBINE_FASTAS_BY_CLASS {
    tag {"Class $meta_class.meta_class"}
    cpus = params.cpu_num
    publishDir "${params.output_dir}/modified_fasta_files", mode: 'copy'

    input:
    tuple val(meta_class), path(fastas)
 
    output:
    path("*modified.fasta.gz")

    script:
        """
        echo ">class_${meta_class.meta_class}" > class_${meta_class.meta_class}_modified.fasta
        for file in $fastas ; do
            if [[ \$file = *.gz ]]; then
                gzip -c -d \$file | \\
                    grep -v ">" | \\
                    sed 's/^/N/g' >> class_${meta_class.meta_class}_modified.fasta
            else 
                cat \$file | \\
                    grep -v ">" | \\
                    sed 's/^/N/g' >> class_${meta_class.meta_class}_modified.fasta
            fi
        done
        gzip class_${meta_class.meta_class}_modified.fasta
        """
}

process COMBINE_FASTQS_BY_CLASS {
    tag {"Class $meta_class.meta_class"}
    cpus = params.cpu_num
    publishDir "${params.output_dir}/modified_fastq_files", mode: 'copy'

    input:
    tuple val(meta_class), val(fastqs)
 
    output:
    path("*modified.fastq.gz")

    script:
        """
        files=\$(echo $fastqs | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/,//g')
        for fileNames in \$files ; do
            file=\$(echo \$fileNames | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/,//g')
            if [[ \$file = *.gz ]]; then
                gzip -c -d \$file >> class_${meta_class.meta_class}_modified.fastq
            else 
                cat \$file >> class_${meta_class.meta_class}_modified.fastq
            fi
        done
        gzip class_${meta_class.meta_class}_modified.fastq
        """
}