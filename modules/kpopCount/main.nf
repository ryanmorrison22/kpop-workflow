process KPOPCOUNT {

    label 'process_low'

    input:
    tuple path(fasta_list), val(prefix)
 
    output:
    path("${prefix}.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for file in $fasta_list ; do
            baseName=\$(basename \$file | \\
                sed 's/\\(.*\\).fasta.*/\\1/' | \\
                sed 's/\\(.*\\).fa.*/\\1/' | \\
                sed 's/_matched//g' | \\
                sed 's/_trimmed//g')
            if [[ \$file = *.gz ]]; then
                open_file=zcat
            else
                open_file=cat
            fi
            KPopCount -l \$baseName -f <(\$open_file \$file) -k ${params.kmer_len} $args
        done | KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}

process KPOPCOUNT_BY_CLASS {
    tag {"Class: $class_name"}

    label 'process_low'

    input:
    tuple path(fasta_list), val(class_name)
    
    output:
    path("*raw_counts.txt.gz")

    script:
        def args = task.ext.args ?: ''
        """
        for file in $fasta_list ; do
            baseName=\$(basename \$file | \\
                sed 's/\\(.*\\).fasta.*/\\1/' | \\
                sed 's/\\(.*\\).fa.*/\\1/' | \\
                sed 's/_matched//g' | \\
                sed 's/_trimmed//g')
            if [[ \$file = *.gz ]]; then
                open_file=zcat
            else
                open_file=cat
            fi
            KPopCount -l \$baseName -f <(\$open_file \$file) -k ${params.kmer_len} $args 
        done | \\
        KPopCountDB -k /dev/stdin -R "~." -A "$class_name" -L "$class_name" -N -D -t /dev/stdout 2> /dev/null > \\
        ${class_name}_raw_counts.txt
        gzip ${class_name}_raw_counts.txt
        """
}

process KPOPCOUNT_READS {

    label 'process_low'

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
                KPopCount -l \$baseName -p <(\$open_file \$R1) <(\$open_file \$R2) -k ${params.kmer_len} $args  
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
                KPopCount -l \$baseName -p <(\$open_file \$R1) <(\$open_file \$R2) -s <(\$open_file \$unmated) -k ${params.kmer_len} $args
            else
                baseName=\$(basename \$file | \\
                    sed 's/\\(.*\\).fastq.*/\\1/' | \\
                    sed 's/\\(.*\\).fq.*/\\1/')
                KPopCount -l \$baseName -s <(\$open_file \$file) -k ${params.kmer_len} $args
            fi
        done | KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}

process KPOPCOUNT_READS_BY_CLASS {
    tag {"Class: $class_name"}
    
    label 'process_low'

    input:
    tuple val(fastq_list), val(class_name)
    
    output:
    path("*raw_counts.txt.gz")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        for fileNames in $fastq_list ; do
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
                KPopCount -l \$baseName -p <(\$open_file \$R1) <(\$open_file \$R2) -k ${params.kmer_len} $args  
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
                KPopCount -l \$baseName -p <(\$open_file \$R1) <(\$open_file \$R2) -s <(\$open_file \$unmated) -k ${params.kmer_len} $args
            else
                baseName=\$(basename \$file | \\
                    sed 's/\\(.*\\).fastq.*/\\1/' | \\
                    sed 's/\\(.*\\).fq.*/\\1/')
                KPopCount -l \$baseName -s <(\$open_file \$file) -k ${params.kmer_len} $args
            fi
        done | KPopCountDB -k /dev/stdin -R "~." -A "$class_name" -L "$class_name" -N -D -t /dev/stdout 2> /dev/null $args2 > \\
        ${class_name}_raw_counts.txt
        gzip ${class_name}_raw_counts.txt
        """
}

process KPOPCOUNT_COMBINE_CLASS_COUNTS {
    
    label 'process_low'

    input:
    tuple path(raw_count_list), val(prefix)
    
    output:
    path("*.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        """
        zcat $raw_count_list | KPopCountDB -k /dev/stdin -o $prefix -v $args
        """
}