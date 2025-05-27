process KPOPCOUNT {

    label 'process_high'

    input:
    tuple path(fasta_list), val(prefix)
 
    output:
    path("${prefix}.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        ls $fasta_list | \\
        Parallel -t $task.cpus -l 1 -- awk '{filename=\$0; basename=filename; gsub(/\\.fasta\\.gz\$|\\.fa\\.gz\$|\\.fna\\.gz\$|\\.fasta\$|\\.fa\$|\\.fna\$|\\.selected|_trimmed/, "", basename);
            if(filename ~ /\\.gz\$/ ) { open_file = "zcat" } else { open_file = "cat" }; 
            system(open_file " " filename " | KPopCount -l " basename " -f /dev/stdin -k "${params.kmer_len}" "$args)}' | \\
        KPopCountDB -k /dev/stdin -o $prefix -v $args2
        """
}

process KPOPCOUNT_BY_CLASS {
    tag {"Class: $class_name"}

    label 'process_high'

    input:
    tuple path(fasta_list), val(class_name)
    
    output:
    path("*raw_counts.txt.gz")

    script:
        def args = task.ext.args ?: ''
        """
        ls $fasta_list | \\
        Parallel -t $task.cpus -l 1 -- awk '{filename=\$0; basename=filename; gsub(/\\.fasta\\.gz\$|\\.fa\\.gz\$|\\.fna\\.gz\$|\\.fasta\$|\\.fa\$|\\.fna\$|\\.selected|_trimmed/, "", basename);
            if(filename ~ /\\.gz\$/ ) { open_file = "zcat" } else { open_file = "cat" }; 
            system(open_file " " filename " | KPopCount -l " basename " -f /dev/stdin -k "${params.kmer_len}" "$args)}' | \\
        KPopCountDB -k /dev/stdin -o _train -R "~." -A "$class_name" -L "$class_name" -N -D -t /dev/stdout 2> /dev/null > \\
        ${class_name}_raw_counts.txt
        gzip ${class_name}_raw_counts.txt
        """
}

process KPOPCOUNT_READS {

    label 'process_high'

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
    
    label 'process_high'

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

process KPOPCOUNT_BY_CLASS_FROM_KPOPCOUNTER {
    tag {"Class: $class_name"}
    
    label 'process_high'

    input:
    tuple val(sample_list), val(class_name), path(input_counter_file)
    
    output:
    path("*raw_counts.txt.gz")

    script:
        def args = task.ext.args ?: ''
        """
        counter_file_prefix=\$(echo $input_counter_file | sed 's/.KPopCounter//')
        corrected_list=\$(echo $sample_list | sed 's/\\[//g' | sed 's/\\]//g' | sed 's/ //g')   
        KPopCountDB -i \$counter_file_prefix -L \$corrected_list -A $class_name -L $class_name -N -D -t $class_name 2> /dev/null $args
        mv ${class_name}.KPopCounter.txt ${class_name}_raw_counts.txt
        gzip ${class_name}_raw_counts.txt
        """
}

process KPOPCOUNT_COMBINE_CLASS_COUNTS {
    
    label 'process_high'

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
