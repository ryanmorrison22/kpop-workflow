process GENERATE_TEST_TWISTED {

    label 'process_high'

    input:
    tuple path(train_twister), path(train_twisted), path(test_fasta_list)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')
        ls $test_fasta_list | \\
        Parallel -t $task.cpus -l 1 -- awk '{filename=\$0; basename=filename; gsub(/\\.(fasta\\.gz|fa\\.gz|fna\\.gz|fasta|fa|fna|_matched|_trimmed)\$/, "", basename);
            if(filename ~ /\\.gz\$/ ) { open_file = "zcat" } else { open_file = "cat" }; 
            system(open_file " " filename " | KPopCount -l " basename " -f /dev/stdin -k "${params.kmer_len}" "$args)}' | \\
        KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t test $args2
        """    
}

process GENERATE_TEST_TWISTED_FROM_READS {

    label 'process_high'

    input:
    tuple path(train_twister), path(train_twisted), val(test_fastq_list)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')

        for fileNames in $test_fastq_list ; do
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
        done | \\
        KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t test $args2
        """    
}

process GENERATE_TEST_TWISTED_FROM_KPOPCOUNTER {
    
    label 'process_high'

    input:
    tuple path(train_twister), path(train_twisted), path(test_counter)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')
        counter_prefix=\$(echo $test_counter | sed 's/.KPopCounter//')
        $projectDir/bin/KPopCountDB -i \$counter_prefix -s /dev/stdout $args | KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t test -v $args2
        """    
}