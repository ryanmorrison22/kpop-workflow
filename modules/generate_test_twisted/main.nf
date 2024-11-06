
process GENERATE_TEST_TWISTED {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/predictions", mode: 'copy'

    input:
    tuple path(train_twister), path(train_twisted), path(test_fasta_list)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')
        for file in $test_fasta_list ; do
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
        done | KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t test $args2
        """    
}