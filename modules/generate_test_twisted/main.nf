
process GENERATE_TEST_TWISTED {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/predictions", mode: 'copy'

    input:
    tuple path(train_twister), path(train_twisted), path(test_fasta)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')
        KPopCount -L -f <(gzip -c -d $test_fasta) -k ${params.kmer_len} $args | \\
            KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t test $args2
        """    
}