
process GENERATE_TEST_TWISTED {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/predictions"

    input:
    tuple path(train_twister), path(train_twisted), path(test_fasta)
 
    output:
    tuple path(train_twister), path(train_twisted), path("test.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        KPopCount -L -f <(gzip -c -d $test_fasta) -k ${params.kmer_len} $args | \\
            KPopTwistDB -i T train -k /dev/stdin -o t test $args2
        """    
}