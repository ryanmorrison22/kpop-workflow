process PREDICT_TEST_SET {
    
    label 'process_medium'

    input:
    tuple path(train_twister), path(train_twisted), path(test_twisted)
 
    output:
    path("${params.output_prefix}.KPopSummary.txt")
    path("${params.output_prefix}.predictions.txt")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $train_twister | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $train_twisted | sed 's/.KPopTwisted//')
        test_twisted_prefix=\$(echo $test_twisted | sed 's/.KPopTwisted//')
        KPopTwistDB \\
            -i T \$twister_prefix \\
            -i t \$twisted_prefix \\
            -K $args \\
            -s \$test_twisted_prefix ${params.output_prefix} \\
            $args2
        cut -f1,6 ${params.output_prefix}.KPopSummary.txt | sed '1i Sample\tPredicted_Class' > ${params.output_prefix}.predictions.txt
        """    
}