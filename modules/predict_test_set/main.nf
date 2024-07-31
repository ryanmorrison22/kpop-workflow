process PREDICT_TEST_SET {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/predictions"

    input:
    tuple path(train_twister), path(train_twisted), path(test_twisted)
 
    output:
    path("${params.output_prefix}.KPopSummary.txt")
    path("${params.output_prefix}.predictions.txt")

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        KPopTwistDB \\
            -i T train \\
            -i t train \\
            -K $args \\
            -s test ${params.output_prefix} \\
            $args2
        cut -f1,6 ${params.output_prefix}.KPopSummary.txt | sed '1i Sample\tPredicted_Class' > ${params.output_prefix}.predictions.txt
        """    
}