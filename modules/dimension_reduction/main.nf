process CLUSTERING {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/dimension_reduction", mode: 'copy'

    input:
    tuple path(retwisted_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    path("*_clusters.txt")

    script:
        def kpopScale_power = task.ext.kpopScale_power ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        retwisted_file_prefix=\$(echo $retwisted_file | sed 's/.KPopTwisted//')
        KPopScale \$twister_prefix \$retwisted_file_prefix \${retwisted_file_prefix} $kpopScale_power
        $projectDir/bin/dimensionality_reduction.R \${retwisted_file_prefix}.embedding.txt ${params.max_dim+1} kmeans ${prefix}_clusters.txt
        """
}