process CLUSTERING {

    label 'process_low'

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

process KPOPCOUNTER_SAMPLE_REDUCTION {

    label 'process_low'

    input:
    tuple path(input_counter_file), val(prefix), val(num_samples)
 
    output:
    path("*.KPopCounter")

    script:
        def args = task.ext.args ?: ''
        """
        counter_file_prefix=\$(echo $input_counter_file | sed 's/.KPopCounter//')
        reduced_samples=\$(KPopCountDB -i \$counter_file_prefix -R "~." -P $args 2>&1 | \\
            sed 's/Currently selected spectra = \\[ //g' | \\
            sed 's/ \\].//g' | \\
            tr ' ' '\n' | \\
            shuf -n $num_samples | \\
            tr '\n' ',' | \\
            sed "s/'//g" | \\
            sed 's/,\$//')
        KPopCountDB -i \$counter_file_prefix -L \$reduced_samples -N -D -o ${prefix} $args
        """
}