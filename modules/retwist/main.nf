process GENERATE_KPOPTWISTED {

    label 'process_high'

    input:
    tuple path(counter_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("*.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file)

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        counter_prefix=\$(echo $counter_file | sed 's/.KPopCounter//')
        $projectDir/bin/KPopCountDB -i \$counter_prefix -s /dev/stdout $args | KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t $prefix -v $args2
        """
}

process KPOPTWIST_UPDATE {

    label 'process_high'

    input:
    tuple path(updating_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("${prefix}.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file), path("${prefix}.KPopTwister"), path("${prefix}.KPopSummary.txt")

    script:
        def args = task.ext.args ?: ''
	def args2 = task.ext.args2 ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updating_prefix=\$(echo $updating_file | sed 's/.KPopTwisted//')
	KPopTwistDB \\
            -i T \$twister_prefix \\
            -i t \$twisted_prefix \\
            -K $args \\
            -s \$updating_prefix ${prefix}
        KPopTwistDB -i T \$twister_prefix -i t \$twisted_prefix -a t \$updating_prefix -o t ${prefix} $args2
        cp $twister_file ${prefix}.KPopTwister 
        """
}

process UPDATE_PLOT {

    label 'process_single'

    input:
    tuple path(updated_twisted_file), val(prefix), path(twister_file), path(twisted_file), path(updated_twister_file), path(kpop_summary_file)
 
    output:
    path("${prefix}_updated_comparison.pdf")

    script:
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updated_twister_prefix=\$(echo $updated_twister_file | sed 's/.KPopTwister//')
        updated_twisted_prefix=\$(echo $updated_twisted_file | sed 's/.KPopTwisted//')
        
        
        $projectDir/bin/KPopPhylo \\
        \$twister_prefix  \\
        \$twisted_prefix \\
        \$twisted_prefix \\
        ${params.kpopPhylo_power} \\
        ${params.kpopPhylo_distance} \\
        ${params.kpopPhylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        
        $projectDir/bin/KPopPhylo \\
        \$updated_twister_prefix  \\
        \$updated_twisted_prefix \\
        \$updated_twisted_prefix \\
        ${params.kpopPhylo_power} \\
        ${params.kpopPhylo_distance} \\
        ${params.kpopPhylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\

        $projectDir/bin/updated_tree.R \${twisted_prefix}.NJ.nwk \${updated_twisted_prefix}.NJ.nwk ${prefix}_updated_comparison.pdf
        """
}
