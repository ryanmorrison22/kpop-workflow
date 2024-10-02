process GENERATE_KPOPTWISTED {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/KPopTwist_files"

    input:
    tuple path(counter_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("*.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file)

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        counter_prefix=\$(echo $counter_file | sed 's/.KPopCounter//')
        KPopCountDB -i \$counter_prefix -t updating
        for i in \$(seq 2 \$((\$(awk '{print NF; exit}' updating.KPopCounter.txt)+1))) ; do cut -f1,\${i} updating.KPopCounter.txt ; done | KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t $prefix -v
        rm updating.KPopCounter.txt
        """
}

process KPOPTWIST_UPDATE {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/updated_KPopTwist_files"

    input:
    tuple path(updating_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("${prefix}.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file), path("${prefix}.KPopTwister")

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updating_prefix=\$(echo $updating_file | sed 's/.KPopTwisted//')
        KPopTwistDB -i T \$twister_prefix -i t \$twisted_prefix -a t \$updating_prefix -o t ${prefix} 
        cp $twister_file ${prefix}.KPopTwister 
        """
}

process UPDATE_PLOT {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/updated_KPopTwist_files"

    input:
    tuple path(updated_twisted_file), val(prefix), path(twister_file), path(twisted_file), path(updated_twister_file)
 
    output:
    path("${prefix}_updated_comparison.pdf")

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updated_twister_prefix=\$(echo $updated_twister_file | sed 's/.KPopTwister//')
        updated_twisted_prefix=\$(echo $updated_twisted_file | sed 's/.KPopTwisted//')
        
        
        $projectDir/bin/KPopPhylo \\
        \$twister_prefix  \\
        \$twisted_prefix \\
        ${params.kpopphylo_power} \\
        ${params.kpopphylo_distance} \\
        ${params.kpopphylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        $args
        
        $projectDir/bin/KPopPhylo \\
        \$updated_twister_prefix  \\
        \$updated_twisted_prefix \\
        ${params.kpopphylo_power} \\
        ${params.kpopphylo_distance} \\
        ${params.kpopphylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        $args

        $projectDir/bin/updated_tree.R \${twisted_prefix}.NJ.nwk \${updated_twisted_prefix}.NJ.nwk  ${prefix}_updated_comparison.pdf
        """
}