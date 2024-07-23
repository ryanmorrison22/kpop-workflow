process GENERATE_KPOPTWISTED {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/KPopTwist_files"

    input:
    tuple path(counter_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("temp.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file)

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        KPopCountDB -i $prefix -t updating
        for i in \$(seq 2 \$(awk '{print NF; exit}' updating.KPopCounter.txt)) ; do cut -f1,\${i} updating.KPopCounter.txt ; done | KPopTwistDB -i T \$twister_prefix -k /dev/stdin -o t temp -v
        rm updating.KPopCounter.txt
        """
}

process KPOPTWIST_UPDATE {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/updated_KPopTwist_files"

    input:
    tuple path(updating_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    tuple path("${prefix}.KPopTwisted"), val(prefix), path(twister_file), path(twisted_file)

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updating_prefix=\$(echo $updating_file | sed 's/.KPopTwisted//')
        KPopTwistDB -i T \$twister_prefix -i t \$twisted_prefix -a t \$updating_prefix -o t ${prefix} 
        """
}

process UPDATE_PLOT {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/updated_KPopTwist_files"

    input:
    tuple path(updated_file), val(prefix), path(twister_file), path(twisted_file)
 
    output:
    path("${prefix}_updated_comparison.pdf")

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        updated_prefix=\$(echo $updated_file | sed 's/.KPopTwisted//')
        
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
        \$twister_prefix  \\
        \$updated_prefix \\
        ${params.kpopphylo_power} \\
        ${params.kpopphylo_distance} \\
        ${params.kpopphylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        $args

        $projectDir/bin/updated_tree.R \${twisted_prefix}.NJ.nwk \${updated_prefix}.NJ.nwk  ${prefix}_updated_comparison.pdf
        """
}