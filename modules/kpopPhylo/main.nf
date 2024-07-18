process KPOPPHYLO {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/trees_and_metrics"

    input:
    tuple path(twister_file), path(twisted_file), val(prefix)
 
    output:
    path("${prefix}.KPopTwisted.txt")
    path("${prefix}.sigm.KPopMetrics.txt")
    path("${prefix}.powx.KPopMetrics.txt")
    path("${prefix}.KPopInertia.txt")
    path("${prefix}.3D.pdf")
    path("${prefix}.2D.pdf")
    path("${prefix}.embedding.txt")
    path("${prefix}.NJ.pdf")
    path("${prefix}.NJ.nwk"), emit: nwk_file
    path("${prefix}.distances.txt")

    script:
        def args = task.ext.args ?: ''
        """
        $projectDir/bin/KPopPhylo \\
        $prefix \\
        $prefix \\
        ${params.kpopphylo_power} \\
        ${params.kpopphylo_distance} \\
        ${params.kpopphylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        $args
        """
}