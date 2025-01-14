process KPOPPHYLO {

    label 'process_single'

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
        """
        twister_prefix=\$(echo $twister_file | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo $twisted_file | sed 's/.KPopTwisted//')
        $projectDir/bin/KPopPhylo \\
        \$twister_prefix \\
        \$twisted_prefix \\
        ${prefix} \\
        ${params.kpopPhylo_power} \\
        ${params.kpopPhylo_distance} \\
        ${params.kpopPhylo_magic} \\
        ${params.tree_type} \\
        ${params.tree_label_size} \\
        """
}