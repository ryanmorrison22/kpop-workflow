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
    path("${prefix}_2.3D.pdf")
    path("${prefix}_2.2D.pdf")
    path("${prefix}_2.embedding.txt")
    path("${prefix}_2.NJ.pdf")
    path("${prefix}_2.NJ.nwk"), emit: nwk_file
    path("${prefix}_2.distances.txt")

    script:
        def args = task.ext.args ?: ''
        """
        $projectDir/bin/KPopPhylo $prefix $prefix $args
        """
}