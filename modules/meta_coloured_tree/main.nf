process META_COLOURED_TREE {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/trees_and_metrics"

    input:
    tuple path(meta_file), path(nwk_file), val(prefix)
 
    output:
    path("${prefix}_meta_tree.pdf")

    script:
        """
        $projectDir/bin/meta_coloured_tree.R $meta_file $nwk_file ${prefix}_meta_tree.pdf
        """
}