process META_COLOURED_TREE {

    label 'process_single'

    input:
    tuple path(meta_file), path(nwk_file), val(prefix)
 
    output:
    path("${prefix}_meta_tree.pdf")

    script:
        """
        $projectDir/bin/meta_coloured_tree.R $meta_file $nwk_file ${prefix}_meta_tree.pdf
        """
}