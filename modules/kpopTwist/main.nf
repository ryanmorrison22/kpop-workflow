process KPOPTWIST {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/KPopTwist_files"

    input:
    tuple path(counter_file), val(prefix)
 
    output:
    tuple path("${prefix}.KPopTwister"), path("${prefix}.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        """
        KPopTwist -i $prefix -o $prefix -v $args
        """
}