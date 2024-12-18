process KPOPTWIST {
    cpus = params.cpu_num
    publishDir "${params.output_dir}/KPopTwist_files", mode: 'copy'

    input:
    tuple path(counter_file), val(prefix)
 
    output:
    tuple path("${prefix}.KPopTwister"), path("${prefix}.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        """
        counter_prefix=\$(echo $counter_file | sed 's/.KPopCounter//')
        KPopTwist -i \$counter_prefix -o $prefix -v $args
        """
}