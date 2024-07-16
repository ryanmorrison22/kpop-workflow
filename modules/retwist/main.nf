process GENERATE_KPOPTWISTED {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/KPopTwist_files"

    input:
    tuple path(counter_file), val(prefix)
 
    output:
    tuple path("temp.KPopTwisted"), val(prefix)

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo ${params.twister_file} | sed 's/.KPopTwister//')
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
    tuple path(updating_file), val(prefix)
 
    output:
    path("${prefix}.KPopTwisted")

    script:
        def args = task.ext.args ?: ''
        """
        twister_prefix=\$(echo ${params.twister_file} | sed 's/.KPopTwister//')
        twisted_prefix=\$(echo ${params.twisted_file} | sed 's/.KPopTwisted//')
        updating_prefix=\$(echo $updating_file | sed 's/.KPopTwisted//')
        KPopTwistDB -i T \$twister_prefix -i t \$twisted_prefix -a t \$updating_prefix -o t ${prefix} 
        """
}

