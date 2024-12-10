process ASSEMBLE_FASTQS {
    tag {"$fileName.fileName"}
    cpus = params.cpu_num
    publishDir "${params.output_dir}/assemblies", mode: 'copy'

    input:
    tuple val(fileName), path(fastq_file)
 
    output:
    tuple val(fileName), path("*.fasta.gz")

    script:
        def is_paired = fastq_file.size() == 2 ? true : false 
        def is_three = fastq_file.size() == 3 ? true : false
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        if [ "$is_paired" == "true" ]; then
            baseName=$fileName.fileName
            R1=\$(echo $fastq_file | cut -d" " -f1)
            R2=\$(echo $fastq_file | cut -d" " -f2)
            
            flash \\
                -m ${params.flash_minOverlap} \\
                -M ${params.flash_maxOverlap} \\
                \$R1 \\
                \$R2 \\
                -o \$baseName \\
                -t ${params.cpu_num} \\
                $args

            megahit \\
                -r \${baseName}.extendedFrags.fastq \\
                -1 \${baseName}.notCombined_1.fastq \\
                -2 \${baseName}.notCombined_2.fastq \\
                -o \$baseName \\
                -t ${params.cpu_num} \\
                $args2
        elif [ "$is_three" == "true" ]; then
            baseName=$fileName.fileName
            unmated=\$(echo $fastq_file | cut -d" " -f3)
            R1=\$(echo $fastq_file | cut -d" " -f1)
            R2=\$(echo $fastq_file | cut -d" " -f2)
            
            flash \\
                -m ${params.flash_minOverlap} \\
                -M ${params.flash_maxOverlap} \\
                \$R1 \\
                \$R2 \\
                -o \$baseName \\
                -t ${params.cpu_num} \\
                $args

            megahit \\
                -r \${baseName}.extendedFrags.fastq,\$unmated \\
                -1 \${baseName}.notCombined_1.fastq \\
                -2 \${baseName}.notCombined_2.fastq \\
                -o \$baseName \\
                -t ${params.cpu_num} \\
                $args2
        else
            baseName=\$(basename $fastq_file | \\
                sed 's/\\(.*\\).fastq.*/\\1/' | \\
                sed 's/\\(.*\\).fq.*/\\1/')
            
            megahit \\
                -r $fastq_file \\
                -o \$baseName \\
                -t ${params.cpu_num} \\
                $args2
        fi
        mv \${baseName}/final.contigs.fa \${baseName}.fasta
        gzip \${baseName}.fasta
        """    
}