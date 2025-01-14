process FASTQ_QC {
    tag {"$fileName.fileName"}

    label 'process_low'

    input:
    tuple val(fileName), path(fastq_file)
 
    output:
    tuple val(fileName), path("*.fq.gz")

    script:
        def args = task.ext.args ?: ''
        def is_paired = fastq_file.size() == 2 ? true : false 
        def is_three = fastq_file.size() == 3 ? true : false
        """
        if [ "$is_paired" == "true" ]; then
            R1=\$(echo $fastq_file | cut -d" " -f1)
            R2=\$(echo $fastq_file | cut -d" " -f2)
            trim_galore --fastqc --paired \$R1 \$R2 --cores $task.cpus --gzip $args
        elif [ "$is_three" == "true" ]; then
            unmated=\$(echo $fastq_file | cut -d" " -f1)
            R1=\$(echo $fastq_file | cut -d" " -f2)
            R2=\$(echo $fastq_file | cut -d" " -f3)
            trim_galore --fastqc --paired \$R1 \$R2 --cores $task.cpus --gzip $args
            trim_galore --fastqc \$unmated --cores $task.cpus --gzip $args
        else
            trim_galore --fastqc $fastq_file --cores $task.cpus --gzip $args
        fi
        """
}