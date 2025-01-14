process INPUT_VALIDATION {
    tag {fast_file}

    label 'process_single'

    input:
    tuple val(fileName), path(fast_file)
 
    output:
    tuple val(fileName), path(fast_file)

    script:
        def is_paired = fast_file.size() == 2 ? true : false 
        def is_three = fast_file.size() == 3 ? true : false
        """
        if [ "$is_paired" == "true" ]; then
            R1=\$(echo $fast_file | cut -d" " -f1)
            R2=\$(echo $fast_file | cut -d" " -f2)
            files=(\$R1 \$R2)
        elif [ "$is_three" == "true" ]; then
            unmated=\$(echo $fast_file | cut -d" " -f1)
            R1=\$(echo $fast_file | cut -d" " -f2)
            R2=\$(echo $fast_file | cut -d" " -f3)
            files=(\$R1 \$R2 \$unmated)
        else
            files=($fast_file)
        fi

        for file in \${files[@]}; do
            val_checkFile=\$(seqkit seq \$file -t dna -s -v  | grep [ACGTacgt] | wc -l)
            if [[ \$val_checkFile == 0 ]]; then 
                exit 1
            fi
        done
        """
}