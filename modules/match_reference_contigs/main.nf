process MATCH_REFERENCE_CONTIGS {
    tag {"Matching with $match_reference"}

    label 'process_low'

    input:
    tuple path(combined_fasta_file), val(prefix), path(match_reference)
 
    output:
    tuple path("*_matched.fasta.gz"), val(prefix)

    script:
        def args = task.ext.args ?: ''
        def args2 = task.ext.args2 ?: ''
        """
        if [[ $combined_fasta_file = *.gz ]]; then
            open_file=zcat
        else 
            open_file=cat
        fi
        lastz ${match_reference}[multiple,unmask] <(\$open_file $combined_fasta_file)[unmask] \\
            --notransition \\
            --step=${params.lastz_step} \\
            --nogapped \\
            --format=general-:name2,size2,strand2,start2,end2,name1,size1,strand1,start1,end1,identity,score \\
            --ambiguous=iupac \\
            $args \\
            | parse_lastz ${params.min_contig_match_len} ${params.min_contig_match_proportion} \\
            | cut -f1 | sort -u > ${prefix}_matching_contigs.txt
        seqkit grep -f ${prefix}_matching_contigs.txt $combined_fasta_file | seqkit seq -w0 -o ${prefix}_matched.fasta.gz
        """
}