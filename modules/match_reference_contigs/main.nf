process MATCH_REFERENCE_CONTIGS {
    tag {"Matching with $match_reference"}

    label 'process_low'

    input:
    tuple path(combined_fasta_file), val(prefix), path(match_reference)
 
    output:
    tuple path("*.selected.fasta"), path("*.selected.txt"), val(prefix)

    script:
        """
        $projectDir/bin/SelectContigs \\
        $combined_fasta_file \\
        $match_reference \\
        $prefix \\
        ${params.lastz_step} \\
        ${params.min_contig_match_len} \\
        ${params.min_contig_match_proportion} 
        """
}
