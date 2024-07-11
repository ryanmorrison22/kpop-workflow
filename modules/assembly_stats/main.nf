process ASSEMBLY_STATS {
    cpus = params.cpu_num
    //conda "${moduleDir}/environment.yml"
    publishDir "${params.output_dir}/assembly_statistics"

    input:
    tuple val(prefix), path(fasta_file)
 
    output:
    path("*assembly_stats.txt")

    script:
        """
        if [[ $fasta_file = *.gz ]]; then
            open_file=zcat
        else 
            open_file=cat
        fi
        num_contigs=\$(\$open_file $fasta_file | grep -c ">")
        seq_len=\$(\$open_file $fasta_file | grep -v ">" | awk '{sum+=length}END{print sum}')
        lar_cont=\$(\$open_file $fasta_file | grep -v ">" | awk '{print length}' | sort -gr | head -1)
        sma_cont=\$(\$open_file $fasta_file | grep -v ">" | awk '{print length}' | sort -g | head -1)
        quants=\$(\$open_file $fasta_file | grep -v ">" | awk '{print length}' | sort -g | \\
        awk '{a[NR]=\$0}END{quart25=a[int(((NR/4)*1)+0.9999)] ; \\
        quart75=a[int(((NR/4)*3)+0.9999)] ; \\
        median=a[int(((NR/4)*2)+0.9999)] ; \\
        print "upperQuartileContig\t"quart75 ; print "medianContig\t"median ; print "lowerQuartileContig\t"quart25 }')

        echo "fileName\t${prefix.fileName}" >> ${prefix.fileName}_assembly_stats.txt
        echo "numContigs\t\$num_contigs" >> ${prefix.fileName}_assembly_stats.txt
        echo "totalLength\t\$seq_len" >> ${prefix.fileName}_assembly_stats.txt
        echo "largestContig\t\$lar_cont" >> ${prefix.fileName}_assembly_stats.txt
        echo "smallestContig\t\$sma_cont" >> ${prefix.fileName}_assembly_stats.txt
        echo "\$quants" >> ${prefix.fileName}_assembly_stats.txt
        """
}

//printf("%d\n",((NR/4)*1)+=((NR/4)*1)<0?0:0.999)