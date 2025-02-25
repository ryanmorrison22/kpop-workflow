#!/usr/bin/env nextflow

// Call processes from modules directory 
include { KPOPCOUNT as KPOPCOUNT1 } from './modules/kpopCount'
include { KPOPCOUNT as KPOPCOUNT2 } from './modules/kpopCount'
include { KPOPCOUNT_BY_CLASS as KPOPCOUNT_BY_CLASS1 } from './modules/kpopCount'
include { KPOPCOUNT_BY_CLASS as KPOPCOUNT_BY_CLASS2 } from './modules/kpopCount'
include { KPOPCOUNT_READS_BY_CLASS as KPOPCOUNT_READS_BY_CLASS1 } from './modules/kpopCount'
include { KPOPCOUNT_READS_BY_CLASS as KPOPCOUNT_READS_BY_CLASS2 } from './modules/kpopCount'
include { KPOPCOUNT_COMBINE_CLASS_COUNTS as KPOPCOUNT_COMBINE_CLASS_COUNTS1 } from './modules/kpopCount'
include { KPOPCOUNT_COMBINE_CLASS_COUNTS as KPOPCOUNT_COMBINE_CLASS_COUNTS2 } from './modules/kpopCount'
include { KPOPCOUNT_READS as KPOPCOUNT_READS1 } from './modules/kpopCount'
include { KPOPCOUNT_READS as KPOPCOUNT_READS2 } from './modules/kpopCount'
include { KPOPCOUNT_BY_CLASS_FROM_KPOPCOUNTER } from './modules/kpopCount'
include { GENERATE_TEST_TWISTED } from './modules/generate_test_twisted'
include { GENERATE_TEST_TWISTED_FROM_READS } from './modules/generate_test_twisted'
include { GENERATE_TEST_TWISTED_FROM_KPOPCOUNTER } from './modules/generate_test_twisted'
include { KPOPPHYLO } from './modules/kpopPhylo'
include { KPOPTWIST as KPOPTWIST1 } from './modules/kpopTwist'
include { KPOPTWIST as KPOPTWIST2 } from './modules/kpopTwist'
include { PREDICT_TEST_SET } from './modules/predict_test_set'
include { ASSEMBLE_FASTQS as ASSEMBLE_FASTQS1 } from './modules/assemble_fastqs'
include { ASSEMBLE_FASTQS as ASSEMBLE_FASTQS2 } from './modules/assemble_fastqs'
include { MATCH_REFERENCE_CONTIGS as MATCH_REFERENCE_CONTIGS1 } from './modules/match_reference_contigs'
include { MATCH_REFERENCE_CONTIGS as MATCH_REFERENCE_CONTIGS2 } from './modules/match_reference_contigs'
include { ASSEMBLY_STATS as ASSEMBLY_STATS1 } from './modules/assembly_stats'
include { ASSEMBLY_STATS as ASSEMBLY_STATS2 } from './modules/assembly_stats'
include { META_COLOURED_TREE } from './modules/meta_coloured_tree'
include { DOWNLOAD_SRAS as DOWNLOAD_SRAS1 } from './modules/download_samples'
include { DOWNLOAD_SRAS as DOWNLOAD_SRAS2 } from './modules/download_samples'
include { FASTERQ_DUMP as FASTERQ_DUMP1 } from './modules/download_samples'
include { FASTERQ_DUMP as FASTERQ_DUMP2 } from './modules/download_samples'
include { KPOPTWIST_UPDATE; UPDATE_PLOT } from './modules/retwist'
include { GENERATE_KPOPTWISTED as GENERATE_KPOPTWISTED1 } from './modules/retwist'
include { GENERATE_KPOPTWISTED as GENERATE_KPOPTWISTED2 } from './modules/retwist'
include { INPUT_VALIDATION as FASTA_VALIDATION } from './modules/input_validation'
include { INPUT_VALIDATION as TEST_FASTA_VALIDATION } from './modules/input_validation'
include { INPUT_VALIDATION as FASTQ_VALIDATION } from './modules/input_validation'
include { INPUT_VALIDATION as TEST_FASTQ_VALIDATION } from './modules/input_validation'
include { CLUSTERING } from './modules/dimension_reduction'
include { KPOPCOUNTER_SAMPLE_REDUCTION } from './modules/dimension_reduction'
include { FASTQ_QC as FASTQ_QC1 } from './modules/quality_control'
include { FASTQ_QC as FASTQ_QC2 } from './modules/quality_control'

// Create a help message
def helpMessage() {
        log.info"""
Usage:
        nextflow run main.nf [--cluster|--classify] --input_dir <path/to/dir> [--test_dir <path/to/dir> (required only if --classify used)]

Required arguments:
    Workflow (At least one of these; if one of these is not supplied the program will download, QC and assemble files then quit)
        --cluster                           Data is run through clustering workflow, starting with an unknown dataset the pipeline produces \
                                            a distance matrix and pseudophylogenetic tree showing relatedness between samples. Creates \
                                            .KPopTwister and KPopTwisted files               
        --classify                          Data is run through classification workflow, starting with separate training and test datasets, \
                                            a model is created using the training dataset and known class metadata> This model is used to \
                                            predict the classes of the unknown test dataset. Requires --test_dir argument. Creates \
                                            .KPopTwister and KPopTwisted files
        --update                            New data is run added to an existing database, creating updated .KPopTwister and KPopTwisted files. \
                                            The data type should match the original data type that the database was built on. \
                                            In other words, avoid updating a database built from fasta files with fastqs and vice versa

    Input (At least one of these)                  
        --input_dir                         Path to directory containing fasta/fastq files. Paired-end fastqs require "R1" and "R2" in filenames. \
                                            Gzipped files are allowed. \
                                            If --classify used this directory is the training dataset. \
                                            If --update used this directory is the new dataset used to update
        --accession_list                    Supply a list of SRA IDs to download as input samples in the form of a text file, with one SRA per line
        --test_dir                          Directory containing unseen test dataset. Only required if --classify workflow invoked 
        --test_accession_list               Supply a list of SRA IDs to download as test samples in the form of a text file, with one SRA per line. \
                                            Only required if --classify workflow invoked
        --twisted_file                      Full path to .KPopTwisted file. Only required for --update workflow
        --twister_file                      Full path to .KPopTwister file. Only required for --update workflow
        --kpopcount_input                   .KPopCounter file as an input. \
                                            Incompatible with --no_assembly, --match_reference, --input_dir, --accession_list and --meta_data options
        --kpopcount_test_input              .KPopCounter file as a test input. Only required if --classify workflow invoked. \
                                            Incompatible with --no_assembly, --match_reference, --input_dir, --accession_list, --max_dim and --meta_data options 

Optional arguments:
    General arguments
        --kmer_len                          Length of k-mer to use when generating spectra. For choosing best k size for DNA we recommend using the base 4 log of the genome size. \
                                            E.g. If the genome size is 5Mbp, then base 4 log of 5000000 is 11.13, so optimal kmer size is likely 11 or 12 [12] 
        --output_dir                        Path to output directory. If directory doesn't exist then a new directory will be created. [projectDir/results]
        --output_prefix                     Prefix for output files [output]
        --no_assembly                       Do not perform assembly on the reads, the workflow will count the number of kmers from the raw reads directly instead of assemblies
        --no_qc                             Do not perform quality control using trim_galore   
        --validate_inputs                   Perform validation check on fastq and fasta inputs to ensure they are formatted correctly, incorrectly formatted files will be skipped                                           
        --meta_data                         Tsv file with two required columns with defined headers; "fileName" and "class". \
                                            "fileName" is file name if a fasta or fasta.gz file, or file prefix if paired-end fastqs. E.g. sample1.fasta.gz if fasta file or \
                                            sample1 if sample1_R1.fastq.gz and sample1_R2.fastq.gz. Additional columns allowed
        --match_reference                   Full path to reference fasta file, used to select contigs that only match the supplied reference
        --max_dim                           Maximum number of dimensions used to separate data. Choosing 0 uses all available dimensions, which will be one less than the number of samples \
                                            for --cluster or one less than the number of classes if --classify. A lower number will reduce memory usage. \
                                            If the data cannot be separated into the number chosen, less dimensions will be chosen automatically. Must not be a number above the maximum number of samples [0]
        --min_contig_match_len              Minimum number of query contig base pairs that match the reference. Only used with --match_reference option [250]
        --min_contig_match_proportion       Minimum fraction of query contig base pairs that match reference. Only used with --match_reference option [0.6]
        --pred_class_num                    Specify the top n number of best predictions to be included in .KPopSummary file. E.g. 2 would choose the top two closest classes [all]
        --tree_type                         Specify the type of tree generated by ggtree - 'rectangular' or 'circular' ['rectangular']
        --tree_label_size                   Specify the size of the labels on the tree generated by ggtree, choose 0 to remove labels [3]
        --max_cpus                          Specify the maximum number of CPUs used during the running of the workflow [maximum available up to 64 CPUs]
        --max_memory                        Specify the maximum number of GBs used during the running of the workflow [maximum available up to 256GB]
        --max_hours                         Specify the maximum number of hours used during a process when running of the workflow [maximum available up to 64 hours]
        --help                              Print help instructions

    Nextflow arguments
        -profile conda                      Install the required conda environment automatically from the environment.yml file found in the same directory as main.nf. Slower than installing it manually.

    Flash arguments (https://pubmed.ncbi.nlm.nih.gov/21903629/)
        --flash_minOverlap                  The minimum required overlap length between two reads to provide a confident overlap. Only used on fastq inputs [20]
        --flash_maxOverlap                  Maximum overlap length expected in approximately 90% of read pairs. Only used on fastq inputs [1000]
        --extra_flash                       Any additional arguments for flash. E.g. --extra_flash '-O -x 0.35'

    Megahit arguments (https://pubmed.ncbi.nlm.nih.gov/25609793/)
        --extra_megahit                     Any additional arguments for Megahit. E.g. --extra_megahit '--k-min 25'

    TrimGalore arguments (https://github.com/FelixKrueger/TrimGalore)
        --extra_trimGalore                  Any additional arguments for TrimGalore. E.g. --extra_trimGalore '--quality 40'

    SRA-toolkit arguments (https://github.com/ncbi/sra-tools)
        --extra_prefetch                    Any additional arguments for fasterq-dump. E.g. --extra_prefetch '--verify'
        --extra_fasterq_dump                Any additional arguments for fasterq-dump. E.g. --extra_fasterq_dump '--concatenate-reads'

    KPop arguments (https://www.biorxiv.org/content/10.1101/2022.06.22.497172v2)
        --extra_kpopCount                   Any additional arguments for KPopCount. E.g. --extra_kpopCount
        --extra_kpopCountDB                 Any additional arguments for KPopCountDB. E.g. --extra_kpopCountDB
        --extra_kpopTwist                   Any additional arguments for KPopTwist. E.g. --extra_kpopTwist
        --extra_kpopTwistDB                 Any additional arguments for KPopTwistDB. E.g. --extra_kpopTwistDB
        --kpopPhylo_power                   Set the external power when computing distances [2]
        --kpopPhylo_distance                Distance measure to be used. This must be one of 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'. ['euclidean']
        --kpopPhylo_magic                   Cluster-related variable (Not currently implemented) ['1.'] 
        --kpopScale_power                   Set the external power when computing distances [2]
         """
         .stripIndent()
}

// Print important information for each run at the start
log.info"""
${workflow.manifest.name} v${workflow.manifest.version}
         ==========================
         input from   : ${params.input_dir}
         output to    : ${params.output_dir}
         kmer length  : ${params.kmer_len}
         --
         run as       : ${workflow.commandLine}
         max cpus     : ${params.max_cpus}
         max memory   : ${params.max_memory}
         max hours    : ${params.max_hours}
         started at   : ${workflow.start}
         config file  : ${workflow.configFiles}
         """
         .stripIndent()

// Print help message if --help invoked
if (params.help){
    helpMessage()
    exit 0
}

// Define warnings
if (!params.cluster && !params.classify && !params.update){
    log.warn"""Current run will only perform steps before KPop analysis. --cluster, --classify or --update required for KPop workflows. Use --help for more information.
    """.stripIndent()
}

if (params.classify && params.meta_data == "" && params.kpopcount_input == ""){
    log.warn"""Each individual file in --input_dir will be treated as a class as --meta_data not specified. If you want to combine files into classes or rename classes please use --meta_data. Use --help for more information.
    """.stripIndent()
}

// Define errors
if (params.input_dir == "" && params.accession_list == "" && params.kpopcount_input == ""){
    log.error"""--input_dir or --accession_list or --kpopcount_input required for workflow. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.classify && params.test_dir == "" && params.test_accession_list == "" && params.kpopcount_test_input == ""){
    log.error"""--test_dir, --test_accession_list or kpopcount_test_input required for --classify workflow. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.kpopcount_input != "" && (params.no_assembly == true || params.match_reference != "" || params.input_dir != "" || params.accession_list != "" || params.meta_data != "")){
    log.error"""--kpopcount_input is incompatible with --no_assembly, --match_reference, --input_dir, --accession_list, and --meta_data options. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.no_assembly == true && params.match_reference != "" ){
    log.error"""--no_assembly is incompatible with --match_reference option as this workflow does not perform fastq mapping. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.input_dir != "") { // If an input directory is supplied
    
    // Channels created for all fastas, paired fastas and all paired fastqs in input_dir
    FASTAS = "${params.input_dir}/*.{fasta,fa,fna,fasta.gz,fa.gz,fna.gz}"
    FASTQS = "${params.input_dir}/*.{fastq,fq,fastq.gz,fq.gz}"
    PAIRED_FASTQS = "${params.input_dir}/*_{R1,R2}.{fastq,fq,fastq.gz,fq.gz}"

    Channel
        .fromPath( FASTAS )
        .map { it -> [[fileName: it.toString().split("/")[-1]], file(it)]}
        .set {fasta_files}

    Channel
        .fromFilePairs( PAIRED_FASTQS )
        .map { it -> [[fileName: it[0].toString()], it[1]]}
        .set {paired_fastq_files}

    Channel
        .fromPath( FASTQS )
        .map { it -> [[fileName: it.toString().split("/")[-1]], file(it)]}
        .filter { !(it.fileName =~ /(R1|R2)/) }
        .concat(paired_fastq_files)
        .set {fastq_files}

    if (!file( FASTAS ) && !file( FASTQS ) && !file( PAIRED_FASTQS)){
        log.error"""--input_dir needs to contain at least one fasta or fastq file. Please check specified directory again. Use --help for more information.
        """.stripIndent()
        exit 0
    }

    fasta_files
        .concat(fastq_files)
        .map {it -> it[0].fileName}
        .filter{ if (it =~/[\[\]\?\*\ ]/) {throw new IllegalArgumentException("Value '$it' contains one of these illegal characters: '[', ']', '?', '*', ' '")} else { it }}

} else {
    // Create empty channels if no fasta or fastq files found in input_dir, this allows for downloading of data from SRA
    Channel
        .empty()
        .set {fasta_files}
    
    Channel
        .empty()
        .set {paired_fastq_files}

    Channel
        .empty()
        .set {fastq_files}
}

// Create channel that includes class meta data for training
if (params.meta_data != "") {
    if (!("fileName" in file(params.meta_data).splitCsv(header:true, sep: "\t").first().keySet()) || !("class" in file(params.meta_data).splitCsv(header:true, sep: "\t").first().keySet())) {
        log.error"""File given with --meta_data needs to have two columns named 'fileName' and 'class' in it. Please check header names in specified file again. Use --help for more information.
        """.stripIndent()
        exit 0  
    }

    Channel
        .fromPath(params.meta_data)
        .splitCsv(header:true, sep: "\t")
        .map { row -> meta = [[fileName: row.fileName.toString().split("/")[-1]], [meta_class: row.class]] }
        .set {meta_file}
    
    if (params.no_assembly) {
        fastq_files
            .join(meta_file)
            .map {it -> [it[2], it[1]]}
            .groupTuple(by: [0])
            .set {meta_fastq_files}
    } else {
        fasta_files
            .join(meta_file)
            .map {it -> [it[2], it[1]]}
            .groupTuple(by: [0])
            .set {meta_fasta_files}
    }
} 

// Channels created for all fastas, paired fastas and all paired fastqs in test_dir
if (params.test_dir != "") {
    TEST_FASTAS = "${params.test_dir}/*.{fasta,fa,fna,fasta.gz,fa.gz,fna.gz}"
    TEST_FASTQS = "${params.test_dir}/*.{fastq,fq,fastq.gz,fq.gz}"
    TEST_PAIRED_FASTQS = "${params.test_dir}/*_{R1,R2}.{fastq,fq,fastq.gz,fq.gz}"
    
    if (!file( TEST_FASTAS ) && !file( TEST_FASTQS ) && !file( TEST_PAIRED_FASTQS ) && params.test_accession_list == ""){
        log.error"""--test_dir needs to contain at least one fasta or fastq file. Please check specified directory again. Use --help for more information.
        """.stripIndent()
        exit 0
    }

    Channel
    .fromPath( TEST_FASTAS )
    .map { it -> [[fileName: it.toString().split("/")[-1]], file(it)]}
    .set {test_fasta_files}

    Channel
    .fromFilePairs( TEST_PAIRED_FASTQS )
    .map { it -> [[fileName: it[0].toString()], it[1]]}
    .set {test_paired_fastq_files}

    Channel
    .fromPath( TEST_FASTQS )
    .map { it -> [[fileName: it.toString().split("/")[-1]], file(it)]}
    .filter { !(it.fileName =~ /(R1|R2)/) }
    .concat(test_paired_fastq_files)
    .set {test_fastq_files}

    test_fasta_files
    .concat(test_fastq_files)
    .map {it -> it[0].fileName}
    .filter{ if (it =~/[\[\]\?\*\ ]/) {throw new IllegalArgumentException("Value '$it' contains one of these illegal characters: '[', ']', '?', '*', ' '")} }

} else {

    // Create empty channels if no fasta or fastq files found in input_dir, this allows for downloading of data from SRA
    Channel
        .empty()
        .set {test_fasta_files}
    
    Channel
        .empty()
        .set {test_paired_fastq_files}

    Channel
        .empty()
        .set {test_fastq_files}
}

workflow {
    /// Download training set
    if (params.accession_list != "") {
        DOWNLOAD_SRAS1(file(params.accession_list))
            .flatten()
            .set {SRA_FILES}
        FASTERQ_DUMP1(SRA_FILES)
            .map {it -> [[fileName: it[0].toString()], it[1]]}
            .concat(fastq_files)
            .set {fastq_files}
    }

    /// Download test set
    if (params.test_accession_list != "") {
        DOWNLOAD_SRAS2(file(params.test_accession_list))
            .flatten()
            .set {TEST_SRA_FILES}
        FASTERQ_DUMP2(TEST_SRA_FILES)
            .map {it -> [[fileName: it[0].toString()], it[1]]}
            .concat(test_fastq_files)
            .set {test_fastq_files}
    }

    /// Input Validation 
    if (params.validate_inputs) {
        FASTQ_VALIDATION(fastq_files)
            .set {fastq_files}
        TEST_FASTQ_VALIDATION(test_fastq_files)
            .set {test_fastq_files}
        FASTA_VALIDATION(fasta_files)
            .set {fasta_files}
        TEST_FASTA_VALIDATION(test_fasta_files)
            .set {test_fasta_files}
    }

    if (!params.no_qc) {
        /// FASTQ Quality Control
        FASTQ_QC1(fastq_files)
            .set { fastq_files }
        FASTQ_QC2(test_fastq_files)
            .set { test_fastq_files }
    }

    /// Either assemble the reads (default) or use fastqs if '--no_assembly' is false
    if (params.no_assembly) {
        assembled_fastas = fastq_files
    } else {
        ASSEMBLE_FASTQS1(fastq_files)
            .map(it -> [it[0], it[1]])
            .set {assembled_fastas}
        ASSEMBLY_STATS1(assembled_fastas)
    }

    /// Clustering workflow
    if (params.cluster) {
        if (params.kpopcount_input != "") { // If .KPopCount is supplied as input
            Channel
                .fromPath( params.kpopcount_input )
                .map(it -> [it, params.output_prefix])
                .set {kpopcount_file}
        } else if (params.no_assembly) { // If starting from reads

            //TO DO: USE READS THAT MAP ONLY TO A REFERENCE

            fastq_files
                .map(it -> it[1].toString().replace(", ", "?")) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                .toSortedList()
                .map(it -> [it, params.output_prefix])
                .set {fastq_list}
            KPOPCOUNT_READS1(fastq_list)
                .map(it -> [it, params.output_prefix])
                .set {kpopcount_file}
        
        } else { // If starting with fastas
            assembled_fastas.concat(fasta_files)
                .set {concat_fasta_files}

            // Only use sequences that match with reference if reference provided
            if (params.match_reference != "") {
                concat_fasta_files
                    .map(it -> [it[1], it[0].fileName.toString().split("/")[-1]
                    .replace(".fasta.gz", "").replace(".fastq.gz", "")
                    .replace(".fasta", "").replace(".fastq", "")
                    .replace(".fna", "").replace(".fna.gz", "")
                    .replace(".fa.gz", "").replace(".fq.gz", "")
                    .replace(".fa", "").replace(".fq", ""), file(params.match_reference, checkIfExists: true)])
                    .set {adjusted_concat_fasta_files}
                MATCH_REFERENCE_CONTIGS1(adjusted_concat_fasta_files)
                    .map(it -> [[fileName: it[1]], it[0]])
                    .set {concat_fasta_files}
                }
            concat_fasta_files
                .map(it -> it[1])
                .toSortedList()
                .map(it -> [it, params.output_prefix])
                .set {fasta_list}
            KPOPCOUNT1(fasta_list)
                .map(it -> [it, params.output_prefix])
                .set {kpopcount_file}
        }

        if (params.max_dim != 0) { // If the dimension size is specified, else run all dimensions available

            // Reduce dimensions from KPopCounts
            if (params.kpopcount_input != "") {
                kpopcount_file
                    | map(it -> [it[0], "reduced_dim", params.max_dim+1])
                    | KPOPCOUNTER_SAMPLE_REDUCTION
                    | map(it -> [it, "reduced_dim"])
                    | KPOPTWIST1
                    | set {red_dim_kpoptwist_files}
            
                // Twisting all samples with reduced twister
                GENERATE_KPOPTWISTED1(kpopcount_file.combine(red_dim_kpoptwist_files))
                    .map(it -> [it[0], "reduced_temp", it[2], it[3]])
                    .set {retwisted_files}

                /// Create embeddings using KPopScale
                CLUSTERING(retwisted_files)
                    .set {cluster_file}

                /// KPopCount per cluster
                cluster_file
                    .splitCsv(header:true, sep: "\t")
                    .map { row -> meta = [[fileName: row.fileName.toString().split("/")[-1]], [new_class: "clusteredClass_${row.class}"]] }
                    .set {cluster_meta_file}
                
                cluster_meta_file
                    .map(it -> [it[1], it[0].fileName])
                    .groupTuple(by: [0])
                    .map(it -> [it[1], it[0].new_class, file(params.kpopcount_input)])
                    .set {cluster_kpopcounter_list}
                
                KPOPCOUNT_BY_CLASS_FROM_KPOPCOUNTER(cluster_kpopcounter_list)
                    .collect()
                    .map(it -> [it, "clustered"])
                    .set {raw_cluster_count_list}

                KPOPCOUNT_COMBINE_CLASS_COUNTS1(raw_cluster_count_list)
                    .map(it -> [it, "clustered"])
                    .set {cluster_kpopcount_file}
            }

            // Generate KPopCounts and KPopTwister for the number of dimensions
            else if (params.no_assembly) { // If starting from reads
                fastq_files
                    | randomSample( params.max_dim+1, 1) 
                    | map(it -> it[1].toString().replace(", ", "?")) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                    | toSortedList
                    | map(it -> [it, "reduced_dim"])
                    | KPOPCOUNT_READS2
                    | map(it -> [it, "reduced_dim"])
                    | KPOPTWIST1
                    | set {red_dim_kpoptwist_files}

                // Twisting all samples with reduced twister
                GENERATE_KPOPTWISTED1(kpopcount_file.combine(red_dim_kpoptwist_files))
                    .map(it -> [it[0], "reduced_temp", it[2], it[3]])
                    .set {retwisted_files}

                /// Create embeddings using KPopScale
                CLUSTERING(retwisted_files)
                    .set {cluster_file}

                // KPopCount per cluster
                cluster_file
                    .splitCsv(header:true, sep: "\t")
                    .map { row -> meta = [[fileName: row.fileName.toString().split("/")[-1]], [new_class: "clusteredClass_${row.class}"]] }
                    .set {cluster_meta_file}

                fastq_files
                    .map( it -> [[fileName: it[0].fileName.toString()
                        .replace(".fastq.gz", "")
                        .replace(".fastq", "")
                        .replace(".fq.gz", "")
                        .replace(".fq", "")
                        .replace("_matched", "")], it[1].toString().replace(", ", "?")])  //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS_BY_CLASS
                    .join(cluster_meta_file, by: [0])
                    .map {it -> [it[2], it[1]]}
                    .groupTuple(by: [0])
                    .map {it -> [it[1], it[0].new_class]}
                    .set {cluster_fastq_list}

                KPOPCOUNT_READS_BY_CLASS1(cluster_fastq_list)
                    .collect()
                    .map(it -> [it, "clustered"])
                    .set {raw_fastq_cluster_count_list}

                KPOPCOUNT_COMBINE_CLASS_COUNTS1(raw_fastq_cluster_count_list)
                    .map(it -> [it, "clustered"])
                    .set {cluster_kpopcount_file}

            } else { 
                concat_fasta_files
                    | randomSample( params.max_dim+1, 1) 
                    | map(it -> it[1])
                    | toSortedList()
                    | map(it -> [it, "reduced_dim"])
                    | KPOPCOUNT2
                    | map(it -> [it, "reduced_dim"])
                    | KPOPTWIST1
                    | set {red_dim_kpoptwist_files}
            
                // Twisting all samples with reduced twister
                GENERATE_KPOPTWISTED1(kpopcount_file.combine(red_dim_kpoptwist_files))
                    .map(it -> [it[0], "reduced_temp", it[2], it[3]])
                    .set {retwisted_files}

                /// Create embeddings using KPopScale
                CLUSTERING(retwisted_files)
                    .set {cluster_file}

                // KPopCount per cluster
                cluster_file
                    .splitCsv(header:true, sep: "\t")
                    .map { row -> meta = [[fileName: row.fileName.toString().split("/")[-1]], [new_class: "clusteredClass_${row.class}"]] }
                    .set {cluster_meta_file}
                
                concat_fasta_files
                    .map( it -> [[fileName: it[1].toString().split("/")[-1]
                        .replace(".fasta.gz", "")
                        .replace(".fasta", "")
                        .replace(".fa.gz", "")
                        .replace(".fa", "")
                        .replace(".fna", "")
                        .replace(".fna.gz", "")
                        .replace("_matched", "")], it[1]])
                    .join(cluster_meta_file, by: [0])
                    .map {it -> [it[2], it[1]]}
                    .groupTuple(by: [0])
                    .map {it -> [it[1], it[0].new_class]}
                    .set {cluster_fasta_list}
                
                KPOPCOUNT_BY_CLASS1(cluster_fasta_list)
                    .collect()
                    .map(it -> [it, "clustered"])
                    .set {raw_cluster_count_list}

                KPOPCOUNT_COMBINE_CLASS_COUNTS1(raw_cluster_count_list)
                    .map(it -> [it, "clustered"])
                    .set {cluster_kpopcount_file}
            }

            // Generate Twister file based on clusters
            KPOPTWIST2(cluster_kpopcount_file)
                .set {clustered_kpoptwist_files}

            // Twist all samples using cluster twister
            GENERATE_KPOPTWISTED2(kpopcount_file.combine(clustered_kpoptwist_files))
                .map(it -> [it[2], it[0], params.output_prefix])
                .set {kpoptwist_files}

        } else { // Use all samples and dimensions instead
            KPOPTWIST1(kpopcount_file)
                .map(it -> [it[0], it[1], params.output_prefix])
                .set {kpoptwist_files}
        }

        // Generate tree and additional output files produced by KPopPhylo
        KPOPPHYLO(kpoptwist_files).nwk_file
            .set {phylo_nwk_file}

        // Generate tree with leaves coloured by metadata 
        if (params.meta_data != "") {
            phylo_nwk_file
                .map(it -> [file(params.meta_data), it, params.output_prefix])
                .set {meta_ch}
            META_COLOURED_TREE(meta_ch)
        }
    }

    /// Classification workflow
    if (params.classify) {

        //TO DO: USE READS THAT MAP ONLY TO A REFERENCE

        // Linking data to supplied metadata file if available, if not then each sample is treated as own class
        if (params.meta_data != "") {
            if (params.no_assembly) {
                fastq_files
                    .join(meta_file)
                    .map {it -> [it[2], it[1]]}
                    .concat(meta_fastq_files.transpose())
                    .ifEmpty { 
                        log.error "There are no matches in the meta file to the file names found in the input directory. Please remember to include the extensions for non-paired read inputs." 
                        exit 1}
                    .set {concat_meta_fastq_files}
            } else {
                assembled_fastas
                    .join(meta_file)
                    .map {it -> [it[2], it[1]]}
                    .concat(meta_fasta_files.transpose())
                    .ifEmpty { 
                        log.error "There are no matches in the meta file to the file names found in the input directory. Please remember to include the extensions for non-paired read inputs." 
                        exit 1}
                    .set {concat_meta_fasta_files}
            }
        } else {
            if (params.no_assembly) {
                fastq_files
                    .map {it -> [[meta_class: it[0].fileName], it[1]]}
                    .set {concat_meta_fastq_files}
            } else {
                assembled_fastas
                    .concat(fasta_files)
                    .map {it -> [[meta_class: it[0].fileName], it[1]]}
                    .set {concat_meta_fasta_files}
            }
        }

        // Only use train sequences that match with reference if reference provided
        if (params.match_reference != "") {
            concat_meta_fasta_files
                .map(it -> [it[1], "${it[0].meta_class}_${it[1].getBaseName(file(it[1]).name.endsWith('.gz')? 2: 1)}", file(params.match_reference, checkIfExists: true)])
                .set {adjusted_concat_meta_fasta_files}
            MATCH_REFERENCE_CONTIGS1(adjusted_concat_meta_fasta_files)
                .map(it -> [[meta_class: it[1].toString().split("_")[0]], it[0]])
                .set {concat_meta_fasta_files}
        }

        // Assemble test samples, calculate assembly stats and combine all files 
        if (params.no_assembly) {
            test_fastq_files
                .map(it -> it[1].toString().replace(", ", "?")) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                .toSortedList()
                .map(it -> [it, "test"])
                .set {test_fastq_list}

        } else {
            ASSEMBLE_FASTQS2(test_fastq_files)
                .set {test_assembled_fastas}
            ASSEMBLY_STATS2(test_assembled_fastas)
            test_assembled_fastas.concat(test_fasta_files)
                .set {test_concat_fasta_files}  
            // Only use test sequences that match with reference if reference provided
            if (params.match_reference != "") {
                test_concat_fasta_files
                    .map(it -> [it[1], it[0].fileName.toString().split("/")[-1]
                    .replace(".fasta.gz", "").replace(".fastq.gz", "")
                    .replace(".fasta", "").replace(".fastq", "")
                    .replace(".fa.gz", "").replace(".fq.gz", "")
                    .replace(".fna", "").replace(".fna.gz", "")
                    .replace(".fa", "").replace(".fq", ""), file(params.match_reference, checkIfExists: true)])
                    .set {test_adjusted_concat_fasta_files}
                MATCH_REFERENCE_CONTIGS2(test_adjusted_concat_fasta_files)
                    .map(it -> [[fileName: it[1]], it[0]])
                    .set {test_concat_fasta_files}
            }
            
            test_concat_fasta_files
                .map(it -> it[1])
                .toSortedList()
                .map(it -> [it, "test"])
                .set {test_fasta_list}
        }    

        // Combine and count kmers training data
        if (params.kpopcount_input != "") { // If .KPopCount is supplied as input
            Channel
                .fromPath( params.kpopcount_input )
                .map(it -> [it, "train"])
                .set {train_kpopcount_file}
        } else if (params.no_assembly) { // If starting from reads
            concat_meta_fastq_files
                .map(it -> [it[0], it[1].toString().replace(", ", "?")]) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                .groupTuple(by: [0])
                .map(it -> [it[1], it[0].meta_class])
                .set {train_fastq_list}

            KPOPCOUNT_READS_BY_CLASS1(train_fastq_list)
                .collect()
                .map(it -> [it, "train"])
                .set {raw_class_count_list}

            KPOPCOUNT_COMBINE_CLASS_COUNTS1(raw_class_count_list)
                .map(it -> [it, "train"])
                .set {train_kpopcount_file}
        } else {
            concat_meta_fasta_files
                .groupTuple(by: [0])
                .map(it -> [it[1], it[0].meta_class])
                .set {train_fasta_list}

            KPOPCOUNT_BY_CLASS1(train_fasta_list)
                .collect()
                .map(it -> [it, "train"])
                .set {raw_class_count_list}

            KPOPCOUNT_COMBINE_CLASS_COUNTS1(raw_class_count_list)
                .map(it -> [it, "train"])
                .set {train_kpopcount_file}
        }

        if (params.max_dim != 0) { // If the dimension size is specified, else run all dimensions available   

            if (params.no_assembly) { // If starting from reads
                concat_meta_fastq_files
                    | randomSample( params.max_dim+1, 1) 
                    | map(it -> it[1].toString().replace(", ", "?")) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                    | toSortedList
                    | map(it -> [it, "reduced_dim"])
                    | KPOPCOUNT_READS2
                    | map(it -> [it, "reduced_dim"])
                    | KPOPTWIST1
                    | set {train_red_dim_kpoptwist_files}

                // Twisting all samples with reduced twister
                GENERATE_KPOPTWISTED1(train_kpopcount_file.combine(train_red_dim_kpoptwist_files))
                    .map(it -> [it[0], "reduced_temp", it[2], it[3]])
                    .set {train_retwisted_files}

                /// Create embeddings using KPopScale
                CLUSTERING(train_retwisted_files)
                    .set {train_cluster_file}
                
                // KPopCount per cluster
                train_cluster_file
                    .splitCsv(header:true, sep: "\t")
                    .map { row -> meta = [[meta_class: row.fileName], [new_class: "clusteredClass_${row.class}"]] }
                    .set {train_cluster_meta_file}

                concat_meta_fastq_files
                    .map( it -> [it[0], it[1].toString().replace(", ", "?")])  //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS_BY_CLASS
                    .groupTuple(by: [0])
                    .join(train_cluster_meta_file, by: [0][0])
                    .map {it -> [it[1], it[2]]}
                    .groupTuple(by: [1])
                    .map {it -> [it[0].flatten().toList(), it[1].new_class]}
                    .set {train_cluster_fastq_list}

                KPOPCOUNT_READS_BY_CLASS2(train_cluster_fastq_list)
                    .collect()
                    .map(it -> [it, "clustered"])
                    .set {raw_fastq_cluster_count_list}

                KPOPCOUNT_COMBINE_CLASS_COUNTS2(raw_fastq_cluster_count_list)
                    .map(it -> [it, "clustered"])
                    .set {train_cluster_kpopcount_file}


            } else {

                concat_meta_fasta_files
                    | randomSample( params.max_dim+1, 1) 
                    | map(it -> it[1])
                    | toSortedList()
                    | map(it -> [it, "reduced_dim"])
                    | KPOPCOUNT2
                    | map(it -> [it, "reduced_dim"])
                    | KPOPTWIST1
                    | set {train_red_dim_kpoptwist_files}
            
                // Twisting all samples with reduced twister
                GENERATE_KPOPTWISTED1(train_kpopcount_file.combine(train_red_dim_kpoptwist_files))
                    .map(it -> [it[0], "reduced_temp", it[2], it[3]])
                    .set {train_retwisted_files}

                /// Create embeddings using KPopScale
                CLUSTERING(train_retwisted_files)
                    .set {train_cluster_file}

                // KPopCount per cluster
                train_cluster_file
                    .splitCsv(header:true, sep: "\t")
                    .map { row -> meta = [[meta_class: row.fileName], [new_class: "clusteredClass_${row.class}"]] }
                    .set {train_cluster_meta_file}

                concat_meta_fasta_files
                    .groupTuple(by: [0])
                    .join(train_cluster_meta_file, by: [0][0])
                    .map {it -> [it[1], it[2]]}
                    .groupTuple(by: [1])
                    .map {it -> [it[0].flatten().toList(), it[1].new_class]}
                    .set {train_cluster_fasta_list}

                KPOPCOUNT_BY_CLASS2(train_cluster_fasta_list)
                    .collect()
                    .map(it -> [it, "clustered"])
                    .set {raw_cluster_count_list}

                KPOPCOUNT_COMBINE_CLASS_COUNTS2(raw_cluster_count_list)
                    .map(it -> [it, "clustered"])
                    .set {train_cluster_kpopcount_file}
            }

            // Generate Twister file based on clusters
            KPOPTWIST2(train_cluster_kpopcount_file)
                .set {train_clustered_kpoptwist_files}

            // Twist all samples using cluster twister
            GENERATE_KPOPTWISTED2(train_kpopcount_file.combine(train_clustered_kpoptwist_files))
                .map(it -> [[it[2], it[0]], "train"])
                .set {train_kpoptwist_files}

        } else {
            // Twist training data
            KPOPTWIST1(train_kpopcount_file)
                .map(it -> [it, "train"])
                .set {train_kpoptwist_files}
        }
        
        // Twist test files based on training twister and use the positional information within multidimensional space to predict classes

        if (params.kpopcount_input != "") { // If .KPopCount is supplied as input
            train_test_files = train_kpoptwist_files
                .map(it -> [it[0][0], it[0][1], file(params.kpopcount_test_input)])
            train_test_kpop_files = GENERATE_TEST_TWISTED_FROM_KPOPCOUNTER(train_test_files)
        } else if (params.no_assembly) { // If starting from reads
            train_test_files = train_kpoptwist_files
                .combine(test_fastq_list)
                .map(it -> [it[0][0], it[0][1], it[2]])
            train_test_kpop_files = GENERATE_TEST_TWISTED_FROM_READS(train_test_files)
        } else {
            train_test_files = train_kpoptwist_files
                .combine(test_fasta_list)
                .map(it -> [it[0][0], it[0][1], it[2]])
            train_test_kpop_files = GENERATE_TEST_TWISTED(train_test_files)
        }
        PREDICT_TEST_SET(train_test_kpop_files)
    }

    /// Update workflow

    if (params.update) {
        if (params.no_assembly) { // If starting from reads

            //TO DO: USE READS THAT MAP ONLY TO A REFERENCE

            fastq_files
                .map(it -> it[1].toString().replace(", ", "?")) //Need to replace the ", " with something unique so we can separate the files in KPOPCOUNT_READS 
                .toSortedList()
                .map(it -> [it, params.output_prefix])
                .set {fastq_list}
            KPOPCOUNT_READS1(fastq_list)
                .map(it -> [it, "temp", file(params.twister_file), file(params.twisted_file)])
                .set {kpopcount_file}
        
        } else {
            
            assembled_fastas.concat(fasta_files)
                .set {concat_fasta_files}
            
            if (params.match_reference != "") {
                concat_fasta_files
                    .map(it -> [it[1], it[0].fileName, file(params.match_reference, checkIfExists: true)])
                    .set {adjusted_concat_fasta_files}
                MATCH_REFERENCE_CONTIGS1(adjusted_concat_fasta_files)
                    .map(it -> [[fileName: it[1]], it[0]])
                    .set {concat_fasta_files}
            }

            concat_fasta_files
                .map(it -> it[1])
                .toSortedList()
                .map(it -> [it, params.output_prefix])
                .set {fasta_list}
            KPOPCOUNT1(fasta_list)
                .map(it -> [it, "temp", file(params.twister_file), file(params.twisted_file)])
                .set {kpopcount_file}
        }
        GENERATE_KPOPTWISTED1(kpopcount_file)
            .map(it -> [it[0], params.output_prefix, it[2], it[3]])
            .set {updating_file}
        KPOPTWIST_UPDATE(updating_file)
            .set {updated_files}
        UPDATE_PLOT(updated_files)
    }
}
