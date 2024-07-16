#!/usr/bin/env nextflow

// Call processes from modules directory 
include { COMBINE_FILES; COMBINE_FILES_BY_CLASS } from './modules/combine_files'
include { KPOPCOUNT; KPOPCOUNT_BY_CLASS } from './modules/kpopCount'
include { GENERATE_TEST_TWISTED } from './modules/generate_test_twisted'
include { KPOPPHYLO } from './modules/kpopPhylo'
include { KPOPTWIST } from './modules/kpopTwist'
include { PREDICT_TEST_SET } from './modules/predict_test_set'
include { ASSEMBLE_FASTQS as ASSEMBLE_FASTQS1 } from './modules/assemble_fastqs'
include { ASSEMBLE_FASTQS as ASSEMBLE_FASTQS2 } from './modules/assemble_fastqs'
include { MATCH_REFERENCE_CONTIGS as MATCH_REFERENCE_CONTIGS1} from './modules/match_reference_contigs'
include { MATCH_REFERENCE_CONTIGS as MATCH_REFERENCE_CONTIGS2} from './modules/match_reference_contigs'
include { ASSEMBLY_STATS as ASSEMBLY_STATS1} from './modules/assembly_stats'
include { ASSEMBLY_STATS as ASSEMBLY_STATS2} from './modules/assembly_stats'
include { META_COLOURED_TREE } from './modules/meta_coloured_tree'
include { DOWNLOAD_SRAS as DOWNLOAD_SRAS1} from './modules/download_samples'
include { DOWNLOAD_SRAS as DOWNLOAD_SRAS2} from './modules/download_samples'
include { FASTERQ_DUMP as FASTERQ_DUMP1} from './modules/download_samples'
include { FASTERQ_DUMP as FASTERQ_DUMP2} from './modules/download_samples'
include { GENERATE_KPOPTWISTED; KPOPTWIST_UPDATE } from './modules/retwist'

// Create a help message
def helpMessage() {
        log.info"""
Usage:
        nextflow run main.nf [--cluster|--classify] --input_dir <path/to/dir> [--test_dir <path/to/dir> (required only if --classify used)]

Required arguments:
    Workflow (At least one of these)
        --cluster                           Data is run through clustering workflow, starting with an unknown dataset the pipeline produces \
                                            a distance matrix and pseudophylogenetic tree showing relatedness between samples.                 
        --classify                          Data is run through classification workflow, starting with separate training and test datasets, \
                                            a model is created using the training dataset and known class metadata> This model is used to \
                                            predict the classes of the unknown test dataset. Requires --test_dir argument.  
        --update                            New data is run added to the existing database, creating updated .KPopTwister and KPopTwisted files.

    Input (At least one of these)                  
        --input_dir                         Path to directory containing fasta/fastq files. Paired-end fastqs require "R1" and "R2" in filenames. \
                                            Gzipped files are allowed. \
                                            If --classify used this directory is the training dataset. \
                                            If --update used this directory is the new dataset used to update.
        --accession_list                    Supply a list of SRA IDs to download as input samples in the form of a text file, with one SRA per line.

        --test_dir                          Directory containing unseen test dataset. Only required if --classify workflow invoked. 
        --test_accession_list               Supply a list of SRA IDs to download as test samples in the form of a text file, with one SRA per line. Only required if --classify workflow invoked. 
        --twisted_file                      Full path to .KPopTwisted file. Only required for --update workflow.
        --twister_file                      Full path to .KPopTwister file. Only required for --update workflow.

Optional arguments:
    General arguments
        --kmer_len                          Length of k-mer to use when generating spectra [12] 
        --output_dir                        Path to output directory. If directory doesn't exist then a new directory will be created. [projectDir/results]
        --output_prefix                     Prefix for output files [output]
        --cpu_num                           Number of CPUs used per process [8] 
        --meta_data                         Tsv file with two required columns with defined headers; "fileName" and "class". \
                                            "fileName" is file name if a fasta or fasta.gz file, or file prefix if paired-end fastqs. E.g. sample1.fasta.gz if fasta file or \
                                            sample1 if sample1_R1.fastq.gz and sample1_R2.fastq.gz. Additional columns allowed
        --match_reference                   Full path to reference fasta file, used to select contigs that only match the supplied reference
        --min_contig_match_len              Minimum number of query contig base pairs that match the reference. Only used with --match_reference option [250]
        --min_contig_match_proportion       Minimum fraction of query contig base pairs that match reference. Only used with --match_reference option [0.6]
        --pred_class_num                    Specify the top n number of best predictions to be included in .KPopSummary file. E.g. 2 would choose the top two closest classes [all]
        --help                              Print help instructions

    Flash arguments (https://pubmed.ncbi.nlm.nih.gov/21903629/)
        --flash_minOverlap                  The minimum required overlap length between two reads to provide a confident overlap. Only used on fastq inputs. [20]
        --flash_maxOverlap                  Maximum overlap length expected in approximately 90% of read pairs. Only used on fastq inputs. [1000]
        --extra_flash                       Any additional arguments for flash. E.g. --extra_flash '-O -x 0.35'

    Megahit arguments (https://pubmed.ncbi.nlm.nih.gov/25609793/)
        --extra_megahit                     Any additional arguments for Megahit. E.g. --extra_megahit '--k-min 25'

    KPop arguments (https://www.biorxiv.org/content/10.1101/2022.06.22.497172v2)
        --extra_kpopCount                   Any additional arguments for KPopCount. E.g. --extra_kpopCount
        --extra_kpopCountDB                 Any additional arguments for KPopCountDB. E.g. --extra_kpopCountDB
        --extra_kpopTwist                   Any additional arguments for KPopTwist. E.g. --extra_kpopTwist
        --extra_kpopTwistDB                 Any additional arguments for KPopTwistDB. E.g. --extra_kpopTwistDB
        --extra_kpopPhylo                   Any additional arguments for KPopPhylo. E.g. --extra_kpopPhylo
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
         cpus used    : ${params.cpu_num}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         """
         .stripIndent()

// Print help message if --help invoked
if (params.help){
    helpMessage()
    exit 0
}

// Define errors
if (!params.cluster && !params.classify && !params.update){
    log.error"""--cluster, --classify or --update required for workflow. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.input_dir == "" && params.accession_list == ""){
    log.error"""--input_dir or --accession_list required for workflow. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.classify && params.test_dir == "" && params.test_accession_list == ""){
    log.error"""--test_dir or --test_accession_list required for --classify workflow. Use --help for more information.
    """.stripIndent()
    exit 0
}

if (params.input_dir != "") {
    // Channels created for all fastas, paired fastas and all paired fastqs in input_dir, error if no fasta or fastq files ound in input_dir
    FASTAS = "${params.input_dir}/*.{fasta,fa,fasta.gz,fa.gz}"
    FASTQS = "${params.input_dir}/*.{fastq,fq,fastq.gz,fq.gz}"
    PAIRED_FASTQS = "${params.input_dir}/*_{R1,R2}.{fastq,fq,fastq.gz,fq.gz}"

    Channel
        .fromPath( FASTAS )
        .map { it -> [[fileName: it.toString().split("/")[-1]], [file(it)]]}
        .set {fasta_files}

    Channel
        .fromFilePairs( PAIRED_FASTQS )
        .map { it -> [[fileName: it[0].toString()], it[1]]}
        .set {paired_fastq_files}

    Channel
        .fromPath( FASTQS )
        .map { it -> [[fileName: it.toString().split("/")[-1]], [file(it)]]}
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
        .filter{ if (it =~/[\[\]\?\*\ ]/) {throw new IllegalArgumentException("Value '$it' contains illegal character")} }

} else {
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
        .ifEmpty {exit 1, log.info "Cannot find path file ${meta_data}"}
        .splitCsv(header:true, sep: "\t")
        .map { row -> meta = [[fileName: row.fileName.toString().split("/")[-1]], [meta_class: row.class]] }
        .set {meta_file}
    fasta_files
        .join(meta_file)
        .map {it -> [it[2], it[1][0]]}
        .groupTuple(by: [0])
        .set {meta_fasta_files}
}

// Channels created for all fastas, paired fastas and all paired fastqs in test_dir, error if no files found in test_dir
if (params.test_dir != "") {
    TEST_FASTAS = "${params.test_dir}/*.{fasta,fa,fasta.gz,fa.gz}"
    TEST_FASTQS = "${params.test_dir}/*.{fastq,fq,fastq.gz,fq.gz}"
    TEST_PAIRED_FASTQS = "${params.test_dir}/*_{R1,R2}.{fastq,fq,fastq.gz,fq.gz}"
    
    if (!file( TEST_FASTAS ) && !file( TEST_FASTQS ) && !file( TEST_PAIRED_FASTQS && params.test_accession_list == "")){
        log.error"""--test_dir needs to contain at least one fasta or fastq file. Please check specified directory again. Use --help for more information.
        """.stripIndent()
        exit 0
    }

    Channel
    .fromPath( TEST_FASTAS )
    .map { it -> [[fileName: it.toString().split("/")[-1]], [file(it)]]}
    .set {test_fasta_files}

    Channel
    .fromFilePairs( TEST_PAIRED_FASTQS )
    .map { it -> [[fileName: it[0].toString()], it[1]]}
    .set {test_paired_fastq_files}

    Channel
    .fromPath( TEST_FASTQS )
    .map { it -> [[fileName: it.toString().split("/")[-1]], [file(it)]]}
    .filter { !(it.fileName =~ /(R1|R2)/) }
    .concat(test_paired_fastq_files)
    .set {test_fastq_files}

    test_fasta_files
    .concat(test_fastq_files)
    .map {it -> it[0].fileName}
    .filter{ if (it =~/[\[\]\?\*\ ]/) {throw new IllegalArgumentException("Value '$it' contains illegal character")} }

} else {
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
    // Download training set
    if (params.accession_list != "") {
        DOWNLOAD_SRAS1(params.accession_list)
            .flatten()
            .set {SRA_FILES}
        FASTERQ_DUMP1(SRA_FILES)
            .map {it -> [[fileName: it[0].toString()], it[1]]}
            .concat(fastq_files)
            .set {fastq_files}
    }

    // Download test set
    if (params.test_accession_list != "") {
        DOWNLOAD_SRAS2(params.test_accession_list)
            .flatten()
            .set {TEST_SRA_FILES}
        FASTERQ_DUMP2(TEST_SRA_FILES)
            .map {it -> [[fileName: it[0].toString()], it[1]]}
            .concat(test_fastq_files)
            .set {test_fastq_files}
    }

    // Clustering workflow
    if (params.cluster) {
        ASSEMBLE_FASTQS1(fastq_files)
            .map(it -> [it[0], [it[1]]])
            .set {assembled_fastas}
        ASSEMBLY_STATS1(assembled_fastas)
        assembled_fastas.concat(fasta_files)
            .set {concat_fasta_files}
        if (params.match_reference != "") {
            concat_fasta_files
                .map(it -> [it[1][0], it[0].fileName])
                .set {adjusted_concat_fasta_files}
            MATCH_REFERENCE_CONTIGS1(adjusted_concat_fasta_files)
                .map(it -> [[fileName: it[1]], [it[0]]])
                .set {concat_fasta_files}
        }
        COMBINE_FILES(concat_fasta_files)
            .collectFile(name: "${params.output_dir}/modified_fasta_files/${params.output_prefix}_combined.fasta.gz", newLine: false)
            .map(it -> [it, params.output_prefix])
            .set {combined_fasta}
        KPOPCOUNT(combined_fasta)
            .map(it -> [it, params.output_prefix])
            .set {kpopcount_file}
        KPOPTWIST(kpopcount_file)
            .map(it -> [it[0], it[1], params.output_prefix])
            .set {kpoptwist_files}
        KPOPPHYLO(kpoptwist_files).nwk_file
            .set {phylo_nwk_file}
        if (params.meta_data != "") {
            phylo_nwk_file
                .map(it -> [file(params.meta_data), it, params.output_prefix])
                .set {meta_ch}
            META_COLOURED_TREE(meta_ch)
        }
    }

    // Classification workflow
    if (params.classify) {
        ASSEMBLE_FASTQS1(fastq_files)
            .map(it -> [it[0], [it[1]]])
            .set {train_assembled_fastas}
        ASSEMBLY_STATS1(train_assembled_fastas)
        train_assembled_fastas
            .join(meta_file)
            .map {it -> [it[2], it[1][0]]}
            .concat(meta_fasta_files.transpose())
            .set {concat_meta_fasta_files}
        if (params.match_reference != "") {
            concat_meta_fasta_files
                .map(it -> [it[1], "${it[0].meta_class}_${it[1].getBaseName(file(it[1]).name.endsWith('.gz')? 2: 1)}"])
                .set {adjusted_concat_meta_fasta_files}
            MATCH_REFERENCE_CONTIGS1(adjusted_concat_meta_fasta_files)
                .map(it -> [[meta_class: it[1].toString().split("_")[0]], it[0]])
                .set {concat_meta_fasta_files}
        }
        concat_meta_fasta_files
            .groupTuple(by: [0])
            .set {assembled_meta_fasta_files}
        train_fasta_list = COMBINE_FILES_BY_CLASS(assembled_meta_fasta_files)
            .toSortedList()
        KPOPCOUNT_BY_CLASS(train_fasta_list)
            .map(it -> [it, "train"])
            .set {train_kpopcount_file}
        KPOPTWIST(train_kpopcount_file)
            .map(it -> [it, "train"])
            .set {train_kpoptwist_files}
        ASSEMBLE_FASTQS2(test_fastq_files)
            .set {test_assembled_fastas}
        ASSEMBLY_STATS2(test_assembled_fastas)
        test_assembled_fastas.concat(test_fasta_files)
            .set {test_concat_fasta_files}      
        COMBINE_FILES(test_concat_fasta_files)
            .collectFile(name: "${params.output_dir}/modified_fasta_files/test_combined.fasta.gz", newLine: false)
            .map(it -> [it, "test"])
            .set {test_combined_fasta}
        if (params.match_reference != "") {
            MATCH_REFERENCE_CONTIGS2(test_combined_fasta)
                .set {test_combined_fasta}
        }
        train_test_files = train_kpoptwist_files
            .combine(test_combined_fasta) // ####add meta file to ensure this bit doesn't go wrong####
            .map(it -> [it[0][0], it[0][1], it[2]])
        train_test_kpop_files = GENERATE_TEST_TWISTED(train_test_files)
        PREDICT_TEST_SET(train_test_kpop_files)
    }

        // Update workflow
    if (params.update) {
        ASSEMBLE_FASTQS1(fastq_files)
            .map(it -> [it[0], [it[1]]])
            .set {assembled_fastas}
        ASSEMBLY_STATS1(assembled_fastas)
        assembled_fastas.concat(fasta_files)
            .set {concat_fasta_files}
        if (params.match_reference != "") {
            concat_fasta_files
                .map(it -> [it[1][0], it[0].fileName])
                .set {adjusted_concat_fasta_files}
            MATCH_REFERENCE_CONTIGS1(adjusted_concat_fasta_files)
                .map(it -> [[fileName: it[1]], [it[0]]])
                .set {concat_fasta_files}
        }
        COMBINE_FILES(concat_fasta_files)
            .collectFile(name: "${params.output_dir}/modified_fasta_files/${params.output_prefix}_combined.fasta.gz", newLine: false)
            .map(it -> [it, params.output_prefix])
            .set {combined_fasta}
        KPOPCOUNT(combined_fasta)
            .map(it -> [it, params.output_prefix])
            .set {kpopcount_file}
        GENERATE_KPOPTWISTED(kpopcount_file)
            .set {updating_file}
        KPOPTWIST_UPDATE(updating_file)
    }
}
