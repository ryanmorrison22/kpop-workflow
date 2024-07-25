# kpop-workflow
A simplified workflow to run KPop (https://github.com/PaoloRibeca/KPop) for clustering and classification

To run the workflow you'll first need to create the kpop environment using the enviroment.yml file and activate the kpop environment 
```
conda env create -f environment.yml
conda activate kpop​
```

Currently, there are three main workflows:
1. Clustering workflow - starting with an unknown dataset, running KPop to generate transformational database files (.KPopTwister and .KPopTwisted) and a pseudo-phylogenetic tree showing relatedness:
```
nextflow run main.nf --cluster --input_dir /dir/containing/fastqAndFastas
```

2. Classifying workflow - starting with a known training dataset with meta data, such as classes/species, and predicting the class/species of an unknown test set:
```
nextflow run main.nf --classify --input_dir /dir/containing/training/fastqAndFastas --meta_data /path/to/meta/file --test_dir /dir/containing/test/fastqAndFastas
```

3. Update transformational database workflow - starting with an unknown dataset and previous transformational database files, the new data is added to the existing database, creating updated .KPopTwister and KPopTwisted files and a new pseudo-phylogenetic tree.
```
nextflow run main.nf --update --test_dir /dir/containing/new/fastqAndFastas --twister_file /path/to/KPopTwister/file --twisted_file /path/to/KPopTwisted/file
```

All workflows create lots of files, but the most important outputs are the database files (.KPopTwister and .KPopTwisted in KPopTwist_files) and KPop pseudo-phylogenetic tree found in results/trees_and_metrics/output_2.<nwk/pdf> (--cluster), the class/species predictions found in results/predictions/output.<predictions/KPopSummary>.txt (--classify) and the updated database files (<output_prefix>.KPopTwisted and <output_prefix>.KPopTwisted) and updated_comparison.pdf plot in updated_KPopTwist_files (--update).

**Workflows**
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--cluster` |  | Data is run through clustering workflow, starting with an unknown dataset the pipeline produces a distance matrix and pseudophylogenetic tree showing relatedness between samples. |  |         
| `--classify` |  | Data is run through classification workflow, starting with separate training and test datasets, a model is created using the training dataset and known class metadata> This model is used to predict the classes of the unknown test dataset. Requires --test_dir argument. |  |
| `--update` |  | New data is run added to the existing database, creating updated .KPopTwister and KPopTwisted files. |  |

**Input**  
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--input_dir` | _directory\_name_ | Path to directory containing fasta/fastq files. Paired-end fastqs require "R1" and "R2" in filenames. Gzipped files are allowed. If --classify used this directory is the training dataset. |  |
| `--test_dir` | _directory\_name_ | Directory containing unseen test dataset. Only required if --classify workflow invoked. |  |
| `--accession_list` | _txt\_filename_ | Supply a list of SRA IDs to download as input samples in the form of a text file, with one SRA per line. |  |
| `--test_accession_list` | _txt\_filename_ | Supply a list of SRA IDs to download as test samples in the form of a text file, with one SRA per line. Only required if --classify workflow invoked. |  |
| `--meta_data` | _TSV\_filename_ | Tsv file with two required columns with defined headers; "fileName" and "class". "fileName" is file name if a fasta or fasta.gz file, or file prefix if paired-end fastqs. E.g. sample1.fasta.gz if fasta file or sample1 if sample1_R1.fastq.gz and sample1_R2.fastq.gz. Only required if --classify workflow invoked. If used with --cluster workflow, it will generate a pseudo-phylogenetic tree coloured with metadata. Additional columns allowed. |  | 
| `--twisted_file` | _.KPopTwisted\_file_ | Full path to .KPopTwisted file. Only required for --update workflow. |  |
| `--twister_file` | _.KPopTwister\_file_ | Full path to .KPopTwister file. Only required for --update workflow. |  |

**Output**
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--output_dir` | _directory\_name_ | Path to output directory. If directory doesn't exist then a new directory will be created. | <ins>default=_projectDir\/results_</ins> |
| `--output_prefix` | _prefix\_name_ | Prefix for output files | <ins>default=_output_</ins> |
| `--pred_class_num` | _positive\_integer \| all_ | Specify the top n number of best predictions to be included in .KPopSummary file. E.g. 2 would choose the top two closest classes | <ins>default=_all_</ins> |

**Reference matching**
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--match_reference` | _fasta\_filename_ | Full path to reference fasta file. Used to select contigs that only match the supplied reference. |  |
| `--min_contig_match_len` | _positive\_integer_ | Minimum number of query contig base pairs that match the reference. Only used with --match_reference option | <ins>default=_250_</ins> |
| `--min_contig_match_proportion` | _fractional\_float_ | Minimum fraction of query contig base pairs that match reference. Only used with --match_reference option | <ins>default=_0.6_</ins> |
        
**General arguments**
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--kmer_len` | _positive\_integer_ | Length of k-mer to use when generating spectra | <ins>default=_12_</ins> |
| `--cpu_num` | _positive\_integer_ | Number of CPUs used per process | <ins>default=_8_</ins> |
| `--help` |  | Print help instructions |  |

**Additional arguments**
| Option | Argument(s) | Effect | Note(s) |
|-|-|-|-|
| `--tree_type` | _string_ | Specify the type of tree generated by ggtree - 'rectangular' or 'circular' | <ins>default=_rectangular_</ins> |
| `--tree_label_size` | _positive\_integer_ | Specify the size of the labels on the tree generated by ggtree, choose 0 to remove labels | <ins>default=_3_</ins> |
| `--extra_flash` | _string_ | Any additional arguments for flash (https://pubmed.ncbi.nlm.nih.gov/21903629/). E.g. --extra_flash '-O -x 0.35' |  |
| `--flash_minOverlap` | _positive\_integer_ | The minimum required overlap length between two reads to provide a confident overlap. Only used on fastq inputs. | <ins>default=_20_</ins> |
| `--flash_maxOverlap` | _positive\_integer_ | Maximum overlap length expected in approximately 90% of read pairs. Only used on fastq inputs. | <ins>default=_1000_</ins> |
| `--extra_megahit` | _string_ | Any additional arguments for Megahit (https://pubmed.ncbi.nlm.nih.gov/25609793/). E.g. --extra_megahit '--k-min 25' |  |
| `--extra_kpopCount` | _string_ | Any additional arguments for KPopCount ((https://www.biorxiv.org/content/10.1101/2022.06.22.497172v2). E.g. --extra_kpopCount |  |
| `--extra_kpopCountDB` | _string_ | Any additional arguments for KPopCountDB. E.g. --extra_kpopCountDB |  |
| `--extra_kpopTwist` | _string_ | Any additional arguments for KPopTwist. E.g. --extra_kpopTwist |  |
| `--extra_kpopTwistDB` | _string_ | Any additional arguments for KPopTwistDB. E.g. --extra_kpopTwistDB |  |
| `--kpopphylo_power` | _positive\_integer_ | Set the external power when computing distances. | <ins>default=_2_</ins> |
| `--kpopphylo_distance` | _string_ | Distance measure to be used. This must be one of 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary' or 'minkowski'. | <ins>default=_euclidean_</ins> |
| `--kpopphylo_magic` | _string_ | Cluster-related variable (Not currently implemented). | <ins>default=_1._</ins> |
| `--extra_kpopPhylo` | _string_ | Any additional arguments for KPopPhylo. E.g. --extra_kpopPhylo |  |
