```
 ______   _____  _  _  _ _______ _____ _______       _____  _____  _____  _______        _____ __   _ _______
 |_____] |     | |  |  |    |      |   |______      |_____]   |   |_____] |______ |        |   | \  | |______
 |_____] |_____| |__|__|    |    __|__ |______      |       __|__ |       |______ |_____ __|__ |  \_| |______
                                                                                                             
```
### Pipeline:
1. Run `kneaddata` to quality trim via `Trimommatic` and remove contamination by mapping to `kneaddata_contamination_db` with `Bowtie2`
2. Clean up any mispaired reads with `repair.sh`
3. Map to reference using `Bowtie2` and store unmapped reads
4. Convert `.sam` to sorted `.bam` via `samtools`
5. Count reads per gene or transcript using `featureCounts`


### Requirements:
* [BBMap](https://anaconda.org/bioconda/bbmap) (repair.sh) | v38.75
* [Kneaddata](https://anaconda.org/bioconda/kneaddata) | v0.6.1
* [Bowtie2](https://anaconda.org/bioconda/bowtie2) | v2.3.5
* [samtools ≥ v1.9](https://anaconda.org/bioconda/samtools) | v1.10
* [Subread](https://anaconda.org/bioconda/subread) (featureCounts) | v1.6.4
* [soothsayer_utils](https://github.com/jolespin/soothsayer_utils) | v2020.03.29
* [genopype](https://github.com/jolespin/genopype) | v2020.03.27

### Set up environment: 

```
conda env create -n bowtie2_env -f environment.yml -y
```

### Inputs:
* R1/R2 (fastq[.gz])
* Reference assembly (.fasta)
* Reference annotation (.gtf)
* Bowtie2 Index (Prefix)

### Accessory Scripts: 
* `fasta_to_saf.py` - Converts fasta to SAF format for mapping without an annotation (e.g. contigs)
* `merge_featurecounts.py` - Merges featureCounts counts or summary tables

### Config: 
Tab-separated file that contains name and executable. 
Must contain at least 2 columns: `["name", "executable"]` but can contain more.  

For example, this is verbose but is version agnostic (if one program needs Python 2 while another needs Python 3, etc):

```
step	category	name	executable
1	preprocessing	repair.sh	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env  &&  repair.sh
1	preprocessing	clumpify.sh	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env  &&  clumpify.sh
1	preprocessing	kneaddata	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env && kneaddata
2	mapping	bowtie2	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env  && STAR
2	mapping	samtools	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env  && samtools
3	quantify	featureCounts	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env  && featureCounts
```
Though, if the above config formatting is used then the script must be run from an environment that contains all of the required executables. (e.g. `conda activate bowtie2_env` or `source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env` to be safe)

These are all in the same `conda` environment so if this is the case, please just use `--path_config CONDA_PREFIX` to use the `bin` of the current `conda` environment. 



### Usage: 
```
(bowtie2_env) -bash-4.2$ bowtie2_pipeline.py -h
usage: bowtie2_pipeline.py -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --bowtie2_index <bowtie2_index/>

    Running: bowtie2_pipeline.py v2021.06.15 via Python v3.8.10 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -1 R1, --r1 R1        path/to/r1.fq
  -2 R2, --r2 R2        path/to/r2.fq
  -12 INTERLEAVED_READS, --interleaved_reads INTERLEAVED_READS
                        path/to/interleaved.fq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: ./bowtie2_output]
  -R REF_ASSEMBLY, --ref_assembly REF_ASSEMBLY
                        path/to/reference.fasta
  -A REF_ANNOTATION, --ref_annotation REF_ANNOTATION
                        path/to/reference.gtf
  -I BOWTIE2_INDEX, --bowtie2_index BOWTIE2_INDEX
                        path/to/bowtie2_index

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/config.tsv [Default: CONDA_PREFIX]
  --n_jobs N_JOBS       Number of threads [Default: 1]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --preprocess_only     Only run preprocess steps (kneaddata)
  --skip_preprocess     Don't run preprocess steps (kneaddata)
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Kneaddata arguments:
  --kneaddata_contamination_db KNEADDATA_CONTAMINATION_DB
                        Kneaddata | path/to/contamination_database [Default: '']
                        For human at JCVI, use the following: /usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/
  --kneaddata_options KNEADDATA_OPTIONS
                        Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home
  --kneaddatabowtie2_options KNEADDATABOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

Bowtie2 arguments:
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: ''] | http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

featureCounts arguments:
  -g ATTRIBUTE_TYPE, --attribute_type ATTRIBUTE_TYPE
                        Attribute type in GTF/GFF file. Use 'ID' for prodigal. [Default: gene_id]
  -t FEATURE_TYPE, --feature_type FEATURE_TYPE
                        Feature type in GTF/GFF file. Use 'CDS' for prodigal. Use 'gene' for prokaryotic genomes from NCBI. [Default: exon]
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)
```

### Directory Structure:
`--project_directory mapping_output`

`--name 5UM_FE_CITRATE_A-PE-P65-1-579_S1`

```
# View directory structure
tree mapping_output
mapping_output
├── 5UM_FE_CITRATE_A-PE-P65-1-579_S1
│   ├── checkpoints
│   │   ├── 1__kneaddata
│   │   ├── 2__bowtie2
│   │   ├── 3__featurecounts
│   │   └── 4__symlink
│   ├── commands.sh
│   ├── intermediate
│   │   ├── bowtie2_output
│   │   │   ├── mapped.sorted.bam
│   │   │   ├── unmapped_1.fastq.gz
│   │   │   └── unmapped_2.fastq.gz
│   │   └── featurecounts_output
│   │       ├── featurecounts.tsv.gz
│   │       └── featurecounts.tsv.summary
│   ├── log
│   │   ├── 1__kneaddata.e
│   │   ├── 1__kneaddata.o
│   │   ├── 1__kneaddata.returncode
│   │   ├── 2__bowtie2.e
│   │   ├── 2__bowtie2.o
│   │   ├── 2__bowtie2.returncode
│   │   ├── 3__featurecounts.e
│   │   ├── 3__featurecounts.o
│   │   ├── 3__featurecounts.returncode
│   │   ├── 4__symlink.e
│   │   ├── 4__symlink.o
│   │   └── 4__symlink.returncode
│   ├── output
│   │   └── featurecounts.tsv.gz -> /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz
│   ├── preprocessing
│   │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_paired_contam_1.fastq.gz
│   │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_paired_contam_2.fastq.gz
│   │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_unmatched_1_contam.fastq.gz
│   │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_unmatched_2_contam.fastq.gz
│   │   ├── kneaddata.log
│   │   ├── kneaddata_repaired_1.fastq.gz
│   │   ├── kneaddata_repaired_2.fastq.gz
│   │   ├── kneaddata_repaired_singletons.fastq.gz
│   │   ├── kneaddata_unmatched_1.fastq.gz
│   │   └── kneaddata_unmatched_2.fastq.gz
│   └── tmp
```

### Output log file:
```
================
Bowtie2 Pipeline
================
--------------
Configuration:
--------------
......................................
Name: 5UM_FE_CITRATE_A-PE-P65-1-579_S1
......................................
Python version: 3.8.6 | packaged by conda-forge | (default, Nov 27 2020, 19:31:52)  [GCC 9.3.0]
Python path: /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/python
Script version: 2021.04.04
Moment: 2021-04-05 17:52:29
Directory: /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas
Commands:
['/home/jespinoz/Algorithms/Pipelines/bowtie2_pipeline/bowtie2_pipeline.py', '-1', 'Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R1_001.fastq.gz', '-2', 'Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R2_001.fastq.gz', '-n', '5UM_FE_CITRATE_A-PE-P65-1-579_S1', '-o', 'mapping_output3', '--ref_assembly', 'Genome/GCF_000172635.2_ASM17263v2_genomic.fna', '--ref_annotation', 'Genome/GCF_000172635.2_ASM17263v2_genomic.gtf', '--n_jobs', '48', '--feature_type', 'gene']
------------------------------------------------------------------
Adding executables to path from the following source: CONDA_PREFIX
------------------------------------------------------------------
samtools --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/samtools
featureCounts --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/featureCounts
bowtie2 --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/bowtie2
repair.sh --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/repair.sh
kneaddata --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/kneaddata

Executing pipeline:   0%|          | 0/4 [00:00<?, ?it/s]===========================
. .. ... Compiling ... .. .
===========================
Step: 1, kneaddata | log_prefix = 1__kneaddata | Quality trimming, removing contaminated reads, and optimizing file compression
Step: 2, bowtie2 | log_prefix = 2__bowtie2 | Aligning reads to reference
Step: 3, featurecounts | log_prefix = 3__featurecounts | Counting reads
Step: 4, symlink | log_prefix = 4__symlink | Symlinking relevant output files
______________________________________________________________________________
. .. ... Bowtie2 Mapping Pipeline || 5UM_FE_CITRATE_A-PE-P65-1-579_S1 ... .. .
______________________________________________________________________________

=============
. kneaddata .
=============
Input: ['Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R1_001.fastq.gz', 'Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R2_001.fastq.gz']
Output: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/kneaddata --input Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R1_001.fastq.gz --input Fastq/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R2_001.fastq.gz  --output mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing --log mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata.log --threads 48 --output-prefix kneaddata --bowtie2-options="--seed 0"  ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/repair.sh in1=mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata.trimmed.1.fastq in2=mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata.trimmed.2.fastq out1=mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz out2=mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz outs=mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_singletons.fastq.gz overwrite=t ) && ( rm -f mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_paired* mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/*trimmed* && pigz -f -p 48 mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/*.fastq )

Validating the following input files:
[=] File exists (1328 MB): /usr/local/JTC/data/prod/solexa_depot/210324_A00588_0064_BHWF22DRXX/Unaligned_Lane1_25-03-2021_2/SOLEXASEQ-1563/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R1_001.fastq.gz
[=] File exists (1442 MB): /usr/local/JTC/data/prod/solexa_depot/210324_A00588_0064_BHWF22DRXX/Unaligned_Lane1_25-03-2021_2/SOLEXASEQ-1563/5UM_FE_CITRATE_A-PE-P65-1-579_S1_R2_001.fastq.gz

Running. .. ... .....
Executing pipeline:  25%|██▌       | 1/4 [06:04<18:14, 364.72s/it]
Log files:
mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/log/1__kneaddata.*

Validating the following output files:
[=] File exists (882 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (957 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz

Duration: 00:06:04

===========
. bowtie2 .
===========
Input: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz']
Output: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam']

Command:
rm -rf mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/tmp/* && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/bowtie2 -x Genome/GCF_000172635.2_ASM17263v2_genomic.fna -1 mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz -2 mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz --threads 48 --un-conc-gz mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/unmapped_%.fastq.gz --seed 0  ) | ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/samtools sort --threads 48 --reference Genome/GCF_000172635.2_ASM17263v2_genomic.fna -T mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/tmp/samtools_sort > mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam )

Validating the following input files:
[=] File exists (882 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (957 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/preprocessing/kneaddata_repaired_2.fastq.gz

Running. .. ... .....
Executing pipeline:  50%|█████     | 2/4 [26:18<20:38, 619.45s/it]
Log files:
mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/log/2__bowtie2.*

Validating the following output files:
[=] File exists (1184 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam

Duration: 00:26:18

=================
. featurecounts .
=================
Input: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam']
Output: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/featureCounts -G Genome/GCF_000172635.2_ASM17263v2_genomic.fna -a Genome/GCF_000172635.2_ASM17263v2_genomic.gtf -o mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv -F GTF --tmpDir mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/tmp -T 48 -g gene_id -t gene  mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam ) && gzip mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv

Validating the following input files:
[=] File exists (1184 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam

Running. .. ... .....
Executing pipeline:  75%|███████▌  | 3/4 [26:22<07:14, 434.73s/it]
Log files:
mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/log/3__featurecounts.*

Validating the following output files:
[=] File exists (63508 bytes): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz

Duration: 00:26:22

===========
. symlink .
===========
Input: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.summary']
Output: ['mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output/mapped.sorted.bam', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output/featurecounts.tsv.gz', 'mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output/featurecounts.tsv.summary']

Command:
( ln -f -s /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output && ln -f -s /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output && ln -f -s /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.summary /local/ifs3_scratch/METAGENOMICS/jespinoz/Alteromonas/mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/output )

Validating the following input files:
[=] File exists (1184 MB): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/bowtie2_output/mapped.sorted.bam
[=] File exists (63508 bytes): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.gz
[=] File exists (455 bytes): mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/intermediate/featurecounts_output/featurecounts.tsv.summary

Running. .. ... .....
Executing pipeline: 100%|██████████| 4/4 [26:22<00:00, 395.57s/it]

Log files:
mapping_output3/5UM_FE_CITRATE_A-PE-P65-1-579_S1/log/4__symlink.*

Duration: 00:26:22


........................
Total duration: 00:26:22
........................
```

### JCVI Internal Users:
```
# Activate conda environment
source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate bowtie2_env

# Run executable
bowtie2_pipeline.py -h
```

### Pending: 
* Take in manifest file for batch jobs
* Cannot handle interleaved reads
* No option to use unpaired reads