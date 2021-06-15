```
 _     _ _____ _______ _______ _______       _____  _____  _____  _______        _____ __   _ _______
 |_____|   |   |______ |_____|    |         |_____]   |   |_____] |______ |        |   | \  | |______
 |     | __|__ ______| |     |    |         |       __|__ |       |______ |_____ __|__ |  \_| |______
                                                                                                                                                                                                                  
```
### Pipeline:
1. Run `kneaddata` to quality trim via `Trimommatic` and remove contamination by mapping to `kneaddata_contamination_db` with `Bowtie2`
2. Clean up any mispaired reads with `repair.sh`
3. Map to reference using `HISAT2` and store unmapped reads
4. Convert `.sam` to sorted `.bam` via `samtools`
5. Count reads per gene or transcript using `featureCounts`


### Requirements:
* [BBMap](https://anaconda.org/bioconda/bbmap) (repair.sh) ≥ v38.75
* [Kneaddata](https://anaconda.org/bioconda/kneaddata) ≥ v0.6.1
* [HISAT2](https://anaconda.org/bioconda/hisat2) ≥ v2.2.1
* [samtools ≥ v1.9](https://anaconda.org/bioconda/samtools) ≥ v1.10
* [Subread](https://anaconda.org/bioconda/subread) (featureCounts) ≥ v1.6.4
* [soothsayer_utils](https://github.com/jolespin/soothsayer_utils) ≥ v2020.03.29
* [genopype](https://github.com/jolespin/genopype) ≥ v2020.03.27

### Set up environment: 

```
conda env create -n bowtie2_env -f environment.yml -y
```

### Inputs:
* R1/R2 (fastq[.gz])
* Reference assembly (.fasta)
* Reference annotation (.gtf)
* HISAT2 Index (Prefix)

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
(bowtie2_env) -bash-4.1$ bowtie2_pipeline.py -h
usage: bowtie2_pipeline.py -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --bowtie2_index <bowtie2_index/>

    Running: bowtie2_pipeline.py v2021.04.04 via Python v3.8.8 | /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/bowtie2_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -1 R1, --r1 R1        path/to/r1.fq
  -2 R2, --r2 R2        path/to/r2.fq
  -12 INTERLEAVED_READS, --interleaved_reads INTERLEAVED_READS
                        path/to/interleaved.fq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: ./mapping_output]
  --ref_assembly REF_ASSEMBLY
                        path/to/reference.fasta
  --ref_annotation REF_ANNOTATION
                        path/to/reference.gtf
  --bowtie2_index BOWTIE2_INDEX
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
  --attribute_type ATTRIBUTE_TYPE
                        Attribute type in GTF/GFF file. Use 'ID' for prodigal. [Default: gene_id]
  --feature_type FEATURE_TYPE
                        Feature type in GTF/GFF file. Use 'CDS' for prodigal. Use 'gene' for prokaryotic genomes from NCBI. [Default: exon]
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)
```

### Directory Structure:
`--project_directory hisat2_output`

`--name test`

```
# View directory structure (Note, this run had no annotation file hence the genome.saf file)
(hisat2_env) -bash-4.2$ tree hisat2_output/test/
hisat2_output/test/
├── checkpoints
│   ├── 0__kneaddata
│   ├── 1__hisat2
│   ├── 2__featurecounts
│   └── 3__symlink
├── commands.sh
├── intermediate
│   ├── featurecounts_output
│   │   ├── featurecounts.tsv.gz
│   │   ├── featurecounts.tsv.summary
│   │   └── genome.saf
│   └── hisat2_output
│       ├── mapped.sorted.bam
│       ├── unmapped_1.fastq.gz
│       └── unmapped_2.fastq.gz
├── log
│   ├── 0__kneaddata.e
│   ├── 0__kneaddata.o
│   ├── 0__kneaddata.returncode
│   ├── 1__hisat2.e
│   ├── 1__hisat2.o
│   ├── 1__hisat2.returncode
│   ├── 2__featurecounts.e
│   ├── 2__featurecounts.o
│   ├── 2__featurecounts.returncode
│   ├── 3__symlink.e
│   ├── 3__symlink.o
│   └── 3__symlink.returncode
├── output
│   ├── featurecounts.tsv.gz -> /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz
│   ├── featurecounts.tsv.summary -> /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.summary
│   └── mapped.sorted.bam -> /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam
├── preprocessing
│   ├── kneaddata.log
│   ├── kneaddata_repaired_1.fastq.gz
│   ├── kneaddata_repaired_2.fastq.gz
│   └── kneaddata_repaired_singletons.fastq.gz
└── tmp

8 directories, 30 files
```

### Output log file:
```
(hisat2_env) -bash-4.2$ hisat2_pipeline.py -1 reads_1.fq.gz -2 reads_2.fq.gz -n test --ref_assembly test.fa --n_jobs 4
================
HISAT2 Pipeline
================
--------------
Configuration:
--------------
..........
Name: test
..........
Python version: 3.8.10 | packaged by conda-forge | (default, May 11 2021, 07:01:05)  [GCC 9.3.0]
Python path: /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/python
Script version: 2021.06.15
Moment: 2021-06-15 17:24:17
Directory: /local/ifs3_scratch/CORE/jespinoz/Oral/testing
Commands:
['/usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/hisat2_pipeline.py', '-1', 'reads_1.fq.gz', '-2', 'reads_2.fq.gz', '-n', 'test', '--ref_assembly', 'test.fa', '--n_jobs', '4']
No --ref_annotation has been provided.  Converting FASTA -> SAF and aggregating counts per sequence record
------------------------------------------------------------------
Adding executables to path from the following source: CONDA_PREFIX
------------------------------------------------------------------
hisat2 --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/hisat2
featureCounts --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/featureCounts
kneaddata --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/kneaddata
repair.sh --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/repair.sh
samtools --> /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/samtools

===========================
. .. ... Compiling ... .. .
===========================
Step: 0, kneaddata | log_prefix = 0__kneaddata | Quality trimming, removing contaminated reads, and optimizing file compression
Step: 1, hisat2 | log_prefix = 1__hisat2 | Aligning reads to reference
Step: 2, featurecounts | log_prefix = 2__featurecounts | Counting reads
Step: 3, symlink | log_prefix = 3__symlink | Symlinking relevant output files
_________________________________________________
. .. ... HISAT2 Mapping Pipeline || test ... .. .
_________________________________________________

Executing pipeline:   0%|                                                                                                                                                | 0/4 [00:00<?, ?it/s]=============
. kneaddata .
=============
Input: ['reads_1.fq.gz', 'reads_2.fq.gz']
Output: ['./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz', './hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz']

Command:
( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/kneaddata --input reads_1.fq.gz --input reads_2.fq.gz  --output ./hisat2_output/test/preprocessing --log ./hisat2_output/test/preprocessing/kneaddata.log --threads 4 --output-prefix kneaddata --bowtie2-options="--seed 0"  ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/repair.sh in1=./hisat2_output/test/preprocessing/kneaddata.trimmed.1.fastq in2=./hisat2_output/test/preprocessing/kneaddata.trimmed.2.fastq out1=./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz out2=./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz outs=./hisat2_output/test/preprocessing/kneaddata_repaired_singletons.fastq.gz overwrite=t ) && ( rm -f ./hisat2_output/test/preprocessing/kneaddata_paired* ./hisat2_output/test/preprocessing/*trimmed* )

Validating the following input files:
[=] File exists (8 MB): reads_1.fq.gz
[=] File exists (8 MB): reads_2.fq.gz

Running. .. ... .....

Log files:
./hisat2_output/test/log/0__kneaddata.*

Validating the following output files:
[=] File exists (8 MB): ./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (8 MB): ./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz

Duration: 00:00:05

Executing pipeline:  25%|██████████████████████████████████                                                                                                      | 1/4 [00:05<00:16,  5.62s/it]==========
. hisat2 .
==========
Input: ['./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz', './hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz']
Output: ['./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam']

Command:
rm -rf ./hisat2_output/test/tmp/* && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/hisat2 -x test.fa -1 ./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz -2 ./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz --threads 4 --un-conc-gz ./hisat2_output/test/intermediate/hisat2_output/unmapped_%.fastq.gz --seed 0  ) | ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/samtools sort --threads 4 --reference test.fa -T ./hisat2_output/test/tmp/samtools_sort > ./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam )

Validating the following input files:
[=] File exists (8 MB): ./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (8 MB): ./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz

Running. .. ... .....

Log files:
./hisat2_output/test/log/1__hisat2.*

Validating the following output files:
[=] File exists (16 MB): ./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam

Duration: 00:00:12

Executing pipeline:  50%|████████████████████████████████████████████████████████████████████                                                                    | 2/4 [00:12<00:12,  6.48s/it]=================
. featurecounts .
=================
Input: ['./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam']
Output: ['./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz']

Command:
(python /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/fasta_to_saf.py -i test.fa > ./hisat2_output/test/intermediate/featurecounts_output/genome.saf ) && ( /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/bin/featureCounts -G test.fa -a ./hisat2_output/test/intermediate/featurecounts_output/genome.saf -o ./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv -F SAF --tmpDir ./hisat2_output/test/tmp -T 4  ./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam ) && gzip ./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv

Validating the following input files:
[=] File exists (16 MB): ./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam

Running. .. ... .....

Log files:
./hisat2_output/test/log/2__featurecounts.*

Validating the following output files:
[=] File exists (4494 bytes): ./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz

Duration: 00:00:14

Executing pipeline:  75%|██████████████████████████████████████████████████████████████████████████████████████████████████████                                  | 3/4 [00:14<00:04,  4.20s/it]===========
. symlink .
===========
Input: ['./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam', './hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz', './hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.summary']
Output: ['./hisat2_output/test/output/mapped.sorted.bam', './hisat2_output/test/output/featurecounts.tsv.gz', './hisat2_output/test/output/featurecounts.tsv.summary']

Command:
( ln -f -s /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/output && ln -f -s /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/output && ln -f -s /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.summary /local/ifs3_scratch/CORE/jespinoz/Oral/testing/hisat2_output/test/output )

Validating the following input files:
[=] File exists (16 MB): ./hisat2_output/test/intermediate/hisat2_output/mapped.sorted.bam
[=] File exists (4494 bytes): ./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.gz
[=] File exists (412 bytes): ./hisat2_output/test/intermediate/featurecounts_output/featurecounts.tsv.summary

Running. .. ... .....

Log files:
./hisat2_output/test/log/3__symlink.*

Duration: 00:00:14

Executing pipeline: 100%|████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████| 4/4 [00:14<00:00,  3.55s/it]

........................
Total duration: 00:00:14
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