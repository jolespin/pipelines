```
 _______ _______ _______  ______       _____  _____  _____  _______        _____ __   _ _______
 |______    |    |_____| |_____/      |_____]   |   |_____] |______ |        |   | \  | |______
 ______|    |    |     | |    \_      |       __|__ |       |______ |_____ __|__ |  \_| |______
 
```
### Pipeline:
1. Run `kneaddata` to quality trim via `Trimommatic` and remove contamination by mapping to `kneaddata_contamination_db` with `Bowtie2`
2. Clean up any mispaired reads with `repair.sh`
3. Map to exons using `STAR` and store unmapped reads
4. Convert `.sam` to `.bam`
5. Count reads per gene or transcript using `featureCounts`


### Requirements:
* [STAR](https://anaconda.org/bioconda/star) | v2.7.3a
* [BBMap](https://anaconda.org/bioconda/bbmap) (repair.sh) | v38.75
* [Kneaddata](https://anaconda.org/bioconda/kneaddata) | v0.6.1
* [Bowtie2](https://anaconda.org/bioconda/bowtie2) | v2.3.5
* [samtools ≥ v1.9](https://anaconda.org/bioconda/samtools) | v1.10
* [Subread](https://anaconda.org/bioconda/subread) (featureCounts) | v1.6.4
* [soothsayer_utils](https://github.com/jolespin/soothsayer_utils) | v2020.03.29
* [genopype](https://github.com/jolespin/genopype) | v2020.03.27

### Inputs:
* R1/R2 (fastq[.gz])
* Reference assembly (.fasta)
* Reference annotation (.gtf)
* STAR Index (Directory)
* Config [Default: Same location as script. Just edit this]

### Config: 
Tab-separated file that contains name and executable. 
Must contain at least 2 columns: `["name", "executable"]` but can contain more.  

For example, this is verbose but is version agnostic (if one program needs Python 2 while another needs Python 3, etc):

```
step	category	name	executable
1	preprocessing	repair.sh	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  &&  repair.sh
1	preprocessing	clumpify.sh	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  &&  clumpify.sh
1	preprocessing	kneaddata	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env && kneaddata
1	preprocessing	kneaddata_contamination_db	/usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/
2	mapping	STAR	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && STAR
2	mapping	samtools	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && samtools
3	quantify	featureCounts	source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && featureCounts
```

These are all in the same `conda` environment so a more concise version would be the following: 

```
name	executable
repair.sh	repair.sh
kneaddata	kneaddata
kneaddata_contamination_db	/usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/
STAR	STAR
samtools	samtools
featureCounts	featureCounts
```
Though, if the above config formatting is used then the script must be run from an environment that contains all of the required executables. (e.g. `conda activate star_env` or `source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env` to be safe)

### Usage: 
```
python /Users/jespinoz/Google/Algorithms/Pipelines/star_pipeline/star_pipeline.py -h
usage: star_pipeline.py -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --star_index <star_index/>

    Running: star_pipeline.py v2020.01.15 via Python v3.6.9 | /Users/jespinoz/anaconda3/envs/µ_env/bin/python

optional arguments:
  -h, --help            show this help message and exit

Required arguments:
  -1 R1, --r1 R1        path/to/r1.fq
  -2 R2, --r2 R2        path/to/r2.fq
  -n NAME, --name NAME  Name of sample
  -o PROJECT_DIRECTORY, --project_directory PROJECT_DIRECTORY
                        path/to/project_directory [Default: ./mapping_output]
  --ref_assembly REF_ASSEMBLY
                        path/to/reference.fasta
  --ref_annotation REF_ANNOTATION
                        path/to/reference.gtf
  --star_index STAR_INDEX
                        path/to/STAR_index

Utility arguments:
  --path_config PATH_CONFIG
                        path/to/star_pipeline_config.tsv
                        [Default: /Users/jespinoz/Google/Algorithms/Pipelines/star_pipeline/star_pipeline_config.tsv]
  --n_jobs N_JOBS       Number of threads [Default: 4]
  --random_state RANDOM_STATE
                        Random state [Default: 0]
  --preprocess_only     Only run preprocess steps (kneaddata)
  --restart_from_checkpoint RESTART_FROM_CHECKPOINT
                        Restart from a particular checkpoint [Default: None]
  -v, --version         show program's version number and exit

Kneaddata arguments:
  --kneaddata_options KNEADDATA_OPTIONS
                        Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home
  --bowtie2_options BOWTIE2_OPTIONS
                        Bowtie2 | More options (e.g. --arg 1 ) [Default: '']
                        http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

STAR arguments:
  --star_options STAR_OPTIONS
                        STAR | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/alexdobin/STAR

featureCounts arguments:
  --featurecounts_options FEATURECOUNTS_OPTIONS
                        featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/

Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org)
```

### Directory Structure:
`--project_directory Testing`
`--name Exp264_B9_FDApos_Infneg24_NTCERCCs_S15`

```
# View directory structure
tree Testing/
Testing/
└── Exp264_B9_FDApos_Infneg24_NTCERCCs_S15
    ├── checkpoints
    │   ├── 1__kneaddata
    │   ├── 2__star
    │   └── 3__featurecounts
    ├── commands.sh
    ├── intermediate
    │   ├── featurecounts_output
    │   │   ├── featurecounts.tsv.gz
    │   │   └── featurecounts.tsv.summary
    │   └── star_output
    │       ├── star_Aligned.out.bam
    │       ├── star_Log.final.out
    │       ├── star_Log.out
    │       ├── star_Log.progress.out
    │       ├── star_SJ.out.tab
    │       ├── star__STARtmp
    │       ├── star_Unmapped.out.mate1.gz
    │       └── star_Unmapped.out.mate2.gz
    ├── log
    │   ├── 1__kneaddata.e
    │   ├── 1__kneaddata.o
    │   ├── 1__kneaddata.returncode
    │   ├── 2__star.e
    │   ├── 2__star.o
    │   ├── 2__star.returncode
    │   ├── 3__featurecounts.e
    │   ├── 3__featurecounts.o
    │   └── 3__featurecounts.returncode
    ├── output
    ├── preprocessing
    │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_paired_contam_1.fastq.gz
    │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_paired_contam_2.fastq.gz
    │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_unmatched_1_contam.fastq.gz
    │   ├── kneaddata_GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index_bowtie2_unmatched_2_contam.fastq.gz
    │   ├── kneaddata.log
    │   ├── kneaddata_repaired_1.fastq.gz
    │   ├── kneaddata_repaired_2.fastq.gz
    │   ├── kneaddata_repaired_singletons.fastq.gz
    │   ├── kneaddata_unmatched_1.fastq.gz
    │   └── kneaddata_unmatched_2.fastq.gz
    └── tmp

10 directories, 32 files
```

### Output log file:
```
=============
STAR Pipeline
=============
--------------
Configuration:
--------------
............................................
Name: Exp264_B9_FDApos_Infneg24_NTCERCCs_S15
............................................
Python version: 3.7.6 | packaged by conda-forge | (default, Jan  7 2020, 22:33:48)  [GCC 7.3.0]
Python path: /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/star_env/bin/python
Script version: 2020.01.15
Moment: 2020-01-27 20:39:21
Directory: /local/ifs2_projdata/0568/projects/PLANKTON/illumina_aallen/jespinoz/Transcriptomes
Commands:
['/home/jespinoz/Algorithms/Pipelines/star_pipeline/star_pipeline.py', '--name', 'Exp264_B9_FDApos_Infneg24_NTCERCCs_S15', '-1', '../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R1_001.fastq.gz', '-2', '../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R2_001.fastq.gz', '-o', './Testing', '--ref_assembly', '../Modeling/funannotate_output/annotate_results/Cylindrotheca_closterium.scaffolds.ERCC92.virus.fa', '--ref_annotation', '../Modeling/funannotate_output/annotate_results/Cylindrotheca_closterium.ERCC92.virus.gtf', '--star_index', '../Modeling/funannotate_output/annotate_results/star_index_genome_ercc92_virus', '--n_jobs', '32']
-----------------------------------------------------------------------------------------------------------------------------------------------------
Adding executables to path from config (Overriding the conda executables): /home/jespinoz/Algorithms/Pipelines/star_pipeline/star_pipeline_config.tsv
-----------------------------------------------------------------------------------------------------------------------------------------------------
repair.sh --> source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  &&  repair.sh
kneaddata --> source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env && kneaddata
kneaddata_contamination_db --> /usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/
STAR --> source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && STAR
samtools --> source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && samtools
featureCounts --> source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && featureCounts

Executing pipeline:   0%|          | 0/3 [00:00<?, ?it/s]===========================
. .. ... Compiling ... .. .
===========================
Step: 1, kneaddata | log_prefix = 1__kneaddata | Quality trimming, removing contaminated reads, and optimizing file compression
Step: 2, star | log_prefix = 2__star | Aligning reads to reference
Step: 3, featurecounts | log_prefix = 3__featurecounts | Counting mapped reads
_________________________________________________________________________________
. .. ... STAR Mapping Pipeline || Exp264_B9_FDApos_Infneg24_NTCERCCs_S15 ... .. .
_________________________________________________________________________________

=============
. kneaddata .
=============
Input: ['../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R1_001.fastq.gz', '../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R2_001.fastq.gz']
Output: ['./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz', './Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz']

Command:
( source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env && kneaddata --input ../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R1_001.fastq.gz --input ../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R2_001.fastq.gz --reference-db /usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/ --output ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing --log ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata.log --threads 32 --output-prefix kneaddata --bowtie2-options="--seed 0"  && source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  &&  repair.sh in1=./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_paired_1.fastq in2=./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_paired_2.fastq out1=./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz out2=./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz outs=./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_singletons.fastq.gz overwrite=t && rm ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_paired* ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/*trimmed* && pigz -f -p 32 ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/*.fastq )

Validating the following input files:
[=] File exists (111615 bytes): ../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R1_001.fastq.gz
[=] File exists (116779 bytes): ../Fastq/SOLEXASEQ-1411/seq_raw/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15_R2_001.fastq.gz

Running. .. ... .....
Executing pipeline:  33%|███▎      | 1/3 [00:10<00:21, 10.71s/it]Log files:
./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/log/1__kneaddata.*

Validating the following output files:
[=] File exists (34228 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (35793 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz

Duration: 00:00:10

========
. star .
========
Input: ['./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz', './Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz']
Output: ['./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam']

Command:
source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && STAR --genomeDir ../Modeling/funannotate_output/annotate_results/star_index_genome_ercc92_virus --readFilesCommand zcat --readFilesIn ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz --outFileNamePrefix ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_ --runThreadN 32 --outReadsUnmapped Fastx  && source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && samtools view -h -b ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.sam > ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam && rm ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.sam && pigz -f -p 32 ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Unmapped.out.mate1 && pigz -f -p 32 ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Unmapped.out.mate2

Validating the following input files:
[=] File exists (34228 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_1.fastq.gz
[=] File exists (35793 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/preprocessing/kneaddata_repaired_2.fastq.gz

Running. .. ... .....
Executing pipeline:  67%|██████▋   | 2/3 [00:23<00:11, 11.32s/it]Log files:
./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/log/2__star.*

Validating the following output files:
[=] File exists (4869 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam

Duration: 00:00:23

=================
. featurecounts .
=================
Input: ['./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam']
Output: ['./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/featurecounts_output/featurecounts.tsv.gz']

Command:
source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env  && featureCounts -G ../Modeling/funannotate_output/annotate_results/Cylindrotheca_closterium.scaffolds.ERCC92.virus.fa -a ../Modeling/funannotate_output/annotate_results/Cylindrotheca_closterium.ERCC92.virus.gtf -o ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/featurecounts_output/featurecounts.tsv -F GTF --tmpDir ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/tmp -T 32  ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam && gzip ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/featurecounts_output/featurecounts.tsv

Validating the following input files:
[=] File exists (4869 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/star_output/star_Aligned.out.bam

Running. .. ... .....
Executing pipeline: 100%|██████████| 3/3 [00:25<00:00,  8.61s/it]
Log files:
./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/log/3__featurecounts.*

Validating the following output files:
[=] File exists (478810 bytes): ./Testing/Exp264_B9_FDApos_Infneg24_NTCERCCs_S15/intermediate/featurecounts_output/featurecounts.tsv.gz

Duration: 00:00:26


........................
Total duration: 00:00:26
........................
```

### JCVI Internal Users:
```
# Activate conda environment
source /usr/local/devel/ANNOTATION/jespinoz/anaconda3/bin/activate star_env

# Run executable
star_pipeline.py -h
```
