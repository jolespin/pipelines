java -ea -Xmx1584m -cp /usr/local/devel/ANNOTATION/jespinoz/anaconda3/envs/hisat2_env/opt/bbmap-38.90-1/current/ jgi.SplitPairsAndSingles rp in1=./hisat2_output/test/preprocessing/kneaddata.trimmed.1.fastq in2=./hisat2_output/test/preprocessing/kneaddata.trimmed.2.fastq out1=./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz out2=./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz outs=./hisat2_output/test/preprocessing/kneaddata_repaired_singletons.fastq.gz overwrite=t
Executing jgi.SplitPairsAndSingles [rp, in1=./hisat2_output/test/preprocessing/kneaddata.trimmed.1.fastq, in2=./hisat2_output/test/preprocessing/kneaddata.trimmed.2.fastq, out1=./hisat2_output/test/preprocessing/kneaddata_repaired_1.fastq.gz, out2=./hisat2_output/test/preprocessing/kneaddata_repaired_2.fastq.gz, outs=./hisat2_output/test/preprocessing/kneaddata_repaired_singletons.fastq.gz, overwrite=t]

Set INTERLEAVED to false
Started output stream.

Input:                  	200000 reads 		29641443 bases.
Result:                 	200000 reads (100.00%) 	29641443 bases (100.00%)
Pairs:                  	200000 reads (100.00%) 	29641443 bases (100.00%)
Singletons:             	0 reads (0.00%) 	0 bases (0.00%)

Time:                         	1.192 seconds.
Reads Processed:        200k 	167.82k reads/sec
Bases Processed:      29641k 	24.87m bases/sec
