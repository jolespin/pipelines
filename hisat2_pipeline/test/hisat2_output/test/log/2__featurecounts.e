Reading sequences [test.fa]: 0it [00:00, ?it/s]Reading sequences [test.fa]: 185it [00:00, 7476.98it/s]

        ==========     _____ _    _ ____  _____  ______          _____  
        =====         / ____| |  | |  _ \|  __ \|  ____|   /\   |  __ \ 
          =====      | (___ | |  | | |_) | |__) | |__     /  \  | |  | |
            ====      \___ \| |  | |  _ <|  _  /|  __|   / /\ \ | |  | |
              ====    ____) | |__| | |_) | | \ \| |____ / ____ \| |__| |
        ==========   |_____/ \____/|____/|_|  \_\______/_/    \_\_____/
	  v2.0.1

//========================== featureCounts setting ===========================\\
||                                                                            ||
||             Input files : 1 BAM file                                       ||
||                           o mapped.sorted.bam                              ||
||                                                                            ||
||             Output file : featurecounts.tsv                                ||
||                 Summary : featurecounts.tsv.summary                        ||
||              Annotation : genome.saf (SAF)                                 ||
||      Dir for temp files : ./hisat2_output/test/tmp                         ||
||                                                                            ||
||                 Threads : 4                                                ||
||                   Level : meta-feature level                               ||
||              Paired-end : no                                               ||
||      Multimapping reads : not counted                                      ||
|| Multi-overlapping reads : not counted                                      ||
||   Min overlapping bases : 1                                                ||
||                                                                            ||
\\============================================================================//

//================================= Running ==================================\\
||                                                                            ||
|| Load annotation file genome.saf ...                                        ||
||    Features : 185                                                          ||
||    Meta-features : 185                                                     ||
||    Chromosomes/contigs : 185                                               ||
||                                                                            ||
|| Load FASTA contigs from test.fa...                                         ||
||    184 contigs were loaded                                                 ||
||                                                                            ||
|| Process BAM file mapped.sorted.bam...                                      ||
||    WARNING: Paired-end reads were found.                                   ||
||    Total alignments : 200072                                               ||
||    Successfully assigned alignments : 7195 (3.6%)                          ||
||    Running time : 0.01 minutes                                             ||
||                                                                            ||
|| Write the final count table.                                               ||
|| Write the read assignment summary.                                         ||
||                                                                            ||
|| Summary of counting results can be found in file "./hisat2_output/test/in  ||
|| termediate/featurecounts_output/featurecounts.tsv.summary"                 ||
||                                                                            ||
\\============================================================================//

