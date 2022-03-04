```
 _______ _______  ______ _____  _____  _______ _______
 |______ |       |_____/   |   |_____]    |    |______
 ______| |_____  |    \_ __|__ |          |    ______|

```
### Contents:
Collection of bioinformatics that can be used to speed up your workflow.

* `scaffolds_to_bins.py`
	- Generate a scaffolds to bins table for programs such as [DAS Tool](https://github.com/cmks/DAS_Tool).
	- Dependencies: `pandas`
	
* `fastani_to_clusters.py`
	- Cluster [FastANI](https://github.com/ParBLiSS/FastANI) output using [NetworkX](https://github.com/networkx/networkx).
	- Dependencies: `pandas`, `networkx`

* `fasta_to_saf.py`
	- Generate a [Simplified Annotation Format](https://rdrr.io/bioc/Rsubread/man/featureCounts.html) [SAF] table from Fasta for use with programs such as [featureCounts](http://subread.sourceforge.net/).
	- Dependencies: `pandas`, `BioPython`, [`soothsayer_utils`](https://github.com/jolespin/soothsayer_utils)

* `append_geneid_to_prodigal_gff.py`
	- Appends `gene_id` to [Prodigal](https://github.com/hyattpd/Prodigal) GFF output.  Useful when combining results from multiple Prodigal runs (e.g. metagenomic assembled genomes). 
	- Dependencies: None

* `convert_metaeuk.py`
	- Simplify the identifiers for [MetaEuk](https://github.com/soedinglab/metaeuk). Input is `--cds` and `--protein` then it generates: 1) gene\_models.ffn; 2) gene\_models.faa; 3) gene\_models.gff; and 4) identifier\_mapping.tsv
	- Dependencies: `pandas`, `BioPython`