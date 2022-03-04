#!/usr/bin/env python
import sys, os, glob, argparse 
from itertools import groupby 
from collections import OrderedDict
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser

__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2022.03.03"

def simplify_metaeuk_identifiers(identifiers):
    def f_sort(id):
        fields = id.split("|")
        start = int(fields[6])
        end = int(fields[7])
        return (start, end)
    
    metaeuk_to_simple = dict()
    for id_contig, contig_orf_ids in groupby(identifiers, key=lambda x: x.split("|")[1]):
        metaeuk_identifiers = sorted(contig_orf_ids, key=f_sort)
        simple_identifiers = map(lambda i: "{}_{}".format(id_contig, i), range(1, len(metaeuk_identifiers)+1))
        metaeuk_to_simple.update(dict(zip(metaeuk_identifiers, simple_identifiers)))
    return pd.Series(metaeuk_to_simple)

def main(argv=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -d <output.codon.fsa> -a <output.fsa> -o <output_directory>".format(__program__)
    epilog = "Copyright 2022 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-d","--cds", type=str, required=True, help = "path/to/output.codons.fsa")
    parser.add_argument("-a","--protein", type=str, required=True, help = "path/to/output.fsa")
    parser.add_argument("-o","--output_directory", type=str, default="metaeuk_output", help = "path/to/output_directory [Default: metaeuk_output]")
    parser.add_argument("-b", "--basename", type=str, default="gene_models", help = "Base name for gene models [Default: gene_models]")

    parser.add_argument("--no_header", action="store_true", help="Specify if header should not be in fasta files")
    parser.add_argument("--cds_extension", type=str, default="ffn", help = "CDS fasta extension [Default: ffn]")
    parser.add_argument("--protein_extension", type=str, default="faa", help = "Protein fasta extension [Default: faa]")

    # Options
    opts = parser.parse_args(argv)
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Output directory
    os.makedirs(opts.output_directory, exist_ok=True)

    # Load in CDS sequences
    id_to_seq = OrderedDict()
    with open(opts.cds, "r") as f:
        for header, seq in SimpleFastaParser(f):
            id = header.split(" ")[0]
            id_to_seq[id] = seq
    cds_sequences = pd.Series(id_to_seq)

    # Load in protein sequences
    id_to_seq = OrderedDict()
    with open(opts.protein, "r") as f:
        for header, seq in SimpleFastaParser(f):
            id = header.split(" ")[0]
            id_to_seq[id] = seq
    protein_sequences = pd.Series(id_to_seq)

    # Convert MetaEuk to simple identifiers
    metaeuk_to_simple = simplify_metaeuk_identifiers(cds_sequences.index)

    # Relabel identifiers
    if opts.no_header:
        cds_sequences.index = cds_sequences.index.map(lambda id_metaeuk: metaeuk_to_simple[id_metaeuk])
        protein_sequences.index = protein_sequences.index.map(lambda id_metaeuk: metaeuk_to_simple[id_metaeuk])
    else:
        cds_sequences.index = cds_sequences.index.map(lambda id_metaeuk: "{} {}".format(metaeuk_to_simple[id_metaeuk], id_metaeuk))
        protein_sequences.index = protein_sequences.index.map(lambda id_metaeuk: "{} {}".format(metaeuk_to_simple[id_metaeuk], id_metaeuk))


    gff_output = list()
    for id_metaeuk, id_gene in metaeuk_to_simple.items():
        metaeuk_fields = id_metaeuk.split("|")
        
        # Gene
        id_target = metaeuk_fields[0]
        id_contig = metaeuk_fields[1]
        strand = metaeuk_fields[2]
        id_tcs = "|".join([id_target, id_contig, strand])

        start_gene = int(metaeuk_fields[6]) + 1
        end_gene = int(metaeuk_fields[7]) + 1
        bitscore = float(metaeuk_fields[3])
        gene_description = "target_id={};tcs_id={};contig_id={};gene_id={};".format(
            id_target,
            id_tcs,
            id_contig,
            id_gene,
        )
        gene_fields = [id_contig, "MetaEuk", "gene", start_gene, end_gene, bitscore, strand, ".", gene_description]

        # gff_output.append("\t".join(map(str, gene_fields)))
        gff_output.append(gene_fields)

        # mRNA
        mrna_description = "target_id={};tcs_id={};contig_id={};gene_id={};Parent={};".format(
            id_target,
            id_tcs,
            id_contig,
            id_gene,
            id_gene,
        )
        mrna_fields = [id_contig, "MetaEuk", "mRNA", start_gene, end_gene, bitscore, strand, ".", mrna_description]

        # gff_output.append("\t".join(map(str, mrna_fields)))
        gff_output.append(mrna_fields)

        # Exons
        for i, exon in enumerate(metaeuk_fields[8:], start=1):
            exon_fields = exon.split(":")
            start_exon = int(exon_fields[0].split("[")[0]) + 1
            end_exon = int(exon_fields[1].split("[")[0]) + 1
            id_exon = "{}.{}".format(id_gene, i)

            exon_description = "target_id={};tcs_id={};contig_id={};gene_id={};Parent={};exon_id={}".format(
                id_target,
                id_tcs,
                id_contig,
                id_gene,
                id_gene,
                id_exon,
            )
            exon_fields = [id_contig, "MetaEuk", "exon", start_exon, end_exon, bitscore, strand, ".", exon_description]
            gff_output.append(exon_fields)
            # output.append("\t".join(map(str, exon_fields)))

    df_gff = pd.DataFrame(gff_output)

    # Write output

    # Identifiers
    metaeuk_to_simple.to_frame().to_csv(os.path.join(opts.output_directory,"identifier_mapping.tsv"), sep="\t", header=False)

    # CDS
    with open(os.path.join(opts.output_directory,"{}.{}".format(opts.basename, opts.cds_extension)) , "w") as f:
        f.writelines(">{}\n{}\n".format(id, seq) for id, seq in cds_sequences.items())

    # Protein
    with open(os.path.join(opts.output_directory,"{}.{}".format(opts.basename, opts.protein_extension)) , "w") as f:
        f.writelines(">{}\n{}\n".format(id, seq) for id, seq in protein_sequences.items())
    # GFF
    df_gff.to_csv(os.path.join(opts.output_directory,"{}.gff".format(opts.basename)), sep="\t", index=False, header=False)
    
if __name__ == "__main__":
    main()
    
                

