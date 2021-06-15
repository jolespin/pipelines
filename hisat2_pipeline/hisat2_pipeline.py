#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, subprocess, datetime, time, glob, shutil, pickle
from collections import OrderedDict, defaultdict
# Scandir
try:
    from os import scandir
except ImportError:
    from scandir import scandir
# Pathlib
try:
    import pathlib
except ImportError:
    import pathlib2 as pathlib

import pandas as pd
import numpy as np

# Soothsayer Ecosystem
from genopype import *
from soothsayer_utils import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.06.15"

# .............................................................................
# Notes
# .............................................................................
# * Make batch version that takes in a manifest file
# .............................................................................
# Primordial
# .............................................................................

# Kneaddata
def get_kneaddata_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    "(",
    os.environ["kneaddata"],
    "--input {}".format(input_filepaths[0]),
    "--input {}".format(input_filepaths[1]),
    {True:"--reference-db {}".format(opts.kneaddata_contamination_db), False:""}[bool(opts.kneaddata_contamination_db)],
    "--output {}".format(output_directory),
    "--log {}".format(os.path.join(output_directory, "kneaddata.log")),
    "--threads {}".format(opts.n_jobs),
    "--output-prefix kneaddata",
    '--bowtie2-options="--seed {}"'.format(opts.random_state, {True:" {}".format(opts.kneaddatabowtie2_options), False:""}[bool(opts.kneaddatabowtie2_options)]), # Work around to ensure there's not an extra space for the options
    opts.kneaddata_options,
    ")",
    "&&",
    "(",
    os.environ["repair.sh"],
    ]
    if bool(opts.kneaddata_contamination_db):
        cmd += [
        "in1={}".format(os.path.join(output_directory, "kneaddata_paired_1.fastq")),
        "in2={}".format(os.path.join(output_directory, "kneaddata_paired_2.fastq")),
        ]
    else:
        cmd += [
        "in1={}".format(os.path.join(output_directory, "kneaddata.trimmed.1.fastq")),
        "in2={}".format(os.path.join(output_directory, "kneaddata.trimmed.2.fastq")),
        ]

    cmd += [
        "out1={}".format(os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")),
        "out2={}".format(os.path.join(output_directory, "kneaddata_repaired_2.fastq.gz")),
        "outs={}".format(os.path.join(output_directory, "kneaddata_repaired_singletons.fastq.gz")),
        "overwrite=t",
        ")",
        "&&",
        "(",
        "rm -f {} {}".format(
            os.path.join(output_directory, "kneaddata_paired*"),
            os.path.join(output_directory, "*trimmed*"),
            ),
        # "&&",
        # "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "*.fastq")),
        ")",
    ]


    return cmd


# HISAT2
def get_hisat2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Is fasta or fastq?
    ignore_quals = False
    unmapped_ext = "fastq"
    fp, ext = os.path.splitext(input_filepaths[0])
    if ext in {".gz",".bz2"}:
        fp, ext = os.path.splitext(fp)
        ext = ext[1:]
    if ext in {"fa", "fasta", "fna"}:
        unmapped_ext = "fasta"
        ignore_quals = True
    
    # Command
    cmd = [
    # Clear temporary directory just in case
    "rm -rf {}".format(os.path.join(directories["tmp"], "*")),
    "&&",
    # HISAT2
    "(",
    os.environ["hisat2"],
   "-x {}".format(opts.hisat2_index),
    ]
    if len(input_filepaths) == 1:
        cmd.append("-U {}".format(input_filepaths[0]))
    else:
        cmd.append("-1 {}".format(input_filepaths[0]))
        cmd.append("-2 {}".format(input_filepaths[1]))


    cmd += [
   "--threads {}".format(opts.n_jobs),
   # Do something with reads that don't map
    "--un-conc-gz {}".format(os.path.join(output_directory, "unmapped_%.{}.gz".format(unmapped_ext))),
    # "--un-gz {}".format(os.path.join(output_directory, "unmapped_singletons_%{}.gz".format(unmapped_ext))),
    "--seed {}".format(opts.random_state),
    ]
    if ignore_quals:
        cmd.append("-f")
    cmd += [
        opts.hisat2_options,
    ")",
    # Convert to sorted BAM
    "|",
    "(",
    os.environ["samtools"],
    "sort",
    "--threads {}".format(opts.n_jobs),
    "--reference {}".format(opts.ref_assembly),
    "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    ">",
    output_filepaths[0],
    ")",
    ]
    return cmd



# featureCounts
def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    if opts.ref_annotation:

        # Command
        cmd = [
            "(",
            os.environ["featureCounts"],
        "-G {}".format(opts.ref_assembly),
        "-a {}".format(opts.ref_annotation),
        "-o {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        "-F GTF",
        "--tmpDir {}".format(directories["tmp"]),
        "-T {}".format(opts.n_jobs),
        "-g {}".format(opts.attribute_type),
        "-t {}".format(opts.feature_type),
            opts.featurecounts_options,
            input_filepaths[0],
            ")",
            "&&",
            "gzip {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        ]
    else:
        cmd = [
        "("
        "python",
        os.path.join(opts.script_directory, "fasta_to_saf.py"),
        "-i {}".format(opts.ref_assembly),
        ">",
        os.path.join(output_directory, "genome.saf"),
        ")",
        "&&",
        "(",
        os.environ["featureCounts"],
        "-G {}".format(opts.ref_assembly),
        "-a {}".format(os.path.join(output_directory, "genome.saf")),
        "-o {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        "-F SAF",
        "--tmpDir {}".format(directories["tmp"]),
        "-T {}".format(opts.n_jobs),
        opts.featurecounts_options,
        input_filepaths[0],
        ")",
        "&&",
        "gzip {}".format(os.path.join(output_directory, "featurecounts.tsv")),
        ]
    return cmd

# Symlink
def get_symlink_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    # Command
    cmd = ["("]
    for filepath in input_filepaths:
        cmd.append("ln -f -s {} {}".format(os.path.realpath(filepath), os.path.realpath(output_directory)))
        cmd.append("&&")
    cmd[-1] = ")"
    return cmd

# ============
# Run Pipeline
# ============
# Set environment variables
def add_executables_to_environment(opts):
    """
    Adapted from Soothsayer: https://github.com/jolespin/soothsayer
    """

    required_executables={
                # 1
                "repair.sh",
                # "clumpify.sh",
                "kneaddata",
                # "kneaddata_contamination_db",
                # 2
                "hisat2",
                "samtools",
                # 3
                "featureCounts",
     }

    if opts.path_config == "CONDA_PREFIX":
        executables = dict()
        for name in required_executables:
            executables[name] = os.path.join(os.environ["CONDA_PREFIX"], "bin", name)
    else:
        if opts.path_config is None:
            opts.path_config = os.path.join(opts.script_directory, "hisat2_pipeline_config.tsv")
        opts.path_config = format_path(opts.path_config)
        assert os.path.exists(opts.path_config), "config file does not exist.  Have you created one in the following directory?\n{}\nIf not, either create one, check this filepath:{}, or give the path to a proper config file using --path_config".format(opts.script_directory, opts.path_config)
        assert os.stat(opts.path_config).st_size > 1, "config file seems to be empty.  Please add 'name' and 'executable' columns for the following program names: {}".format(required_executables)
        df_config = pd.read_csv(opts.path_config, sep="\t")
        assert {"name", "executable"} <= set(df_config.columns), "config must have `name` and `executable` columns.  Please adjust file: {}".format(opts.path_config)
        df_config = df_config.loc[:,["name", "executable"]].dropna(how="any", axis=0).applymap(str)
        # Get executable paths
        executables = OrderedDict(zip(df_config["name"], df_config["executable"]))
        assert required_executables <= set(list(executables.keys())), "config must have the required executables for this run.  Please adjust file: {}\nIn particular, add info for the following: {}".format(opts.path_config, required_executables - set(list(executables.keys())))

    # Display
    accessory_scripts = ["fasta_to_saf.py"]
    for name in accessory_scripts:
        executables[name] = "python " + os.path.join(opts.script_directory, name)
    print(format_header( "Adding executables to path from the following source: {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()
    print("", file=sys.stdout)





def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name="HISAT2 Mapping Pipeline", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    if not opts.skip_preprocess:
        # =========
        # Kneaddata
        # =========
        program = "kneaddata"
        # Add to directories
        output_directory = directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))

        # Info
        step = 0
        description = "Quality trimming, removing contaminated reads, and optimizing file compression"
        # i/o
        # if not opts.bypass_decontamination:
        input_filepaths = [opts.r1, opts.r2]
        # else:
            # assert opts.single_reads is not None, "If `bypass_decontamination` then the reads must be single-ended"
            # input_filepaths = [os.path.join(output_directory, "reads.subsampled.clumpify.fastq.gz")]
            # os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")
        output_filenames = ["kneaddata_repaired_1.fastq.gz", "kneaddata_repaired_2.fastq.gz"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        # Parameters
        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }
        # Command
        cmd = get_kneaddata_cmd(**params)
        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
        )
    else:
        if opts.unpaired_reads:
            output_filepaths = [opts.unpaired_reads]
        else:
            output_filepaths = [opts.r1, opts.r2]


    if not opts.preprocess_only:
        # ==========
        # HISAT2
        # ==========
        program = "hisat2"
        # Add to directories
        output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

        # Info
        step = 1
        description = "Aligning reads to reference"

        # i/o
        input_filepaths = output_filepaths
        assert len(input_filepaths) in {1,2}, "`input_filepaths` at this stage must have 1 file for unpaired or 2 files for paired:\n{}".format(input_filepaths)

        output_filenames = ["mapped.sorted.bam"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_hisat2_cmd(**params)
        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
    )

        # ==========
        # featureCounts
        # ==========
        program = "featurecounts"
        # Add to directories
        output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

        # Info
        step = 2
        description = "Counting reads"

        # i/o
        input_filepaths = output_filepaths

        output_filenames = ["featurecounts.tsv.gz"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_featurecounts_cmd(**params)
        pipeline.add_step(
                    id=program,
                    description = description,
                    step=step,
                    cmd=cmd,
                    input_filepaths = input_filepaths,
                    output_filepaths = output_filepaths,
                    validate_inputs=True,
                    validate_outputs=True,
        )

        # =============
        # Symlink
        # =============
        program = "symlink"
        # Add to directories
        output_directory = directories["output"]

        # Info
        step = 3
        description = "Symlinking relevant output files"


        # i/o
        input_filepaths = [
            os.path.join(directories[("intermediate", "hisat2")], "mapped.sorted.bam"),
            # os.path.join(directories[("intermediate", "hisat2")], "*"),
            os.path.join(directories[("intermediate", "featurecounts")], "featurecounts.tsv.gz"),
            os.path.join(directories[("intermediate", "featurecounts")], "featurecounts.tsv.summary"),

        ]

        output_filenames =  map(lambda fp: fp.split("/")[-1], input_filepaths)
        output_filepaths = list(map(lambda fn:os.path.join(directories["output"], fn), output_filenames))
            # Prodigal
            # os.path.join(directories["output"], "*"),
        
        params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
        }

        cmd = get_symlink_cmd(**params)
        pipeline.add_step(
                id=program,
                description = description,
                step=step,
                cmd=cmd,
                input_filepaths = input_filepaths,
                output_filepaths = output_filepaths,
                validate_inputs=True,
                validate_outputs=False,
        )

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    # os.environ[]
    if opts.preprocess_only:
        assert not opts.skip_preprocess, "Conflicting logic with `--preprocess_only` and `--skip_preprocess`"
    if bool(opts.r1):
        assert opts.r1 != opts.r2, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.r1)
        assert not bool(opts.unpaired_reads), "Cannot have --unpaired_reads if --r1.  Note, this behavior may be changed in the future but it's an adaptation of interleaved reads."
    if not opts.hisat2_index:
        opts.hisat2_index = opts.ref_assembly
    if not opts.preprocess_only:
        assert opts.ref_assembly is not None, "Please provide --ref_assembly for mapping"
        if opts.ref_annotation is None:
            print("No --ref_annotation has been provided.  Converting FASTA -> SAF and aggregating counts per sequence record", file=sys.stdout) #assert opts.ref_annotation is not None, "Please provide --ref_annotation for feature counting"
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --hisat2_index <hisat2_index/>".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required arguments')
    parser_required.add_argument("-1","--r1", type=str, help = "path/to/r1.fq")
    parser_required.add_argument("-2","--r2", type=str, help = "path/to/r2.fq")
    parser_required.add_argument("-U","--unpaired_reads", type=str, help = "path/to/unpaired_reads.fq. Can be comma separated list")


    parser_required.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_required.add_argument("-o","--project_directory", type=str, default="./hisat2_output", help = "path/to/project_directory [Default: ./hisat2_output]")
    parser_required.add_argument("-R", "--ref_assembly", type=str, required=False, help = "path/to/reference.fasta" ) # ; or (2) a directory of fasta files [Must all have the same extension.  Use `query_ext` argument]
    parser_required.add_argument("-A", "--ref_annotation",type=str, required=False, help="path/to/reference.gtf")
    parser_required.add_argument("-I", "--hisat2_index",type=str, required=False, help="path/to/hisat2_index")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  default="CONDA_PREFIX", help="path/to/config.tsv [Default: CONDA_PREFIX]")  #site-packges in future
    parser_utility.add_argument("--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--preprocess_only", action="store_true", help = "Only run preprocess steps (kneaddata)") #!
    parser_utility.add_argument("--skip_preprocess", action="store_true", help = "Don't run preprocess steps (kneaddata)") #!
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Kneaddata
    parser_kneaddata = parser.add_argument_group('Kneaddata arguments')
    parser_kneaddata.add_argument("--kneaddata_contamination_db", type=str, default="", help="Kneaddata | path/to/contamination_database [Default: '']\nFor human at JCVI, use the following: /usr/local/scratch/METAGENOMICS/jespinoz/db/genomes/human/GRCh38.p13/")
    parser_kneaddata.add_argument("--kneaddata_options", type=str, default="", help="Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home")
    parser_kneaddata.add_argument("--kneaddatabowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # HISAT2
    parser_hisat2 = parser.add_argument_group('HISAT2 arguments')
    parser_hisat2.add_argument("--hisat2_options", type=str, default="", help="HISAT2 | More options (e.g. --arg 1 ) [Default: ''] | http://daehwankimlab.github.io/hisat2/")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    # parser_featurecounts.add_argument("--annotation_type", type=str, default="GTF", help = "Annotation file type. Use 'SAF' when there are no features (fasta_to_saf.py helper script converts fasta to saf). [Default: GTF]")
    parser_featurecounts.add_argument("-g", "--attribute_type", type=str, default="gene_id", help = "Attribute type in GTF/GFF file. Use 'ID' for prodigal. [Default: gene_id]")
    parser_featurecounts.add_argument("-t", "--feature_type", type=str, default="exon", help = "Feature type in GTF/GFF file. Use 'CDS' for prodigal. Use 'gene' for prokaryotic genomes from NCBI. [Default: exon]")
    parser_featurecounts.add_argument("--featurecounts_options", type=str, default="", help="featureCounts | More options (e.g. --arg 1 ) [Default: ''] | http://bioinf.wehi.edu.au/featureCounts/")

    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename


    # Potential names
    {}

    # Directories
    directories = dict()
    directories["project"] = create_directory(opts.project_directory)
    directories["sample"] = create_directory(os.path.join(directories["project"], opts.name))
    directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))
    directories["output"] = create_directory(os.path.join(directories["sample"], "output"))
    directories["log"] = create_directory(os.path.join(directories["sample"], "log"))
    directories["tmp"] = create_directory(os.path.join(directories["sample"], "tmp"))
    directories["checkpoints"] = create_directory(os.path.join(directories["sample"], "checkpoints"))
    directories["intermediate"] = create_directory(os.path.join(directories["sample"], "intermediate"))

    # Configure parameters


    # Info
    print(format_header("HISAT2 Pipeline", "="), file=sys.stdout)
    print(format_header("Configuration:", "-"), file=sys.stdout)
    print(format_header("Name: {}".format(opts.name), "."), file=sys.stdout)
    print("Python version:", sys.version.replace("\n"," "), file=sys.stdout)
    print("Python path:", sys.executable, file=sys.stdout) #sys.path[2]
    print("Script version:", __version__, file=sys.stdout)
    print("Moment:", get_timestamp(), file=sys.stdout)
    print("Directory:", os.getcwd(), file=sys.stdout)
    print("Commands:", list(filter(bool,sys.argv)),  sep="\n", file=sys.stdout)
    configure_parameters(opts, directories)
    sys.stdout.flush()

    # Run pipeline
    with open(os.path.join(directories["sample"], "commands.sh"), "w") as f_cmds:
        pipeline = create_pipeline(
                     opts=opts,
                     directories=directories,
                     f_cmds=f_cmds,
        )
        pipeline.compile()
        pipeline.execute(restart_from_checkpoint=opts.restart_from_checkpoint)

if __name__ == "__main__":
    main()
