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
__version__ = "2021.04.04"

# .............................................................................
# Notes
# .............................................................................
# * Still needs unpaired
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
    {True:"--reference-db {}".format(os.environ["kneaddata_contamination_db"]), False:""}[bool(os.environ["kneaddata_contamination_db"])],
    "--output {}".format(output_directory),
    "--log {}".format(os.path.join(output_directory, "kneaddata.log")),
    "--threads {}".format(opts.n_jobs),
    "--output-prefix kneaddata",
    '--bowtie2-options="--seed {}"'.format(opts.random_state, {True:" {}".format(opts.kneaddatabowtie2_options), False:""}[bool(opts.kneaddatabowtie2_options)]), # Work around to ensure there's not an extra space for the options
    opts.kneaddata_options,
    "&&",
    os.environ["repair.sh"],
    "in1={}".format(os.path.join(output_directory, "kneaddata_paired_1.fastq")),
    "in2={}".format(os.path.join(output_directory, "kneaddata_paired_2.fastq")),
    "out1={}".format(os.path.join(output_directory,"kneaddata_repaired_1.fastq.gz")),
    "out2={}".format(os.path.join(output_directory, "kneaddata_repaired_2.fastq.gz")),
    "outs={}".format(os.path.join(output_directory, "kneaddata_repaired_singletons.fastq.gz")),
    "overwrite=t",
    "&&",
    "rm {} {}".format(
    os.path.join(output_directory, "kneaddata_paired*"),
    os.path.join(output_directory, "*trimmed*"),
    ),
    "&&",
    "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "*.fastq")),
    ")",
    # clumpify.sh ?
    ]
    return cmd


# Bowtie2
def get_bowtie2_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
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

    # Bowtie2
    "(",
    os.environ["bowtie2"],
   "-x {}".format(opts.bowtie2_index),
    ]
    if len(input_filepaths) == 1:
        cmd.append("--interleaved {}".format(input_filepaths[0]))
    else:
        cmd.append("-1 {}".format(input_filepaths[0]))
        cmd.append("-2 {}".format(input_filepaths[1]))

    # Do something with unpaired reads eventually

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
        opts.bowtie2_options,
    ")",
    # "|",
    # # Convert to BAM
    # os.environ["samtools"],
    # "view -h -b",
    # "|",
    # Convert to sorted BAM
    "|",
    "(",
    os.environ["reformat.sh"],
    "in=stdin.sam",
    "out=stdout.bam",
    "mappedonly=t",
    ")",

    # "sort",
    # "--threads {}".format(opts.n_jobs),
    # "--reference {}".format(opts.ref_assembly),
    # "--input-fmt SAM",
    # "--output-fmt BAM",
    # "-T {}".format(os.path.join(directories["tmp"], "samtools_sort")),
    ">",
    output_filepaths[0],
    # "&&",
    # "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "*{}".format(ext))), # May use need to use regular gzip if this throws an error
    ]
    return cmd



# featureCounts
def get_featurecounts_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
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
                "reformat.sh",
                # "clumpify.sh",
                "kneaddata",
                "kneaddata_contamination_db",
                # 2
                "bowtie2",
                "samtools",
                # 3
                "featureCounts",
     }


    if opts.path_config is None:
        opts.path_config = os.path.join(opts.script_directory, "bowtie2_pipeline_config.tsv")
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
    print(format_header( "Adding executables to path from config (Overriding the conda executables): {}".format(opts.path_config), "-"), file=sys.stdout)
    for name, executable in executables.items():
        if name in required_executables:
            print(name, executable, sep = " --> ", file=sys.stdout)
            os.environ[name] = executable.strip()


    # Expand real path of taxa.sqlite for ete[2/3]
    # os.environ["taxa.sqlite"] = os.path.expanduser(os.environ["taxa.sqlite"])
    print("", file=sys.stdout)
# Configure parameters
def configure_parameters(opts, directories):
    # os.environ[]
    # R1 != R2
    assert opts.r1 != opts.r2, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.r1)
    # Set environment variables
    add_executables_to_environment(opts=opts)




def create_pipeline(opts, directories, f_cmds):

    # .................................................................
    # Primordial
    # .................................................................
    # Commands file
    pipeline = ExecutablePipeline(name="Bowtie2 Mapping Pipeline", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

    if not opts.skip_preprocess:
        # =========
        # Kneaddata
        # =========
        program = "kneaddata"
        # Add to directories
        output_directory = directories["preprocessing"] = create_directory(os.path.join(directories["sample"], "preprocessing"))

        # Info
        step = 1
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
        if opts.interleaved_reads:
            output_filepaths = [opts.interleaved_reads]
        else:
            output_filepaths = [opts.r1, opts.r2]


    if not opts.preprocess_only:
        # ==========
        # Bowtie2
        # ==========
        program = "bowtie2"
        # Add to directories
        output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

        # Info
        step = 2
        description = "Aligning reads to reference"

        # i/o
        input_filepaths = output_filepaths
        assert len(input_filepaths) in {1,2}, "`input_filepaths` at this stage must have 1 file for interleaved or 2 files for paired:\n{}".format(input_filepaths)

        output_filenames = ["mapped.bam"]
        output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

        params = {
            "input_filepaths":input_filepaths,
            "output_filepaths":output_filepaths,
            "output_directory":output_directory,
            "opts":opts,
            "directories":directories,
        }

        cmd = get_bowtie2_cmd(**params)
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
        step = 3
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
    step = 4
    description = "Symlinking relevant output files"


    # i/o
    input_filepaths = [
        # os.path.join(directories[("intermediate", "bowtie2")], "mapped.sorted.bam"),
        # os.path.join(directories[("intermediate", "bowtie2")], "*"),
        os.path.join(directories[("intermediate", "featurecounts")], "featurecounts.tsv.gz"),
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
        assert not bool(opts.interleaved_reads), "Cannot have --interleaved_reads if --r1"
    if not opts.bowtie2_index:
        opts.bowtie2_index = opts.ref_assembly
    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: bowtie2_pipeline.py v{} via Python v{} | {}""".format(__version__, sys.version.split(" ")[0], sys.executable)
    usage = "bowtie2_pipeline.py -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --bowtie2_index <bowtie2_index/>"
    epilog = "Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required arguments')
    parser_required.add_argument("-1","--r1", type=str, help = "path/to/r1.fq")
    parser_required.add_argument("-2","--r2", type=str, help = "path/to/r2.fq")
    parser_required.add_argument("-12","--interleaved_reads", type=str, help = "path/to/interleaved.fq")
    # parser_required.add_argument("-U","--unpaired_reads", type=str, help = "path/to/unpaired_reads.fq. Can be comma separated list")
    # parser_required.add_argument("--include_unpaired_reads",  action="store_true", help = "Include unpaired reads")


    parser_required.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_required.add_argument("-o","--project_directory", type=str, default="./mapping_output", help = "path/to/project_directory [Default: ./mapping_output]")
    parser_required.add_argument("--ref_assembly", type=str, required=True, help = "path/to/reference.fasta" ) # ; or (2) a directory of fasta files [Must all have the same extension.  Use `query_ext` argument]
    parser_required.add_argument("--ref_annotation",type=str, required=True, help="path/to/reference.gtf")
    parser_required.add_argument("--bowtie2_index",type=str, required=False, help="path/to/bowtie2_index")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  help="path/to/bowtie2_pipeline_config.tsv\n[Default: {}]".format(os.path.join(script_directory, "bowtie2_pipeline_config.tsv")))  #site-packges in future
    parser_utility.add_argument("--n_jobs", type=int, default=1, help = "Number of threads [Default: 1]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--preprocess_only", action="store_true", help = "Only run preprocess steps (kneaddata)") #!
    parser_utility.add_argument("--skip_preprocess", action="store_true", help = "Don't run preprocess steps (kneaddata)") #!
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Kneaddata
    parser_kneaddata = parser.add_argument_group('Kneaddata arguments')
    parser_kneaddata.add_argument("--kneaddata_options", type=str, default="", help="Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home")
    parser_kneaddata.add_argument("--kneaddatabowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # Bowtie2
    parser_bowtie2 = parser.add_argument_group('Bowtie2 arguments')
    parser_bowtie2.add_argument("--bowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: ''] | http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
    parser_featurecounts.add_argument("--attribute_type", type=str, default="gene_id", help = "Attribute type in GTF/GFF file. Use 'ID' for prodigal. [Default: gene_id]")
    parser_featurecounts.add_argument("--feature_type", type=str, default="exon", help = "Feature type in GTF/GFF file. Use 'CDS' for prodigal. [Default: exon]")
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
