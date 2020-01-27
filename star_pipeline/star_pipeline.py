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

# Goomba
from goomba import *

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2020.01.15"


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
    '--bowtie2-options="--seed {}"'.format(opts.random_state, {True:" {}".format(opts.bowtie2_options), False:""}[bool(opts.bowtie2_options)]), # Work around to ensure there's not an extra space for the options
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


# STAR
def get_star_cmd(input_filepaths, output_filepaths, output_directory, directories, opts):
    os.environ["TMPDIR"] = directories["tmp"]

    # Command
    cmd = [
    os.environ["STAR"],
   "--genomeDir {}".format(opts.star_index),
   "--readFilesCommand zcat",
   "--readFilesIn {} {}".format(input_filepaths[0], input_filepaths[1]),
   "--outFileNamePrefix {}/star_".format(output_directory),
   "--runThreadN {}".format(opts.n_jobs),
   # "--outFilterScoreMinOverLread 0",
   # "--outFilterMatchNminOverLread 0",
   # "--outFilterMatchNmin {}".format(CONFIG_PARAMETERS["outFilterMatchNmin"]),
   "--outReadsUnmapped Fastx",
    opts.star_options,
    "&&",
    os.environ["samtools"],
    "view -h -b {}".format(os.path.join(output_directory, "star_Aligned.out.sam")),
    ">",
    os.path.join(output_directory, "star_Aligned.out.bam"),

    # "|",
    # os.environ["samtools"],
    # "sort",
    # "-@ {}".format(opts.n_jobs),
    # ">",
    # os.path.join(output_directory, "star_Aligned.out.sorted.bam"),
    "&&",
    "rm {}".format(os.path.join(output_directory, "star_Aligned.out.sam")),
    "&&",
    "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "star_Unmapped.out.mate1")), # May use need to use regular gzip if this throws an error
    "&&",
    "pigz -f -p {} {}".format(opts.n_jobs, os.path.join(output_directory, "star_Unmapped.out.mate2")),
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
    opts.featurecounts_options,
    input_filepaths[0],
    "&&",
    "gzip {}".format(os.path.join(output_directory, "featurecounts.tsv")),
    ]
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
                "kneaddata_contamination_db",
                # 2
                "STAR",
                "samtools",
                # 3
                "featureCounts",
     }


    if opts.path_config is None:
        opts.path_config = os.path.join(opts.script_directory, "star_pipeline_config.tsv")
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
    pipeline = ExecutablePipeline(name="STAR Mapping Pipeline", description=opts.name, f_cmds=f_cmds, checkpoint_directory=directories["checkpoints"], log_directory=directories["log"])

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


    # ==========
    # STAR
    # ==========
    program = "star"
    # Add to directories
    output_directory = directories[("intermediate",  program)] = create_directory(os.path.join(directories["intermediate"], "{}_output".format(program)))

    # Info
    step = 2
    description = "Aligning reads to reference"

    # i/o
    input_filepaths = output_filepaths

    output_filenames = ["star_Aligned.out.bam"]
    output_filepaths = list(map(lambda filename: os.path.join(output_directory, filename), output_filenames))

    params = {
        "input_filepaths":input_filepaths,
        "output_filepaths":output_filepaths,
        "output_directory":output_directory,
        "opts":opts,
        "directories":directories,
    }

    cmd = get_star_cmd(**params)
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
    description = "Counting mapped reads"

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

    return pipeline

# Configure parameters
def configure_parameters(opts, directories):
    # os.environ[]
    # R1 != R2
    assert opts.r1 != opts.r2, "You probably mislabeled the input files because `r1` should not be the same as `r2`: {}".format(opts.r1)


    # Set environment variables
    add_executables_to_environment(opts=opts)

def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: star_pipeline.py v{} via Python v{} | {}""".format(__version__, sys.version.split(" ")[0], sys.executable)
    usage = "star_pipeline.py -1 <r1.fq> -2 <r2.fq> -n <name> -o <output_directory> --ref_assembly <reference.fa> --ref_annotation <reference.gtf> --star_index <star_index/>"
    epilog = "Copyright 2020 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser_required = parser.add_argument_group('Required arguments')
    parser_required.add_argument("-1","--r1", type=str, help = "path/to/r1.fq", required=True)
    parser_required.add_argument("-2","--r2", type=str, help = "path/to/r2.fq", required=True)
    parser_required.add_argument("-n", "--name", type=str, help="Name of sample", required=True)
    parser_required.add_argument("-o","--project_directory", type=str, required=True, default="./mapping_output", help = "path/to/project_directory [Default: ./mapping_output]")
    parser_required.add_argument("--ref_assembly", type=str, required=True, help = "path/to/reference.fasta" ) # ; or (2) a directory of fasta files [Must all have the same extension.  Use `query_ext` argument]
    parser_required.add_argument("--ref_annotation",type=str, required=True, help="path/to/reference.gtf")
    parser_required.add_argument("--star_index",type=str, required=True, help="path/to/STAR_index")

    # Utility
    parser_utility = parser.add_argument_group('Utility arguments')
    parser_utility.add_argument("--path_config", type=str,  help="path/to/star_pipeline_config.tsv\n[Default: {}]".format(os.path.join(script_directory, "star_pipeline_config.tsv")))  #site-packges in future
    parser_utility.add_argument("--n_jobs", type=int, default=4, help = "Number of threads [Default: 4]")
    parser_utility.add_argument("--random_state", type=int, default=0, help = "Random state [Default: 0]")
    parser_utility.add_argument("--preprocess_only", action="store_true", help = "Only run preprocess steps (kneaddata)") #!
    parser_utility.add_argument("--restart_from_checkpoint", type=str, default=None, help = "Restart from a particular checkpoint [Default: None]")
    parser_utility.add_argument("-v", "--version", action='version', version="{} v{}".format(__program__, __version__))

    # Kneaddata
    parser_kneaddata = parser.add_argument_group('Kneaddata arguments')
    parser_kneaddata.add_argument("--kneaddata_options", type=str, default="", help="Kneaddata | More options (e.g. --arg 1 ) [Default: ''] | https://bitbucket.org/biobakery/kneaddata/wiki/Home")
    parser_kneaddata.add_argument("--bowtie2_options", type=str, default="", help="Bowtie2 | More options (e.g. --arg 1 ) [Default: '']\nhttp://bowtie-bio.sourceforge.net/bowtie2/manual.shtml")

    # STAR
    parser_star = parser.add_argument_group('STAR arguments')
    parser_star.add_argument("--star_options", type=str, default="", help="STAR | More options (e.g. --arg 1 ) [Default: ''] | https://github.com/alexdobin/STAR")

    # featureCounts
    parser_featurecounts = parser.add_argument_group('featureCounts arguments')
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
    print(format_header("STAR Pipeline", "="), file=sys.stdout)
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
