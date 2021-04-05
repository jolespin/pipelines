#!/usr/bin/env python
from __future__ import print_function, division
import sys, os, argparse, glob
import pandas as pd
from collections import OrderedDict

# Soothsayer Ecosystem
from soothsayer_utils import get_file_object, pv, assert_acceptable_arguments, read_dataframe

pd.options.display.max_colwidth = 100
# from tqdm import tqdm
__program__ = os.path.split(sys.argv[0])[-1]
__version__ = "2021.04.04"




#
def main(args=None):
    # Path info
    script_directory  =  os.path.dirname(os.path.abspath( __file__ ))
    script_filename = __program__
    # Path info
    description = """
    Running: {} v{} via Python v{} | {}""".format(__program__, __version__, sys.version.split(" ")[0], sys.executable)
    usage = "{} --feature_axis columns /path/*/featurecounts.tsv.gz".format(__program__)
    epilog = "Copyright 2021 Josh L. Espinoza (jespinoz@jcvi.org)"

    # Parser
    parser = argparse.ArgumentParser(description=description, usage=usage, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
    # Pipeline
    parser.add_argument("-i","--project_directory", type=str, required=True, help = "Project directory for bowtie2 pipeline")
    parser.add_argument("-f","--feature_axis", type=str, default="columns", help = "Axis for features {rows, columns} [Default: columns")
    parser.add_argument("--summary", action="store_true", help = "Use to output summary stats instead of feature counts")



    # Options
    opts = parser.parse_args()
    opts.script_directory  = script_directory
    opts.script_filename = script_filename

    # Outputs
    if opts.summary:
        filename = "featurecounts.tsv.summary"
        feature_label = "MappingGroup"
    else:
        filename = "featurecounts.tsv.gz"
        feature_label = "FeatureID"

    
    merged_data = OrderedDict()
    for filepath in pv(glob.glob(os.path.join(opts.project_directory, "*", "intermediate", "featurecounts_output", filename)), "Reading featureCounts tables"):
        id_sample = filepath.split("/")[-4]
        merged_data[id_sample] = read_dataframe(filepath, skiprows=1, sep="\t").iloc[:,-1]
    df_featurecounts = pd.DataFrame(merged_data)
    df_featurecounts.index.name = feature_label
    df_featurecounts.columns.name = "SampleID"

    if opts.feature_axis == "columns":
        df_featurecounts = df_featurecounts.T 
    
    df_featurecounts.to_csv(sys.stdout, sep="\t")
        

if __name__ == "__main__":
    main()
