#!/usr/bin/env python

import argparse
from functools import partial
from multiprocessing import Pool, cpu_count
from numpy import isnan
import os
import re
import shutil
import subprocess
import sys
# Python 3+
if sys.version_info > (3, 0):
    from urllib.request import urlretrieve
# Python 2.7
else:
    from urllib import urlretrieve

from __init__ import ENCODE, ParseUtils

usage_msg = """
usage: %s --genome STR --feature STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
retrieves data from ENCODE (Encyclopedia of DNA Elements)

  --genome STR        genome assembly
  --feature STR       type of genomic feature ("accessibility",
                      "histone" or "tf")

optional arguments:
  -h, --help          show this help message and exit
  -o, --out-dir DIR   output directory (default = "./")
  --threads INT       number of threads to use (default = %s)
""" % (usage_msg, (cpu_count() - 1))

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an
    {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--genome")
    parser.add_argument("--feature")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-o", "--out-dir", default=".")
    optional_group.add_argument("--threads", default=(cpu_count() - 1))

    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if not args.genome or not args.feature:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--genome\" \"--feature\" are required\n"]
        print(": ".join(error))
        exit(0)

    # Check for invalid feature
    valid_features = ["accessibility", "histone", "tf"]
    if args.feature not in valid_features:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"--feature\"", "invalid choice", "\"%s\" (choose from" % \
            args.feature, "%s)\n" % " ".join(["\"%s\"" % i for i in valid_features])]
        print(": ".join(error))
        exit(0)

    # Check "--threads" argument
    try:
        args.threads = int(args.threads)
    except:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "argument \"-t\" \"--threads\"", "invalid int value", "\"%s\"\n" % \
            args.threads]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Retrieve ENCODE data
    encode(args.genome, args.feature, args.out_dir, args.threads)

def encode(genome, feat_type, out_dir=".", threads=1):
    """
    e.g. ./encode.py --genome hg38 --feature accessibility
    """

    # Initialize
    global encodes

    # Create output dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)   

    # Download metadata
    metadata_file = _download_ENCODE_metadata(genome, feat_type, out_dir)

    # Unless pickle file exists
    pickle_file = "%s.pickle.gz" % metadata_file
    if not os.path.exists(pickle_file):

        # Parse metadata
        encodes = _parse_ENCODE_metadata(genome, metadata_file)
        # Save pickle
        ParseUtils.save_compressed_pickle(pickle_file, encodes)

    else:

        # Load pickle
        encodes = ParseUtils.load_compressed_pickle(pickle_file)

    # Filter ENCODE accessions (i.e. for each experiment, keep the accession of
    # the best type of file), then group ENCODE accessions by experiment target
    # and type
    grouped_accessions = _group_ENCODE_accessions(_filter_ENCODE_accessions(feat_type))

    # For each experiment target/type...
    for experiment_target, experiment_type in sorted(grouped_accessions):

        # Beware, for this should not be possible!
        if experiment_target is not None:
            if feat_type != "histone" and feat_type != "tf":
                continue
        else:
            if feat_type == "histone" or feat_type == "tf":
                continue

        # Prepare data
        subgrouped_accessions = grouped_accessions[(experiment_target, experiment_type)]
        data_file = _preprocess_data(subgrouped_accessions, out_dir, threads)

def _download_ENCODE_metadata(genome, feat_type, out_dir="."):

    # Initialize
    url = "https://www.encodeproject.org/metadata/type=Experiment&status=released"

    # Fix hg38
    if genome == "hg38":
        genome = "GRCh38"

    # Add experiment to URL
    if feat_type == "accessibility":
        experiment_url = "&assay_title=DNase-seq&assay_title=FAIRE-seq&assay_slims=DNA%2Baccessibility&files.file_type=bed%2BnarrowPeak"
    elif feat_type == "histone":
        experiment_url = "&assay_title=Histone%2BChIP-seq&assay_slims=DNA%2Bbinding&files.file_type=bed%2BnarrowPeak"
    else:
        # TF
        experiment_url = "&assay_title=TF%2BChIP-seq&assay_slims=DNA%2Bbinding&files.file_type=bed%2BnarrowPeak"
    url += experiment_url

    # Add genome assembly
    url += "&assembly=%s" % genome

    # Download metadata
    metadata_file = os.path.join(out_dir, "metadata.%s.%s.tsv" % (genome, feat_type))
    if not os.path.exists(metadata_file):
        urlretrieve(os.path.join(url, "metadata.tsv"), metadata_file)

    return(metadata_file)

def _parse_ENCODE_metadata(genome, metadata_file):

    # Initialize
    accession_idx = None
    encode_objects = {}
    regexp = re.compile("^(3xFLAG|eGFP)?-?(.+)-(human|mouse|dmelanogaster)$")

    # For each line...
    for line in ParseUtils.parse_tsv_file(metadata_file):

        # If first line...
        if accession_idx is not None:

            # Initialize
            accession = line[accession_idx]
            biosample_name = line[biosample_name_idx]
            biosample_type = line[biosample_type_idx]
            download_url = line[download_idx]
            experiment_accession = line[experiment_acc_idx]
            experiment_type = line[experiment_type_idx]
            experiment_target = line[experiment_target_idx]
            if type(experiment_target) is float and isnan(experiment_target):
                experiment_target = None
            else:
                m = regexp.search(experiment_target)
                experiment_target = m.group(2)
            genome_assembly = line[genome_assembly_idx]
            output_format = line[output_format_idx]
            output_type = line[output_type_idx]
            status = line[status_idx]
            treatments = line[treatment_idx]
            if type(treatments) is float and isnan(treatments):
                treatments = None
            genetic_modifications = line[genetic_modifications_idx]
            if type(genetic_modifications) is float and isnan(genetic_modifications):
                genetic_modifications = None

            # Add ENCODE object
            encode = ENCODE(accession, biosample_name, biosample_type, download_url,
                experiment_accession, experiment_type, experiment_target, genome_assembly,
                output_format, output_type, status, treatments, genetic_modifications)
            if (encode.genome_assembly == genome and encode.status == "released"):
                encode_objects.setdefault(encode.accession, encode)

        else:
            # Get indices
            accession_idx = line.index("File accession")
            biosample_name_idx = line.index("Biosample term name")
            biosample_type_idx = line.index("Biosample type")
            download_idx = line.index("File download URL")
            experiment_acc_idx = line.index("Experiment accession")
            experiment_type_idx = line.index("Assay")
            experiment_target_idx = line.index("Experiment target")
            genome_assembly_idx = line.index("File assembly")
            output_format_idx = line.index("File format")
            output_type_idx = line.index("Output type")
            status_idx = line.index("File Status") - 1
            treatment_idx = line.index("Biosample treatments")
            genetic_modifications_idx = line.index("Biosample genetic modifications methods")

    return(encode_objects)

def _filter_ENCODE_accessions(feat_type):

    # Initialize
    done = set()
    grouped_accessions = {}
    filtered_accessions = set()

    # Possible output files (from ENCODE pipelines)
    # https://www.encodeproject.org/atac-seq/
    # https://www.encodeproject.org/data-standards/dnase-seq/
    # https://www.encodeproject.org/chip-seq/histone/
    # https://www.encodeproject.org/chip-seq/transcription_factor/
    output_types = {
        "accessibility": ["peaks"],
        "atac-seq": [],
        "histone": ["peaks"],
        "tf": ["optimal idr thresholded peaks", "conservative idr thresholded peaks",
            "pseudoreplicated idr thresholded peaks", "peaks"]
    }

    # For each accession...
    for accession in encodes:

        # Group ENCODE objects by experiment accession
        encode = encodes[accession]
        grouped_accessions.setdefault(encode.experiment_accession, [])
        grouped_accessions[encode.experiment_accession].append(accession)

    # For each output type...
    for out_type in output_types[feat_type]:

        # For each experiment accession...
        for experiment_accession in grouped_accessions:        

            # Skip experiment
            if experiment_accession in done:
                continue

            # For each accession...
            for accession in grouped_accessions[experiment_accession]:

                if encodes[accession].output_type.lower() == out_type:
                    filtered_accessions.add(accession)
                    done.add(experiment_accession)

    return(filtered_accessions)

def _group_ENCODE_accessions(accessions):

    # Initialize
    grouped_accessions = {}

    # For each accession...
    for accession in accessions:

        # Initialize
        experiment_target = encodes[accession].experiment_target
        experiment_type = encodes[accession].experiment_type

        # Group metadata
        grouped_accessions.setdefault((experiment_target, experiment_type), set())
        grouped_accessions[(experiment_target, experiment_type)].add(accession)

    return(grouped_accessions)

def _preprocess_data(accessions, out_dir=".", threads=1):

    # Initialize
    dummy_files = []
    encode_objects = set()

    # Get label
    accession = next(iter(accessions))
    encode = encodes[accession]
    label = encode.experiment_type
    if encode.experiment_target is not None:
        label += ".%s" % encode.experiment_target

    # Skip if BED file exists
    bed_file = os.path.join(out_dir, "%s.bed" % label)
    if not os.path.exists(bed_file):

        # Skip if BED file exists
        dummy_file = os.path.join(out_dir, "dummy.bed")
        if not os.path.exists(dummy_file):

            # For each accession...
            for accession in accessions:
                encode_objects.add(encodes[accession])

            # Get ENCODE BED files
            pool = Pool(processes=threads)
            p = partial(_download_ENCODE_bed_file, out_dir=out_dir)
            for download_file in pool.imap(p, encode_objects):

                # Initialize
                m = re.search("\/(\w+).(bam|bed.gz)$", download_file)
                accession = m.group(1)

                # Concatenate
                cmd = "zless %s | awk '{print $1\"\t\"$2\"\t\"$3\"\t%s\t\"$7\"\t\"$10}' >> %s" % \
                    (download_file, accession, dummy_file)
                subprocess.call(cmd, shell=True)

                # Remove downloaded file
                os.remove(download_file)

            pool.close()

        # Add dummy file
        dummy_files.append(dummy_file)

        # Sort BED
        dummy_file = os.path.join(out_dir, "dummy.sorted.bed")
        if not os.path.exists(dummy_file):

            # UNIX parallel sort is far more efficient than bedtools
            cmd = "LC_ALL=C sort --parallel=%s -T %s -k1,1 -k2,2n %s > %s" % \
                (str(threads), out_dir, dummy_files[0], dummy_file)
            subprocess.call(cmd, shell=True)

        # Add dummy file
        dummy_files.append(dummy_file)

        # Copy file
        shutil.copy(dummy_files[1], bed_file)

        # Remove ALL dummy files
        for dummy_file in dummy_files:
            os.remove(dummy_file)

    return(bed_file)

def _download_ENCODE_bed_file(encode, out_dir="."):

    # Initialize
    download_file = os.path.join(out_dir, encode.accession)

    # Download BED file
    download_file += ".bed.gz"
    if not os.path.exists(download_file):
        urlretrieve(encode.download_url, download_file)

    return(download_file)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
