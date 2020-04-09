#!/usr/bin/env python

import argparse
from functools import partial
from multiprocessing import Pool, cpu_count
import numpy as np
import os
import pandas as pd
import pickle
from pybedtools import BedTool, Interval
import re
import sparse
import subprocess

from __init__ import ENCODE, ParseUtils

usage_msg = """
usage: %s --dhs-dir DIR --encode-dir DIR --remap-dir DIR
          --unibind-dir DIR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
builds a 3D matrix of bound and open regions across TFs and
biosamples from ENCODE, ReMap and UniBind data

  --dhs-file FILE     from get_dhs.sh (e.g. DHS.200bp.bed)
  --encode-dir DIR    output directory from get_encode.py
  --remap-dir DIR     output directory from get_remap.sh
  --unibind-dir DIR   output directory from get_unibind.sh

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
    parser.add_argument("--dhs-file")
    parser.add_argument("--encode-dir")
    parser.add_argument("--remap-dir")
    parser.add_argument("--unibind-dir")

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
    if not args.dhs_file or not args.encode_dir or not args.remap_dir or not args.unibind_dir:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--dhs-file\" \"--encode-dir\" \"--remap-dir\" \"--unibind-dir\" are required\n"]
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

def _create_indices(dhs_file, encode_dir, remap_dir, unibind_dir,
    out_dir=".", threads=1):

    # Get encodes
    pickle_file = os.path.join(encode_dir, "metadata.GRCh38.accessibility.tsv.pickle")
    handle = ParseUtils._get_file_handle(pickle_file, "rb")
    encodes_acc = pickle.load(handle)
    pickle_file = os.path.join(encode_dir, "metadata.GRCh38.tf.tsv.pickle")
    handle = ParseUtils._get_file_handle(pickle_file, "rb")
    encodes_tfs = pickle.load(handle)

    # Get ReMap/UniBind ENCODE assay accessions (ENCSRs) and sample names
    encsr2sample = {}
    file_name = os.path.join(remap_dir, "encsr2sample.txt")
    for line in ParseUtils.parse_file(file_name):
        encsr, sample = line.split(".")
        encsr2sample.setdefault(encsr, sample.upper())
    file_name = os.path.join(unibind_dir, "encsr2sample.txt")
    for line in ParseUtils.parse_file(file_name):
        encsr, sample = line.split(".")
        encsr2sample.setdefault(encsr, sample.upper())

    # Get UniBind TFs 
    tfs = set()
    file_name = os.path.join(unibind_dir, "tfs.txt")
    for tf in ParseUtils.parse_file(file_name):
        tfs.add(tf)

    # Skip if samples index file already exists
    samples_idx_file = os.path.join(out_dir, "samples_idx.pickle")
    xrefs_file = os.path.join(out_dir, ".xrefs.pickle")
    if not os.path.exists(samples_idx_file):
        # Get ENCODE samples with accessibility and TFs
        xrefs = {}
        samples_acc = set()
        for accession in encodes_acc:
            encode = encodes_acc[accession]
            samples_acc.add(encode.biosample_name)
        samples_tfs = set()
        for accession in encodes_tfs:
            encode = encodes_tfs[accession]
            if (encode.experiment_accession in encsr2sample and \
                encode.experiment_target in tfs):
                samples_tfs.add(encode.biosample_name)
                if encode.biosample_name in samples_acc:
                    sample = encsr2sample[encode.experiment_accession]
                    xrefs.setdefault(sample, set())
                    xrefs[sample].add(encode.biosample_name)
        samples_list = sorted(list(samples_acc.intersection(samples_tfs)))
        for v, k in enumerate(samples_list):
            samples_idx.setdefault(v, k)
            samples_idx.setdefault(k, v)
        for sample in xrefs:
            xrefs[sample] = sorted(list(xrefs[sample]))
        # Write
        ParseUtils.write(samples_idx_file, json.dumps(samples_idx, **opts))
        ParseUtils.write(xrefs_file, json.dumps(xrefs, **opts))
    else:
        with open(samples_idx_file) as f:
            samples_idx = json.load(f)
        with open(xrefs_file) as f:
            xrefs = json.load(f)

    # Skip if TFs file already exists
    tfs_idx_file = os.path.join(out_dir, "tfs_idx.json")
    if not os.path.exists(tfs_idx_file):
        # Get TFs in ENCODE samples
        tfs = set()
        for file_name in os.listdir(unibind_dir):
            if not file_name.endswith(".pwm.bed"):
                continue
            file_name = file_name.split(".")
            sample = file_name[1].upper()
            tf = file_name[2]
            if sample in xrefs:
                tfs.add(tf)
        tfs_idx = {k: v for v, k in enumerate(sorted(list(tfs)))}
        # Write
        ParseUtils.write(tfs_idx_file, json.dumps(tfs_idx, **opts))
    else:
        with open(tfs_idx_file) as f:
            tfs_idx = json.load(f)

    # Skip if regions file already exists
    regions_idx_file = os.path.join(out_dir, "regions_idx.json")
    if not os.path.exists(regions_idx_file):
        # Get BED files
        a = BedTool(dhs_file)
        b = BedTool(os.path.join(encode_dir, "DNase-seq.150bp.bed"))
        # Get intersections
        intervals = []
        its = a.intersect(b, wa=True, wb=True, sorted=True, F=0.5, stream=True)
        for it in its:
            encode = encodes_acc[it.fields[-1]]
            if encode.biosample_name in samples_idx:
                chrom = it.chrom
                start = int(it.start)
                end = int(it.end)
                interval = Interval(chrom, start, end)
                if len(intervals) == 0:
                    intervals.append(interval)
                elif interval != intervals[-1]:
                    intervals.append(interval)
        # Get regions
        regions = BedTool("".join(map(str, intervals)), from_string=True)
        regions_idx = {k: v for v, k in enumerate([
                "%s:%s-%s" % (r.chrom, r.start, r.end) for r in regions
            ])
        }
        # Write
        ParseUtils.write(regions_idx_file, json.dumps(regions_idx, **opts))
    else:
        with open(regions_idx_file) as f:
            regions_idx = json.load(f)


def _split_data(data_file, threads=1):
    """
    Adapted from [GUD](https://github.com/wassermanlab/GUD)
    """

    # Initialize
    split_files = []
    split_dir = os.path.dirname(os.path.realpath(data_file))

    # Get number of lines
    output = subprocess.check_output(["wc -l %s" % data_file], shell=True)
    m = re.search("(\d+)", str(output))
    L = float(m.group(1))

    # Split
    prefix = "%s." % data_file.split("/")[-1]
    cmd = "split -d -l %s %s %s" % (int(L/threads)+1, data_file, os.path.join(split_dir, prefix))
    subprocess.run(cmd, shell=True)

    # For each split file...
    for split_file in os.listdir(split_dir):

        # Append split file
        if split_file.startswith(prefix):
            split_files.append(os.path.join(split_dir, split_file))

    return(split_files)

def _upsert_ReMap(B, A):
    """
    Return data frame index locations to be upserted.
    """

    # Initialize
    ilocs = set()

    # Get BED files
    a = BedTool(A)
    b = BedTool(B)

    # Get intersections
    its = a.intersect(b, wa=True, wb=True, sorted=True, stream=True)

    # For each intersection...
    for it in its:
        fields = it.fields[-1].split(".")
        sample = fields[2].upper()
        tf = fields[1]
        if tf not in tfs_idx:
            continue
        if sample not in xrefs:
            continue
        region = "%s:%s-%s" % (it.chrom, it.start, it.end)
        if region not in regions_idx:
            continue
        for s in xrefs[sample]:
            if s not in samples_idx:
                continue
            # Get data frame idx
            x = samples_idx[s]
            y = regions_idx[region]
            z = tfs_idx[tf]
            try:
                ilocs.add(idf.loc[(x, y, z), "Iloc"])
            except:
                pass

    return(ilocs)

def _upsert_UniBind(B, A):
    """
    Return data frame index locations to be upserted.
    """

    # Initialize
    ilocs = set()

    # Get BED files
    a = BedTool(A)
    b = BedTool(B)

    # Get intersections
    its = a.intersect(b, wa=True, wb=True, sorted=True, stream=True)

    # For each intersection...
    for it in its:
        fields = it.fields[-1].split(".")
        sample = fields[1].upper()
        tf = fields[2]
        if tf not in tfs_idx:
            continue
        if sample not in xrefs:
            continue
        region = "%s:%s-%s" % (it.chrom, it.start, it.end)
        if region not in regions_idx:
            continue
        for s in xrefs[sample]:
            if s not in samples_idx:
                continue
            # Get data frame idx
            x = samples_idx[s]
            y = regions_idx[region]
            z = tfs_idx[tf]
            try:
                ilocs.add(idf.loc[(x, y, z), "Iloc"])
            except:
                pass

    return(ilocs)

def main():

    # Parse arguments
    args = parse_args()

    # Build matrix
    build_matrix(args.dhs_file, args.encode_dir, args.remap_dir,
        args.unibind_dir, args.out_dir, args.threads)

def build_matrix(dhs_file, encode_dir, remap_dir, unibind_dir,
    out_dir=".", threads=1):
    """
    e.g. ./build_3d_matrix.py --dhs-file ../DHS/DHS.200bp.bed --encode-dir ../ENCODE/hg38/
                              --remap-dir ../ReMap/ --unibind-dir ../UniBind/
    """

    # Initialize
    global idf
    global regions_idx
    global samples_idx
    global tfs_idx
    global xrefs
    opts = {"indent":4, "separators":(",", ": ")}

    # Get bp extension
    m = re.search("DHS\.(.+)\.bed", dhs_file)
    out_dir = os.path.join(out_dir, m.group(1))

    # Create output dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)   

    # Create matrix indices
    samples_idx, regions_idx, tfs_idx = _create_indices(dhs_file,
        encode_dir, remap_dir, unibind_dir, out_dir, threads)

    #######################################################
    # Since it is not possible to build a 2431675x174x66  #
    # 3D matrix (i.e. requires 223.4028456 GB of memory), #
    # instead build a sparse 3D matrix.                   #
    #######################################################

    # Skip if data frame already initialized
    pandas_file = os.path.join(out_dir, "data_frame.pickle")
    if not os.path.exists(pandas_file):

        # Build a numpy array of sample-specific regions
        numpy_file = os.path.join(out_dir, ".regions2samples.npy")
        if not os.path.exists(numpy_file):
            m = np.zeros((len(regions_idx), len(samples_idx)))
            a = BedTool(dhs_file)
            b = BedTool(os.path.join(encode_dir, "DNase-seq.150bp.bed"))
            its = a.intersect(b, wa=True, wb=True, sorted=True, F=0.5, stream=True)
            for it in its:
                encode = encodes_acc[it.fields[-1]]
                if encode.biosample_name not in samples_idx:
                    continue
                region = "%s:%s-%s" % (it.chrom, it.start, it.end)
                if region not in regions_idx:
                    continue
                i = regions_idx[region]
                j = samples_idx[encode.biosample_name]
                m[i][j] = 1
            np.save(numpy_file, m)
        else:
            m = np.load(numpy_file)

        # Build a numpy array of sample-specific TFs
        numpy_file = os.path.join(out_dir, ".tfs2samples.npy")
        if not os.path.exists(numpy_file):
            n = np.zeros((len(tfs_idx), len(samples_idx)))
            for file_name in os.listdir(unibind_dir):
                if not file_name.endswith(".pwm.bed"):
                    continue
                file_name = file_name.split(".")
                sample = file_name[1].upper()
                tf = file_name[2]
                if tf not in tfs_idx:
                    continue
                if sample not in xrefs:
                    continue
                for s in xrefs[sample]:
                    if s not in samples_idx:
                        continue
                    i = tfs_idx[tf]
                    j = samples_idx[s]
                    n[i][j] = 1
            np.save(numpy_file, n)
        else:
            n = np.load(numpy_file)

        # Initialize data frame w/ open regions
        # across TFs and biosamples (score = 1)
        data = []
        for k1, v1 in samples_idx.items():
            for k2, v2 in regions_idx.items():
                # i.e. closed
                if m[v2][v1] == 0:
                    continue
                # i.e. open
                for k3, v3 in tfs_idx.items():
                    # i.e. not chip'ed
                    if n[v3][v1] == 0:
                        continue
                    # i.e. chip'ed
                    data.append([v1, v2, v3, 1])
        df = pd.DataFrame(data, columns = ["Biosample", "Region", "TF", "Score"])
        df.to_pickle(pandas_file)

    else:
        df = pd.read_pickle(pandas_file)

    # Skip if completed
    completed_file = os.path.join(out_dir, ".remap.completed")
    if not os.path.exists(completed_file):

        # Set index on data frame
        idf = df.set_index(["Biosample", "Region", "TF"])
        idf["Iloc"] = list(df.index)
        idf.loc[(0, 13, 5), "Iloc"]

        # Split data
        data_file = os.path.join(remap_dir, "summits.bed")
        split_files = _split_data(data_file, threads)

        # Upsert data frame w/ ReMap-bound and open
        # regions across TFs and biosamples (score = 2)
        pool = Pool(threads)
        parallel_upsert = partial(_upsert_ReMap, A=dhs_file)
        for ilocs in pool.imap(parallel_upsert, split_files):
            df.Score.iloc[list(ilocs)] = 2
        df.to_pickle(pandas_file)
        open(completed_file, "a").close()

    # Skip if completed
    completed_file = os.path.join(out_dir, ".unibind.completed")
    if not os.path.exists(completed_file):

        # Set index on data frame
        idf = df.set_index(["Biosample", "Region", "TF"])
        idf["Iloc"] = list(df.index)
        idf.loc[(0, 13, 5), "Iloc"]

        # Split data
        data_file = os.path.join(unibind_dir, "tfbs.bed")
        split_files = _split_data(data_file, threads)

        # Upsert data frame w/ UniBind-bound and open
        # regions (score = 3) or w/ ReMap and UniBind
        # co-bound and open regions across TFs and
        # biosamples (score = 4)
        pool = Pool(threads)
        parallel_upsert = partial(_upsert_UniBind, A=dhs_file)
        for ilocs in pool.imap(parallel_upsert, split_files):
            cobound = []
            for iloc in ilocs:
                # i.e. co-bound
                if df.at[iloc, "Score"] == 2:
                    cobound.append(iloc)
            df.Score.iloc[list(ilocs)] = 3
            df.Score.iloc[cobound] = 4
        df.to_pickle(pandas_file)
        open(completed_file, "a").close()

    # Skip if file already exists
    matrix_file = os.path.join(out_dir, "matrix3d.npz")
    if not os.path.exists(matrix_file):

        # Create matrix
        coords = []
        coords.append(df["Biosample"].tolist())
        coords.append(df["Region"].tolist())
        coords.append(df["TF"].tolist())
        data = df["Score"].tolist()
        shape = (len(samples_idx), len(regions_idx), len(tfs_idx))
        s = sparse.COO(coords, data, shape=shape)
        sparse.save_npz(matrix_file, s, compressed=True)

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()
