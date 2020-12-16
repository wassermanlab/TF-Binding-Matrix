#!/usr/bin/env python

import argparse
from functools import partial
import gzip
from multiprocessing import Pool, cpu_count
import numpy as np
import os
import pandas as pd
from pybedtools import BedTool
import re
import sparse
import subprocess as sp

from __init__ import ENCODE, ParseUtils

usage_msg = """
usage: %s --dhs-dir DIR --fasta-file FILE --encode-dir DIR
          --remap-dir DIR --unibind-dir DIR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
builds multiple matrix of bound and open regions across TFs
and cells/tissues from ENCODE, ReMap and UniBind data

  --dhs-file FILE     from get_dhs.sh (e.g. DHS.200bp.bed)
  --encode-dir DIR    output directory from get_encode.py
  --fasta-file FILE   from get_hg38.sh (i.e. hg38.fa)
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
    parser.add_argument("--fasta-file")
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
    if not args.dhs_file or not args.encode_dir or not args.fasta_file or not args.remap_dir or not args.unibind_dir:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error",
            "arguments \"--dhs-file\" \"--encode-dir\" \"--fasta-file\" \"--remap-dir\" \"--unibind-dir\" are required\n"]
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

def _split_data(data_file, threads=1):
    """
    Adapted from [GUD](https://github.com/wassermanlab/GUD)
    """

    # Initialize
    split_files = []
    split_dir = os.path.dirname(os.path.realpath(data_file))

    # Get number of lines
    output = sp.check_output(["wc -l %s" % data_file], shell=True)
    m = re.search("(\d+)", str(output))
    L = float(m.group(1))

    # Split
    prefix = "%s." % data_file.split("/")[-1]
    cmd = "split -d -l %s %s %s" % (int(L/threads)+1, data_file, os.path.join(split_dir, prefix))
    sp.run(cmd, shell=True)

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
        chrom = it.chrom
        start = int(it.start)
        end = int(it.end)
        region = tuple([chrom, start, end])
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
                ilocs.add(ixdf.loc[(x, y, z), "Iloc"])
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
        chrom = it.chrom
        start = int(it.start)
        end = int(it.end)
        region = tuple([chrom, start, end])
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
                ilocs.add(ixdf.loc[(x, y, z), "Iloc"])
            except:
                pass

    return(ilocs)

def main():

    # Parse arguments
    args = parse_args()

    # Build matrices
    build_matrix(args.dhs_file, args.encode_dir, args.fasta_file,
        args.remap_dir, args.unibind_dir, args.out_dir, args.threads)

def build_matrix(dhs_file, encode_dir, fasta_file, remap_dir, unibind_dir,
    out_dir=".", threads=1):
    """
    e.g. ./matrix.py --dhs-file ../../DHS/UCSC/DHS.200bp.bed \
                     --encode-dir ../../ENCODE/hg38/ \
                     --fasta-file ../../Genomes/hg38/hg38.fa \
                     --remap-dir ../../ReMap/ \
                     --unibind-dir ../../UniBind/
    """

    # Globals
    global ixdf
    global regions_idx
    global samples_idx
    global tfs_idx
    global xrefs

    # Get bp extension
    m = re.search("DHS\.(.+)\.bed", dhs_file)
    out_dir = os.path.join(out_dir, m.group(1))
    n = re.search("(\d+)bp", m.group(1))
    bp = int(n.group(1))

    # Create output dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)   

    #######################################################
    # First, figure out the dimensions of the matrix:     #
    # (*) Create indices for x, y and z axis              #
    #######################################################

    # Get encodes
    pkl_file = "metadata.GRCh38.accessibility.tsv.pickle.gz"
    encodes_acc = ParseUtils.load_pickle(os.path.join(encode_dir, pkl_file))
    pkl_file = "metadata.GRCh38.tf.tsv.pickle.gz"
    encodes_tfs = ParseUtils.load_pickle(os.path.join(encode_dir, pkl_file))

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

    # Get UCSC sources
    sources = {}
    ucsc_dir = os.path.dirname(os.path.realpath(dhs_file))
    sources_file = os.path.join(ucsc_dir, "wgEncodeRegDnaseClusteredSources.txt.gz")
    for ix, sample in ParseUtils.parse_tsv_file(sources_file):
        sources.setdefault(int(ix), sample)

    # Get UCSC to ENCODE mappings
    ucsc2encode = {}
    encode2ucsc = {}
    ucsc2encode_file = os.path.join(ucsc_dir, "ucsc2encode.tsv")
    for ucsc, encode in ParseUtils.parse_tsv_file(ucsc2encode_file):
        ucsc2encode.setdefault(ucsc, [])
        ucsc2encode[ucsc].append(encode)
        encode2ucsc.setdefault(encode, [])
        encode2ucsc[encode].append(ucsc)

    # Get UniBind TFs 
    tfs = set()
    file_name = os.path.join(unibind_dir, "tfs.txt")
    for tf in ParseUtils.parse_file(file_name):
        tfs.add(tf)

    # Skip if samples pickle file already exists
    samples_idx_file = os.path.join(out_dir, "samples_idx.pickle.gz")
    xrefs_file = os.path.join(out_dir, ".xrefs.pickle.gz")
    if not os.path.exists(samples_idx_file):
        # Get ENCODE samples with accessibility and TFs
        xrefs = {}
        samples_idx = {}
        samples_acc = set()
        for accession in encodes_acc:
            encode = encodes_acc[accession]
            if encode.biosample_name in encode2ucsc:
                samples_acc.add(encode.biosample_name)
        samples_tfs = set()
        for accession in encodes_tfs:
            encode = encodes_tfs[accession]
            if (encode.experiment_accession in encsr2sample and \
                encode.experiment_target in tfs):
                if encode.biosample_name in encode2ucsc:
                    samples_tfs.add(encode.biosample_name)
                    sample = encsr2sample[encode.experiment_accession]
                    xrefs.setdefault(sample, set())
                    xrefs[sample].add(encode.biosample_name)
        samples = sorted(list(samples_acc.intersection(samples_tfs)))
        for v, k in enumerate(samples):
            samples_idx.setdefault(k, v)
        for sample in xrefs:
            xrefs[sample] = sorted(list(xrefs[sample]))
        # Save pickles
        ParseUtils.save_compressed_pickle(samples_idx_file, samples_idx)
        ParseUtils.save_compressed_pickle(xrefs_file, xrefs)
    else:
        # Load pickles
        samples_idx = ParseUtils.load_compressed_pickle(samples_idx_file)
        xrefs = ParseUtils.load_compressed_pickle(xrefs_file)

    # Skip if regions pickle already exists
    regions_idx_file = os.path.join(out_dir, "regions_idx.pickle.gz")
    if not os.path.exists(regions_idx_file):
        # Get bound and open regions
        regions = []
        a = BedTool(dhs_file)
        for it in a:
            ixs = it.fields[-1].split(",")
            chrom = it.chrom
            start = int(it.start)
            end = int(it.end)
            for ix in ixs[:-1]:
                if ucsc2encode[sources[int(ix)]][0] in samples_idx:
                    regions.append(tuple([chrom, start, end]))
                    break
        regions_idx = {k: v for v, k in enumerate(regions)}
        # Save pickle
        ParseUtils.save_compressed_pickle(regions_idx_file, regions_idx)
    else:
        # Load pickle
        regions_idx = ParseUtils.load_compressed_pickle(regions_idx_file)

    # Skip if TFs pickle already exists
    tfs_idx_file = os.path.join(out_dir, "tfs_idx.pickle.gz")
    if not os.path.exists(tfs_idx_file):
        # Get TFs in ENCODE samples
        tfs = set()
        tfs_idx = {}
        for file_name in os.listdir(unibind_dir):
            if not file_name.endswith(".pwm.bed"):
                continue
            file_name = file_name.split(".")
            sample = file_name[1].upper()
            tf = file_name[2]
            if sample in xrefs:
                tfs.add(tf)
        tfs_idx = {k: v for v, k in enumerate(sorted(list(tfs)))}
        # Save pickle
        ParseUtils.save_compressed_pickle(tfs_idx_file, tfs_idx)
    else:
        # Load pickle
        tfs_idx = ParseUtils.load_compressed_pickle(tfs_idx_file)

    ########################################################
    # Since it is computationally too expensive to build a #
    # 1817918x163x52 matrix (i.e. requires ~123GB of RAM), #
    # instead build a sparse 3D matrix of 89854212 cells.  #
    ########################################################

    # Skip if data frame already initialized
    pandas_file = os.path.join(out_dir, ".data_frame.pickle.gz")
    if not os.path.exists(pandas_file):

        # Build a numpy array of sample-specific regions
        numpy_file = os.path.join(out_dir, ".regions2samples.npz")
        if not os.path.exists(numpy_file):
            m = np.zeros((len(regions_idx), len(samples_idx)))
            a = BedTool(dhs_file)
            for it in a:
                ixs = it.fields[-1].split(",")
                chrom = it.chrom
                start = int(it.start)
                end = int(it.end)
                region = tuple([chrom, start, end])
                if region not in regions_idx:
                    continue
                i = regions_idx[region]
                for ix in ixs[:-1]:
                    if ucsc2encode[sources[int(ix)]][0] in samples_idx:
                        j = samples_idx[ucsc2encode[sources[int(ix)]][0]]
                        m[i][j] = 1
            np.savez_compressed(numpy_file, m)
        else:
            data = np.load(numpy_file)
            for i in data.files:
                m = data[i]

        # Build a numpy array of sample-specific TFs
        numpy_file = os.path.join(out_dir, ".tfs2samples.npz")
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
            np.savez_compressed(numpy_file, n)
        else:
            data = np.load(numpy_file)
            for i in data.files:
                n = data[i]

        # Initialize data frame w/ open regions
        # across TFs and cells/tissues (score = 1)
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
        ParseUtils.save_compressed_pickle(pandas_file, df)
    else:
        df = ParseUtils.load_compressed_pickle(pandas_file)

    # Skip if completed
    completed_file = os.path.join(out_dir, ".remap.completed")
    if not os.path.exists(completed_file):

        # Set index on data frame
        ixdf = df.set_index(["Biosample", "Region", "TF"])
        ixdf["Iloc"] = list(df.index)
        ixdf.loc[(0, 16, 5), "Iloc"]

        # Split data
        data_file = os.path.join(remap_dir, "summits.bed")
        split_files = _split_data(data_file, threads)

        # Upsert data frame w/ ReMap-bound and open
        # regions across TFs and cells/tissues (score = 2)
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
        ixdf = df.set_index(["Biosample", "Region", "TF"])
        ixdf["Iloc"] = list(df.index)
        ixdf.loc[(0, 16, 5), "Iloc"]

        # Split data
        data_file = os.path.join(unibind_dir, "tfbs.bed")
        split_files = _split_data(data_file, threads)

        # Upsert data frame w/ UniBind-bound and open
        # regions (score = 3) or w/ ReMap and UniBind
        # bound and open regions across TFs and
        # cells/tissues (score = 4)
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
        matrix3d = sparse.COO(coords, data, shape=shape)
        sparse.save_npz(matrix_file, matrix3d, compressed=True)

    else:

        # Load sparse matrix
        matrix3d = sparse.load_npz(matrix_file)

    #######################################################
    # Build two matrices for transfer learning.           #
    #######################################################

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap.resolved.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # ReMap bound and open regions (2 and 4) = 1
    #     # Remaining regions (0, 1 and 3) = 0
    #     matrix2d[matrix2d == 0] = 0
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 1
    #     matrix2d[matrix2d == 3] = 0
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists (i.e. for pre-training)
    # matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap.less-sparse.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # ReMap bound and open regions (2 and 4) = 1
    #     # Open regions (1) and UniBind-only bound and open regions (3) = 0
    #     # Remaining regions (0) = 0
    #     matrix2d[matrix2d == 0] = np.nan
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 1
    #     matrix2d[matrix2d == 3] = 0
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap.sparse.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # ReMap bound and open regions (2 and 4) = 1
    #     # Open regions (1) = 0
    #     # Remaining regions (0) and UniBind-only bound and open regions (3) = None
    #     matrix2d[matrix2d == 0] = np.nan
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 1
    #     matrix2d[matrix2d == 3] = np.nan
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.UniBind.resolved.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # UniBind bound and open regions (3 and 4) = 1
    #     # Remaining regions (0, 1 and 2) = 0
    #     matrix2d[matrix2d == 0] = 0
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 0
    #     matrix2d[matrix2d == 3] = 1
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.UniBind.less-sparse.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # UniBind bound and open regions (3 and 4) = 1
    #     # Open regions (1) and ReMap-only bound and open regions (2) = 0
    #     # Remaining regions (0) = 0
    #     matrix2d[matrix2d == 0] = np.nan
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 0
    #     matrix2d[matrix2d == 3] = 1
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.UniBind.sparse.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # UniBind bound and open regions (3 and 4) = 1
    #     # Open regions (1) = 0
    #     # Remaining regions (0) and ReMap-only bound and open regions (2) = None
    #     matrix2d[matrix2d == 0] = np.nan
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = np.nan
    #     matrix2d[matrix2d == 3] = 1
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # # Skip if file already exists
    # matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap+UniBind.resolved.npz")
    # if not os.path.exists(matrix2d_file):

    #     # Collapse into a 1817918x163 matrix
    #     matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
    #     # ReMap AND UniBind bound and open regions (4) = 1
    #     # Remaining regions (0, 1, 2, 3) = 0
    #     matrix2d[matrix2d == 0] = 0
    #     matrix2d[matrix2d == 1] = 0
    #     matrix2d[matrix2d == 2] = 0
    #     matrix2d[matrix2d == 3] = 0
    #     matrix2d[matrix2d == 4] = 1
    #     np.savez_compressed(matrix2d_file, matrix2d)

    # Skip if file already exists
    matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap+UniBind.less-sparse.npz")
    if not os.path.exists(matrix2d_file):

        # Collapse into a 1817918x163 matrix
        matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
        # ReMap AND UniBind bound and open regions (4) = 1
        # Open regions (1), ReMap-/UniBind-only bound and open regions (2 and 3) = 0
        # Remaining regions (0) = None
        matrix2d[matrix2d == 0] = np.nan
        matrix2d[matrix2d == 1] = 0
        matrix2d[matrix2d == 2] = 0
        matrix2d[matrix2d == 3] = 0
        matrix2d[matrix2d == 4] = 1
        np.savez_compressed(matrix2d_file, matrix2d)

    # Skip if file already exists
    matrix2d_file = os.path.join(out_dir, "matrix2d.ReMap+UniBind.sparse.npz")
    if not os.path.exists(matrix2d_file):

        # Collapse into a 1817918x163 matrix
        matrix2d = matrix3d.reduce(np.maximum, axis=0).todense().astype("float")
        # ReMap AND UniBind bound and open regions (4) = 1
        # Open regions (1) = 0
        # Remaining regions (0, 2, 3) = None
        matrix2d[matrix2d == 0] = np.nan
        matrix2d[matrix2d == 1] = 0
        matrix2d[matrix2d == 2] = np.nan
        matrix2d[matrix2d == 3] = np.nan
        matrix2d[matrix2d == 4] = 1
        np.savez_compressed(matrix2d_file, matrix2d)

    # # For each TF...
    # for tf in sorted(tfs_idx, key=lambda x: tfs_idx[x]):

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.ReMap.sparse.npz" % tf)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x52
    #         matrix2d = matrix3d[:,:,tfs_idx[tf]].todense().astype("float")
    #         # ReMap bound and open regions (2 and 4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0) and UniBind-only bound and open regions (3) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = 1
    #         matrix2d[matrix2d == 3] = np.nan
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.UniBind.sparse.npz" % tf)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x52
    #         matrix2d = matrix3d[:,:,tfs_idx[tf]].todense().astype("float")
    #         # UniBind bound and open regions (3 and 4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0) and ReMap-only bound and open regions (2) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = np.nan
    #         matrix2d[matrix2d == 3] = 1
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.ReMap+UniBind.sparse.npz" % tf)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x52
    #         matrix2d = matrix3d[:,:,tfs_idx[tf]].todense().astype("float")
    #         # ReMap AND UniBind bound and open regions (4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0, 2, 3) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = np.nan
    #         matrix2d[matrix2d == 3] = np.nan
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    # # For each sample...
    # for s in sorted(samples_idx, key=lambda x: samples_idx[x]):

    #     # Initialize
    #     biosample_name = s.replace(" ", "_").replace("/", "-")

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.ReMap.sparse.npz" % biosample_name)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x163
    #         matrix2d = matrix3d[samples_idx[s],:,:].todense().astype("float")
    #         # ReMap bound and open regions (2 and 4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0) and UniBind-only bound and open regions (3) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = 1
    #         matrix2d[matrix2d == 3] = np.nan
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.UniBind.sparse.npz" % biosample_name)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x163
    #         matrix2d = matrix3d[samples_idx[s],:,:].todense().astype("float")
    #         # UniBind bound and open regions (3 and 4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0) and ReMap-only bound and open regions (2) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = np.nan
    #         matrix2d[matrix2d == 3] = 1
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    #     # Skip if file already exists
    #     matrix2d_file = os.path.join(out_dir, "matrix2d.%s.ReMap+UniBind.sparse.npz" % biosample_name)
    #     if not os.path.exists(matrix2d_file):
    #         # Select a slice of 3D matrix; size = 1817918x163
    #         matrix2d = matrix3d[samples_idx[s],:,:].todense().astype("float")
    #         # ReMap AND UniBind bound and open regions (4) = 1
    #         # Open regions (1) = 0
    #         # Remaining regions (0, 2, 3) = None
    #         matrix2d[matrix2d == 0] = np.nan
    #         matrix2d[matrix2d == 1] = 0
    #         matrix2d[matrix2d == 2] = np.nan
    #         matrix2d[matrix2d == 3] = np.nan
    #         matrix2d[matrix2d == 4] = 1
    #         np.savez_compressed(matrix2d_file, matrix2d)

    #######################################################
    # Create FASTA files of sequence length 200-1000bp.   #
    #######################################################

    # Initialize
    minlen = 200
    maxlen = 1000
    step = 100

    # For each sequence length...
    for l in range(minlen, maxlen + step, step):

        # Initialize
        chroms = set()
        n = int((l - bp) / 2) # increase n base pairs in each direction

        # Skip if BED file does not exist
        bed_file = os.path.join(out_dir, ".regions.bed")
        if not os.path.exists(bed_file):

            for i in sorted(regions_idx, key=lambda x: regions_idx[x]):
                string = "%s\t%s\t%s\t%s" % (
                    i[0], i[1] - n, i[2] + n, regions_idx[i]
                )
                ParseUtils.write(bed_file, string)
                chroms.add(i[0])

        # Get BED file
        a = BedTool(bed_file)

        for chrom in sorted(chroms):

            # Skip if file already exists
            sequences_file = os.path.join(
                out_dir, "sequences.%sbp.%s.fa.gz" % (l, chrom)
            )
            if not os.path.exists(sequences_file):

                # Get chromosome FASTA file
                c = a.filter(lambda x: x.chrom == chrom)
                c = c.sequence(fi=fasta_file, name=True)

                # Write
                with gzip.open(sequences_file, "wb") as f:
                    f.write(open(c.seqfn).read().encode())

        # Remove
        os.remove(bed_file)

        break

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":
    main()