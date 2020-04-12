from Bio import SeqIO
from io import BytesIO
from ftplib import FTP
import gzip
import numpy as np
import os
import pandas as pd
import pickle
import sys
from zipfile import ZipFile

#-------------#
# Classes     #
#-------------#

class ENCODE:

    def __init__(self, accession, biosample_name, biosample_type, download_url,
        experiment_accession, experiment_type, experiment_target, genome_assembly,
        output_format, output_type, status, treatments, genetic_modifications):

        # Fix hg38
        if genome_assembly == "GRCh38":
            genome_assembly = "hg38"

        self.accession = accession
        self.biosample_name = biosample_name
        self.biosample_type = biosample_type
        self.download_url = download_url
        self.experiment_accession = experiment_accession
        self.experiment_type = experiment_type
        self.experiment_target = experiment_target
        self.genome_assembly = genome_assembly
        self.output_format = output_format
        self.output_type = output_type
        self.status = status
        self.treatments = treatments
        self.genetic_modifications = genetic_modifications

    @property
    def sex(self):
        """
        sex of biosample
        """
        return(self._biosample_sex)

    @sex.setter
    def sex(self, value):
        if value == "female" or value == "male":
            self._biosample_sex = str(value)

    @property
    def summary(self):
        """
        summary of biosample
        """
        return(self._biosample_summary)

    @summary.setter
    def summary(self, value):
        self._biosample_summary = str(value)

    @property
    def X(self):
        """
        number of X chromosomes
        """
        if self.sex == "female":
            return(2)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

    @property
    def Y(self):
        """
        number of Y chromosomes
        """
        if self.sex == "female":
            return(0)
        elif self.sex == "male":
            return(1)
        else:
            return(None)

    @property
    def treatment(self):
        """
        is biosamble treated?
        """
        return(self.treatments is not None)

    @property
    def genetic_modification(self):
        """
        has biosample been genetically modified?
        """
        return(self.genetic_modifications is not None)

class ParseUtililities:

    #++++++++++++++#
    # Input/Output #
    #++++++++++++++#

    def _get_file_handle(self, file_name, mode="r"):

        # Initialize
        raiseValueError = False
        
        # Open file handle
        if file_name.endswith(".gz"):
            try:
                handle = gzip.open(file_name, mode)
            except:
                raiseValueError = True

        elif file_name.endswith(".zip"):
            try:
                zf = ZipFile(file_name, mode)
                for f in zf.infolist():

                    handle = zf.open(f, mode)
                    break
            except:
                raiseValueError = True

        else:
            try:
                handle = open(file_name, mode)
            except:
                raiseValueError = True
        
        if raiseValueError:
            raise ValueError("Could not open file handle: %s" % file_name)

        return(handle)

    def parse_file(self, file_name):
        """
        Parses a file and yields lines one by one.

        @input:
        file_name {str}

        @yield: {str}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each line...
        for line in handle:

            yield(line.strip("\n"))

        handle.close()

    def parse_csv_file(self, file_name, delimiter=","):
        """
        Parses a CSV file and yields lines one by one as a list.

        @input:
        file_name {str}
        delimiter {str} e.g. "\t"; default = ","

        @yield: {list}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # Read in chunks
        for chunk in pd.read_csv(handle, header=None, encoding="utf8", sep=delimiter, chunksize=1024, comment="#"):
            for index, row in chunk.iterrows():
                yield(row.tolist())

        handle.close()

    def parse_tsv_file(self, file_name):
        """
        Parses a TSV file and yields lines one by one as a list.

        @input:
        file_name {str}

        @yield: {list}
        """

        # For each line...
        for line in self.parse_csv_file(file_name, delimiter="\t"):

            yield(line)

    def parse_fasta_file(self, file_name):
        """
        Parses a FASTA file and yields sequences one by one  as a list of
        length 2 (i.e. [{header}, {sequence}]).

        @input:
        file_name {str}

        @yield: {list}
        """

        # Get file handle
        handle = self._get_file_handle(file_name)

        # For each SeqRecord...
        for seq_record in SeqIO.parse(handle, "fasta"):
            # Initialize
            header = seq_record.id
            sequence = str(seq_record.seq).upper()

            yield(header, sequence)

        handle.close()

    def write(self, file_name=None, content=None):
        """
        Writes content to a file or, if no file is provided, to STDOUT.
        Content will be appended at the end of the file.
        """

        if file_name:

            # Get file handle
            handle = self._get_file_handle(file_name, mode="a")

            # Write
            handle.write("%s\n" % content)
            handle.close()

        else:

            sys.stdout.write("%s\n" % content)

    #++++++++++++++++#
    # Pickle recipes #
    #++++++++++++++++#

    # Adapted from "Load Faster in Python With Compressed Pickles"
    # https://medium.com/better-programming/load-fast-load-big-with-compressed-pickles-5f311584507e

    def save_pickle(self, file_name, data_obj):
        """
        Saves a pickle.
        """

        # Get file handle
        handle = self._get_file_handle(file_name, mode="wb")

        # Save
        pickle.dump(data_obj, handle)
        handle.close()

    def load_pickle(self, file_name):
        """
        Loads and returns a pickle.
        """

        # Get file handle
        handle = self._get_file_handle(file_name, mode="rb")

        # Load
        pkl = pickle.load(handle)
        handle.close()

        return(pkl)

    def save_compressed_pickle(self, file_name, data_obj):
        """
        Saves a compressed pickle.
        """

        # Get file handle
        handle = self._get_file_handle(file_name, mode="wb")

        # Save
        pickle.dump(data_obj, handle, protocol=4)
        handle.close()

    def load_compressed_pickle(self, file_name):
        """
        Loads and returns a compressed pickle.
        """

        # Get file handle
        handle = self._get_file_handle(file_name, mode="rb")

        # Load
        pkl = pickle.load(handle)
        handle.close()

        return(pkl)

    #+++++++++++++++#
    # NumPy recipes #
    #+++++++++++++++#

    def load_compressed_npz(self, file_name):
        """
        Loads and returns a compressed numpy array.
        """

        data = np.load(file_name)
        for i in data.files:
            return(data[i])

ParseUtils = ParseUtililities()