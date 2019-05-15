#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import os
import sys
import shutil
import glob
import subprocess
import random
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord

__author__ = "Rick Conzemius"
__copyright__ = "Copyright 2019, Rick Conzemius"
__license__ = "MIT"
__maintainer__ = "Rick Conzemius"
__version__ = "1.0"
__email__ = "rick.conzemius.fl@ait.ac.at"
__status__ = "Production"

class IndexBuilder:
    def __init__(self, fasta_file):
        self.base_folder   = os.getcwd() + "/"
        self.input_folder  = self.base_folder + "input/"
        self.bowtie_folder = self.base_folder + "bowtie/"
        self.genome_folder = self.base_folder + "genomes/"

        self.fasta_file = fasta_file
        self.fasta_file_path = self.input_folder + self.fasta_file
        self.fasta_file_without_ext = self.fasta_file.replace(".fasta", "")
        self.output_fasta_file = self.genome_folder + self.fasta_file

        self.num_threads = 32
        self.bowtie_build_cmd = "/home/development/bigData2/Rick/Deployment/PRIMEval/Apps/bowtie-1.2.2/bowtie-build"
        self.faidx_cmd = "/home/development/bigData2/Rick/samtools/samtools faidx "
        self.separator = "__contigname__"

    def prepare_folders(self):
        if not os.path.exists(self.bowtie_folder):
            os.makedirs(self.bowtie_folder)
        if not os.path.exists(self.genome_folder):
            os.makedirs(self.genome_folder)

    def import_contigs(self):
        seq_records = []
        for entry in SeqIO.parse(self.fasta_file_path, "fasta"):
             new_id = self.fasta_file_without_ext + self.separator + entry.id
             seq_records.append(SeqRecord(entry.seq, id = new_id, description = ""))
        output_handle = open(self.output_fasta_file, "w")
        SeqIO.write(seq_records, output_handle, "fasta")
        output_handle.close()

    def faidx_index(self):
        command = self.faidx_cmd + self.output_fasta_file
        subprocess.call(command, shell = True)

    def create_bowtie_index(self):
        command = self.bowtie_build_cmd + " --threads " + str(self.num_threads) + " -r -f " + self.output_fasta_file + " " + self.bowtie_folder + self.fasta_file_without_ext
        subprocess.call(command, shell = True)

if __name__ == '__main__':
    for file in os.listdir("input/"):
        if file.endswith(".fasta"):
            print("Importing " + str(file))
            index = IndexBuilder(file)
            index.prepare_folders()
            index.import_contigs()
            index.faidx_index()
            index.create_bowtie_index()
            print("Finished importing " + str(file))

