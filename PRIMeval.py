#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ahocorasick
import math
import os
import re
import sys
import shutil
import glob
import xlsxwriter
import subprocess
from functools import partial
from itertools import product, combinations
from subprocess import DEVNULL
from multiprocessing import Pool
from threading import Timer
import random
import pandas as pd
import tqdm
import primer3
from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

__author__ = "Rick Conzemius"
__copyright__ = "Copyright 2019, Rick Conzemius"
__license__ = "MIT"
__maintainer__ = "Rick Conzemius"
__version__ = "1.0"
__email__ = "rick.conzemius.fl@ait.ac.at"
__status__ = "Production"

class PRIMeval:
    def __init__(self, run_id, max_primer_mismatches, max_probe_mismatches, max_product_size, cross_check, probes_only, method, dimer_check,
                 primerMonovalentCations, primerDivalentCations, primerDNTPs, primerConcentration, primerAnnealingTemp,
                 probeMonovalentCations, probeDivalentCations, probeDNTPs, probeConcentration, probeAnnealingTemp, prebuilt = ""):
        # parameters
        self.run_id = run_id
        self.max_primer_mismatches = int(max_primer_mismatches)
        self.max_probe_mismatches = int(max_probe_mismatches)
        self.max_product_size = int(max_product_size)
        self.max_mismatches = max(max_primer_mismatches, max_probe_mismatches)
        if self.max_mismatches == 0:
            self.l, self.e, self.qcov, self.perciden = 5, 10, 100, 100
        elif self.max_mismatches == 1:
            self.l, self.e, self.qcov, self.perciden = 5, 40, 90, 90
        elif self.max_mismatches == 2:
            self.l, self.e, self.qcov, self.perciden = 5, 70, 85, 85
        else:
            self.l, self.e, self.qcov, self.perciden = 5, 100, 80, 80
        self.prebuilt = str(prebuilt)
        self.bowtie_dbs = ["list_of_prebuilt_dbs_in_prebuilt_folder"]
        self.bowtie_runs = []
        self.blast_db_name = "user_db"
        self.bowtie_index_name = "bindex"
        self.num_threads = 48
        self.method = method
        if dimer_check == "True":
            self.dimer_check = True
        else:
            self.dimer_check = False
        if cross_check == "True":
            self.same_package = False
        else:
            self.same_package = True
        if probes_only == "True":
            self.probes_only = True
        else:
            self.probes_only = False

        # Cross dimer check
        self.primer_monovalent_cations = str(primerMonovalentCations)
        self.primer_divalent_cations = str(primerDivalentCations)
        self.primer_dntps = str(primerDNTPs)
        self.primer_annealing_oligo = str(primerConcentration)
        self.primer_annealing_temp = str(primerAnnealingTemp)
        self.probe_monovalent_cations = str(probeMonovalentCations)
        self.probe_divalent_cations = str(probeDivalentCations)
        self.probe_dntps = str(probeDNTPs)
        self.probe_annealing_oligo = str(probeConcentration)
        self.probe_annealing_temp = str(probeAnnealingTemp)
        self.cross_dimer_dfs = []
        self.cross_dimer_dfs_dg = []
        self.hairpin_dfs = []

        # Aho-Corasick Automaton
        self.aho = ahocorasick.Automaton()

        # folders
        self.base_folder = os.getcwd() + "/"
        self.run_folder = self.base_folder + "runs/" + str(self.run_id) + "/"
        self.input_folder = self.run_folder + "input/"
        self.output_folder = self.run_folder + "output/"
        self.tmp_folder = self.run_folder + "tmp/"
        self.input_contigs = self.run_folder + "input/contigs/"
        self.primer_input_folder = self.run_folder + "input/primers/"
        self.probes_input_folder = self.run_folder + "input/probes/"
        self.blast_db_folder = self.run_folder + "tmp/blastdb/"
        self.prebuilt_genomes = self.base_folder + "prebuilt/genomes/"
        self.prebuilt_bowtie = self.base_folder + "prebuilt/bowtie/"

        # files
        self.output_contigs = self.run_folder + "tmp/merged_contigs.fasta"
        self.blast_output_tmp_file = self.run_folder + "tmp/blast_tmp_results.txt"
        self.blast_output_file = self.run_folder + "tmp/blast_results.txt"
        self.bowtie_output_tmp_file = self.run_folder + "tmp/bowtie_tmp_results.txt"
        self.bowtie_output_file = self.run_folder + "tmp/bowtie_results.txt"
        self.bowtie_index_folder = self.run_folder + "tmp/bowtie_index_folder/"
        self.oligo_file = self.run_folder + "output/oligos.fasta"
        self.results_all = self.run_folder + "output/results.csv"
        self.results_wob = self.run_folder + "output/results_wobbled.csv"
        self.results_dimers = self.run_folder + "output/results_dimers.xlsx"

        # settings
        self.blastdb_cmd = "/path/to/makeblastdb"
        self.bowtie_build_cmd = "/path/to/bowtie-build"
        self.blast_cmd = "/path/to/blastn"
        self.bowtie_cmd = "/path/to/bowtie"
        self.faidx_cmd = "/path/to/samtools faidx "
        self.pd_col_hits = ["Sequence", "Type", "Name", "Package", "StartPos", "EndPos", "MismatchesTotal",
                            "Strand", "HitSequence", "Tm", "dG"]
        self.pd_col_results = ["Sequence", "Contig", "Primer1", "Primer2", "Probe", "Primer1Package",
                               "Primer2Package", "ProbePackage", "StartPos1", "EndPos1", "StartPos2", "EndPos2",
                               "StartPos3", "EndPos3", "Primer1Tm", "Primer2Tm", "ProbeTm", "Primer1dG", "Primer2dG", "ProbedG", "ProductSize", "ProductTm", "NoMismatchesLeft", "NoMismatchesRight",
                               "NoMismatchesProbe", "MismatchesLeft", "MismatchesRight", "MismatchesProbe",
                               "Comment", "Product"]
        self.blast_txt_params = "\"6 qseqid sseqid nident qlen length mismatch qstart qend sstart sseq sstrand " \
                                "send\""
        self.blast_txt_fields = ["qseqid", "sseqid", "nident", "qlen", "length", "mismatch", "qstart", "qend",
                                 "sstart", "sseq", "sstrand", "send"]

    # return list of all possible sequences given an ambiguous DNA input
    def _extend_ambiguous_dna(self, seq):
        d = Seq.IUPAC.IUPACData.ambiguous_dna_values
        return list(map("".join, product(*map(d.get, seq))))

    def _get_sequence(self, contig_file, wanted_contig, start, end, strand=1):
        try:
            command = self.faidx_cmd + contig_file + " '" + wanted_contig + ":" + str(start) + "-" + str(end) + "'"
            call = subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL).decode().split("\n", 1)[1]
        except:
            try:
                contig_file = self.prebuilt_genomes + wanted_contig.split("__contigname__", 1)[0] + ".fasta"
                command = self.faidx_cmd + contig_file + " '" + wanted_contig + ":" + str(start) + "-" + str(end) + "'"
                call = subprocess.check_output(command, shell=True, stderr=subprocess.DEVNULL).decode().split("\n", 1)[1]
            except:
                sys.exit("Failed retrieving: " + command)
        call = re.sub("\n|\r", "", call)
        sequence = Seq.Seq(call)
        if strand == 1:
            return sequence.upper()
        else:
            return sequence.reverse_complement().upper()

    # Get a visual representation of mismatches between two sequences
    def _mismatch_visualization(self, seq_a, seq_b):
        seq_a, seq_b = seq_a.upper(), seq_b.upper()
        mismatches = ""
        if (len(seq_a) - len(seq_b) != 0):
            return "Error"
        for pos in range(0, len(seq_a)):
            if seq_a[pos] != seq_b[pos]:
                mismatches += "(" + seq_a[pos] + "/" + seq_b[pos] + ")"
            else:
                mismatches += "="
        return mismatches

    def _prepare_folders(self):
        if os.path.exists(self.output_folder):
            shutil.rmtree(self.output_folder)
        if os.path.exists(self.tmp_folder):
            shutil.rmtree(self.tmp_folder)
        # Create output and tmp folders
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        if not os.path.exists(self.tmp_folder):
            os.makedirs(self.tmp_folder)
        if not os.path.exists(self.bowtie_index_folder):
            os.makedirs(self.bowtie_index_folder)

    def _clean_up_folders(self):
        # if os.path.exists(self.input_folder):
        #   shutil.rmtree(self.input_folder)
        if os.path.exists(self.tmp_folder):
            shutil.rmtree(self.tmp_folder)

    # Rename primer and probes, create new sequences without IUPAC codes and save in file
    # Used for dimer check: self.packages, packages, package, oligo_name
    def _import_oligos(self, folder, oligotype):
        packages = {}
        primer_records = []
        allowed_chars = "[^0-9a-zA-Z()'_\+-]+"
        for file in os.listdir(folder):
            if file.endswith(".fasta"):
                package = file.rsplit(".fasta", 1)[0]
                packages[package] = {}
                sequences = SeqIO.parse(open(folder + file), "fasta")
                for fasta in sequences:
                    m = re.search("[M,R,W,S,Y,K,V,H,D,B,N]", str(fasta.seq))
                    if m:
                        sequence_mutations = self._extend_ambiguous_dna(str(fasta.seq))
                        mutation_count = 0
                        for mutation in sequence_mutations:
                            mutation_count += 1
                            oligo_name = re.sub(allowed_chars, "_", fasta.description) + "_mut" + str(mutation_count)
                            packages[package][oligo_name] = str(mutation)
                            if oligotype == "probe":
                                rec = SeqRecord(Seq.Seq(mutation, IUPAC),
                                                id=package + "^" + re.sub(allowed_chars, "_",
                                                                                        fasta.description) + "_mut" + str(
                                                    mutation_count) + "_probe", description="")
                            else:
                                rec = SeqRecord(Seq.Seq(mutation, IUPAC),
                                                id=package + "^" + re.sub(allowed_chars, "_",
                                                                                        fasta.description) + "_mut" + str(
                                                    mutation_count), description="")
                            primer_records.append(rec)
                    else:
                        oligo_name = re.sub(allowed_chars, "_", fasta.description)
                        packages[package][oligo_name] = str(fasta.seq)
                        if oligotype == "probe":
                            rec = SeqRecord(fasta.seq, id=package + "^" + re.sub(allowed_chars, "_",
                                                                                               fasta.description) + "_probe",
                                            description="")
                        else:
                            rec = SeqRecord(fasta.seq,
                                            id=package + "^" + re.sub(allowed_chars, "_", fasta.description),
                                            description="")
                        primer_records.append(rec)
        output_handle = open(self.oligo_file, "a")
        SeqIO.write(primer_records, output_handle, "fasta")
        output_handle.close()
        if oligotype == "primer":
            self.primer_packages = packages
        else:
            self.probe_packages = packages

    # Rename and merge contigs
    def _import_contigs(self):
        seq_records = []
        for file in os.listdir(self.input_contigs):
            # CHANGE: other file endings should also be possible (see with Django upload permitted filenames)
            if file.endswith(".fasta"):
                base_contig_name = file.replace(".fasta", "")
                for entry in SeqIO.parse(self.input_contigs + file, "fasta"):
                    my_new_id = base_contig_name + "__contigname__" + entry.id
                    seq_records.append(SeqRecord(entry.seq, id=my_new_id, description=""))
        output_handle = open(self.output_contigs, "w")
        SeqIO.write(seq_records, output_handle, "fasta")
        output_handle.close()

        command = self.faidx_cmd + self.output_contigs
        subprocess.call(command, shell=True)

    def _import_sequences(self):
        if self.probes_only == False:
            self._import_oligos(self.primer_input_folder, "primer")
        self._import_oligos(self.probes_input_folder, "probe")
        self._import_contigs()

    def _create_blast_db(self):
        command = self.blastdb_cmd + " -in " + self.output_contigs + " -dbtype nucl -out " + self.blast_db_folder + self.blast_db_name
        subprocess.call(command, shell=True)

    def _create_bowtie_index(self):
        command = self.bowtie_build_cmd + " --threads " + str(
            self.num_threads) + " -f " + self.output_contigs + " " + self.bowtie_index_folder + self.bowtie_index_name
        subprocess.call(command, shell=True)

    def _blast_call(self):
        command = self.blast_cmd + " -db " + self.blast_db_folder + self.blast_db_name + " -query " + self.oligo_file + " -out " + \
                  self.blast_output_tmp_file + " -outfmt " + self.blast_txt_params + " -num_threads " + str(
            self.num_threads) + " -evalue 200000 " \
                                "-qcov_hsp_perc " + str(self.qcov) + " -perc_identity " + str(self.perciden) + " -max_target_seqs 2000000 -word_size 4 -ungapped"
        subprocess.call(command, shell=True)
        with open(self.blast_output_file, "a") as out_file:
            with open(self.blast_output_tmp_file) as in_file:
                out_file.write(in_file.read())

    def _bowtie_call(self, index_folder = "", index_name = ""):
        mismatches = self.max_primer_mismatches if self.max_primer_mismatches >= self.max_probe_mismatches else self.max_probe_mismatches
        if index_folder == "" and index_name == "":
            if os.path.getsize(self.output_contigs) == 0:
                return
            index_folder = self.bowtie_index_folder
            index_name = self.bowtie_index_name
        command = self.bowtie_cmd + " -f -a -p " + str(self.num_threads) + " -n " + str(
            mismatches) + " -l " + str(self.l) + " -e " + str(self.e) + " " + index_folder + index_name + " " + self.oligo_file + " " + self.bowtie_output_tmp_file
        subprocess.call(command, shell=True)
        with open(self.bowtie_output_file, "a") as out_file:
            with open(self.bowtie_output_tmp_file) as in_file:
                out_file.write(in_file.read())

    def _specificity_calls(self):
        for db in self.bowtie_runs:
            self._bowtie_call(self.prebuilt_bowtie, db)

    def _multiprocess_convert_bowtie_to_blast(self):
        # in case no hits are returned
        try:
            df = pd.read_csv(self.bowtie_output_file, sep="\t", header=None)
        except pd.errors.EmptyDataError:
            df = pd.DataFrame(columns = ["", "", "", "", "", "", "", "", "", "", "", ""])
            
        split_df = self.splitDataFrameIntoChunks(df)

        func = partial(self._convert_bowtie_to_blast)
        with Pool(self.num_threads) as p:
            multiprocessing_results = list(tqdm.tqdm(p.imap(func, split_df), total=len(split_df)))

        self.df_bowtie = pd.DataFrame(
            columns=["qseqid", "sseqid", "nident", "qlen", "length", "mismatch", "qstart", "qend", "sstart", "sseq",
                     "sstrand", "send"])
        self.df_bowtie = pd.concat(multiprocessing_results, ignore_index=True)

    def _convert_bowtie_to_blast(self, df):
        df_bowtie = pd.DataFrame(
            columns=["qseqid", "sseqid", "nident", "qlen", "length", "mismatch", "qstart", "qend", "sstart", "sseq",
                     "sstrand", "send"])

        for index, line in df.iterrows():
            mismatch = str(line[7]).count(":")
            if line[0].endswith("_probe") and mismatch > self.max_probe_mismatches:
                continue
            if not line[0].endswith("_probe") and mismatch > self.max_primer_mismatches:
                continue
            sstrand = "plus" if line[1] == "+" else "minus"
            qseqid = line[0]
            sseqid = line[2]
            qlen = len(line[4])
            length = qlen
            qstart = 1
            qend = qlen
            sstart = int(line[3]) + 1
            send = sstart + qlen - 1
            nident = qlen - mismatch

            if sstrand == "minus":
                temp_swap = sstart
                sstart = send
                send = temp_swap
                if mismatch == 0:
                    sseq = str(Seq.Seq(line[4]).reverse_complement())
                else:
                    sseq = self._resolve_bowtie_mismtches(line[4], line[7], -1)
            else:
                if mismatch == 0:
                    sseq = line[4]
                else:
                    sseq = self._resolve_bowtie_mismtches(line[4], line[7], 1)

            df_bowtie.loc[len(df_bowtie)] = [str(qseqid), str(sseqid), str(nident), str(qlen), str(length),
                                             str(mismatch), str(qstart),
                                             str(qend), str(sstart), str(sseq), str(sstrand), str(send)]
        return df_bowtie

    def _resolve_bowtie_mismtches(self, sequence, mismatches, strand):
        sequence = Seq.Seq(sequence) if strand == 1 else Seq.Seq(sequence).reverse_complement()
        mismatches = mismatches.split(",")
        for mismatch in mismatches:
            position, base = mismatch.split(":", 1)
            position = int(position)
            base = base[0] if strand == 1 else Seq.Seq(base[0]).reverse_complement()
            sequence = sequence[:position] + base + sequence[position+1:]
        return str(sequence)

    def _split_output(self):
        if self.method == "blast":
            # in case no hits are returned
            try:
                df_blast = pd.read_csv(self.blast_output_file, sep="\t", header=None)
            except pd.errors.EmptyDataError:
                df_blast = pd.DataFrame(columns = ["", "", "", "", "", "", "", "", "", "", "", ""])
            self.df_blast_split = self.splitDataFrameIntoChunks(df_blast)
        if self.method == "bowtie":
            self.df_bowtie_split = self.splitDataFrameIntoChunks(self.df_bowtie)
        if self.method == "aho-corasick":
            self.df_aho_split = self.splitDataFrameIntoChunks(self.df_aho)

    def splitDataFrameIntoChunks(self, df):
        chunkSize = math.ceil(len(df) / self.num_threads)
        if chunkSize == 0:
            chunkSize = 1
        listOfDf = list()
        numberChunks = len(df) // chunkSize + 1
        for i in range(numberChunks):
            listOfDf.append(df[i * chunkSize:(i + 1) * chunkSize])
        return listOfDf

    def _multiprocess_split_files(self):
        if self.method == "blast":
            input_files = self.df_blast_split
        if self.method == "bowtie":
            input_files = self.df_bowtie_split
        if self.method == "aho-corasick":
            input_files = self.df_aho_split

        func = partial(self._parse_blastlike_results_df)
        with Pool(self.num_threads) as p:
            multiprocessing_results = list(tqdm.tqdm(p.imap(func, input_files), total=len(input_files)))

        self.hits = pd.concat(multiprocessing_results, ignore_index=True)

        hits_output = self.hits.copy()
        if len(hits_output.index) > 0:
            hits_output[['Sequence', 'Contig']] = hits_output['Sequence'].str.split("__contigname__", 1, expand=True)
            hits_output = hits_output[
                ['Sequence', 'Contig', 'Type', 'Name', 'Package', 'StartPos', 'EndPos', 'MismatchesTotal', 'Strand',
                 'HitSequence', 'Tm', 'dG']]
            tmp = hits_output['Name'].str.rsplit("_probe", 1, expand = True)
            hits_output['Name'] = tmp[0]
        hits_output.to_csv(self.output_folder + "all_hits.csv", index=False, sep=";")

    def _process_probes_only(self):
        probes_df = self.hits[(self.hits['Type'] == "Probe")]
        if len(probes_df.index) > 0:
            oligos_full_sequences = SeqIO.index(self.oligo_file, "fasta")

            probes_df = probes_df.drop(columns = ['Type', 'Strand'])
            probes_df = probes_df.rename(columns = {'Name': 'Probe', 'Package': 'ProbePackage', 'MismatchesTotal': 'NoMismatchesProbe'})
            probes_df[['Sequence', 'Contig']] = probes_df['Sequence'].str.split("__contigname__", 1, expand = True)
            probes_df['MismatchesProbe'] = probes_df.apply(lambda x: self._mismatch_visualization(oligos_full_sequences[x['ProbePackage'] + "^" + x['Probe']].seq, x['HitSequence']), axis=1)
            probes_df = probes_df[['Sequence', 'Contig', 'Probe', 'ProbePackage', 'StartPos', 'EndPos', 'NoMismatchesProbe', 'MismatchesProbe', 'HitSequence', 'Tm' ,'dG']]
            tmp = probes_df['Probe'].str.rsplit("_probe", 1, expand = True)
            probes_df['Probe'] = tmp[0]
        probes_df.to_csv(self.results_all, index=False, sep=";")

        # parse wobbled primers
        subset = probes_df[probes_df['Probe'].str.contains("_mut")]
        subset_r = subset.replace(['_mut([0-9])+'], [''], regex=True)

        # hits without mutations
        unique = probes_df.merge(subset, indicator=True, how="outer")
        unique = unique[unique['_merge'] == 'left_only']
        unique = unique.drop("_merge", axis=1)

        results2 = pd.DataFrame(columns=['Sequence', 'Contig', 'Probe', 'ProbePackage', 'StartPos', 'EndPos', 'NoMismatchesProbe', 'MismatchesProbe', 'HitSequence', 'Tm', 'dG'])

        for s in subset_r.groupby(['Sequence', 'Contig', 'Probe', 'ProbePackage', 'StartPos', 'EndPos', 'HitSequence']).groups.items():
            # Fields to be changed: NoMismatchesProbe, MismatchesProbe
            sample = subset_r.loc[s[1]]  # get one set
            first_row = sample.iloc[0]
            if len(sample) < 2:
                results2.loc[len(results2)] = first_row

            else:
                mismatch_min, mismatch_max = min(sample['NoMismatchesProbe']), max(sample['NoMismatchesProbe'])
                mismatches = mismatch_min if mismatch_min == mismatch_max else str(mismatch_min) + "-" + str(mismatch_max)
                
                tm_min, tm_max = min(sample['Tm']), max(sample['Tm'])
                tm = tm_min if tm_min == tm_max else str(tm_min) + "-" + str(tm_max)
                
                dg_min, dg_max = min(sample['dG']), max(sample['dG'])
                dg = dg_min if dg_min == dg_max else str(dg_min) + "-" + str(dg_max)

                # Get first row and then replace values of first row (all other fields are identifical)
                results2.loc[len(results2)] = [first_row['Sequence'], first_row['Contig'], first_row['Probe'], first_row['ProbePackage'],
                                               first_row['StartPos'], first_row['EndPos'], mismatches, '', first_row['HitSequence'], tm, dg]
        wobbled = unique.append(results2)
        wobbled.to_csv(self.results_wob, index=False, sep=";")

    def _parse_blastlike_results_df(self, blast_df):
        oligos_full_sequences = SeqIO.index(self.oligo_file, "fasta")
        hits = pd.DataFrame(columns=self.pd_col_hits)
        for index, line in blast_df.iterrows():
            if self.method == "aho-corasick":
                if line[0].endswith("_probe") and int(line[5]) > self.max_probe_mismatches:
                    continue
                if not line[0].endswith("_probe") and int(line[5]) > self.max_primer_mismatches:
                    continue
                new_package = line[0].split("^", 1)[0]
                new_qresult = line[0].split("^", 1)[1]
                hit_strand = 1 if line[10] == "plus" else -1
                mismatches_total = int(line[5])
                hit_seq = line[9]
                type = "Probe" if line[0].endswith("_probe") == True else "Primer"

                if hit_strand == -1:
                    temp_swap = int(line[8])
                    sstart = int(line[11])
                    send = temp_swap
                else:
                    sstart = int(line[8])
                    send = int(line[11])

                tm, dg = self._calc_thermal_parameters(str(oligos_full_sequences[line[0]].seq.reverse_complement()), hit_seq, type)

                hits.loc[len(hits)] = [line[1], type, new_qresult, new_package, sstart, send,
                                       mismatches_total, hit_strand, hit_seq, tm, dg]
            else:
                mismatches_left = int(line[6])
                mismatches_right = int(line[3]) - int(line[7])
                mismatches_middle = int(line[3]) - int(line[2]) - mismatches_left - mismatches_right
                mismatches_total = mismatches_left + mismatches_right + mismatches_middle

                if line[0].endswith("_probe") and mismatches_total > self.max_probe_mismatches:
                    continue
                if not line[0].endswith("_probe") and mismatches_total > self.max_primer_mismatches:
                    continue

                new_package = line[0].split("^", 1)[0]
                new_qresult = line[0].split("^", 1)[1]
                hit_strand = 1 if line[10] == "plus" else -1
                type = "Probe" if line[0].endswith("_probe") == True else "Primer"

                correct_start = mismatches_left - 1 if hit_strand == 1 else mismatches_right if hit_strand == -1 else 0
                correct_end = mismatches_right if hit_strand == 1 else mismatches_left - 1 if hit_strand == -1 else 0

                if hit_strand == -1:
                    temp_swap = int(line[8])
                    sstart = int(line[11]) - correct_start
                    send = temp_swap + correct_end
                else:
                    sstart = int(line[8]) - correct_start
                    send = int(line[11]) + correct_end

                if mismatches_left > 0 or mismatches_right > 0:
                    hit_seq = self._get_sequence(self.output_contigs, line[1], sstart, send, hit_strand)
                else:
                    hit_seq = line[9]

                tm, dg = self._calc_thermal_parameters(str(oligos_full_sequences[line[0]].seq.reverse_complement()), hit_seq, type)

                hits.loc[len(hits)] = [line[1], type, new_qresult, new_package, sstart, send,
                                       mismatches_total, hit_strand, hit_seq, tm, dg]

        return hits

    def _multiprocess_hits(self):
        objects = []
        groups = []
        num_groups = math.ceil(len(self.hits[self.hits['Type'] == "Primer"].groupby("Sequence").groups.items()) / self.num_threads)
        i = 1
        for s in self.hits[self.hits['Type'] == "Primer"].groupby("Sequence").groups.items():
            if i > num_groups:
                groups.append(objects)
                objects = []
                i = 1
            objects.append(self.hits.loc[s[1]])
            i += 1
        groups.append(objects)

        multiprocessing_results = []
        func = partial(self._parse_hits)
        with Pool(self.num_threads) as p:
            multiprocessing_results = list(tqdm.tqdm(p.imap(func, groups), total=len(groups)))
        self.results = pd.concat(multiprocessing_results, ignore_index=True)

    def _parse_hits(self, groups):
        oligos_full_sequences = SeqIO.index(self.oligo_file, "fasta")
        results = pd.DataFrame(columns=self.pd_col_results)
        for df in groups:
            start = df['StartPos'].values.tolist()
            end = df['EndPos'].values.tolist()
            a = [(y - x, (i, j)) for i, x in enumerate(start) for j, y in enumerate(end) if
                 (y - x) > 0 and (y - x) <= self.max_product_size and i != j]
            if len(a) == 0:
                continue
            for product_size, (left, right) in a:
                if df['Strand'].iloc[left] == df['Strand'].iloc[right] or df['Strand'].iloc[left] == -1:
                    continue
                if self.same_package == True and df['Package'].iloc[left] != df['Package'].iloc[right]:
                    continue
                sequence, contig = df['Sequence'].iloc[left].split("__contigname__", 1)
                product = self._get_sequence(self.output_contigs, df['Sequence'].iloc[left], df['StartPos'].iloc[left],
                                             df['EndPos'].iloc[right], 1)
                tmp_left_name = df["Package"].iloc[left] + "^" + df['Name'].iloc[left]
                mismatches1 = self._mismatch_visualization(oligos_full_sequences[tmp_left_name].seq,
                                                           df['HitSequence'].iloc[left])
                tmp_right_name = df["Package"].iloc[right] + "^" + df['Name'].iloc[right]
                mismatches2 = self._mismatch_visualization(oligos_full_sequences[tmp_right_name].seq,
                                                           df['HitSequence'].iloc[right])

                probe_matches = self.hits[(self.hits['Type'] == "Probe") &
                                          (self.hits['Sequence'] == df['Sequence'].iloc[left]) &
                                          (self.hits['StartPos'] >= df['StartPos'].iloc[left]) &
                                          (self.hits['EndPos'] <= df['EndPos'].iloc[right])]
                if self.same_package == True:
                    probe_package_match = str(df['Package'].iloc[left]) + "_probes"
                    probe_matches = probe_matches[(probe_matches['Package'] == probe_package_match)]
                no_probe_hits = len(probe_matches.index)

                product_tm = self._calc_Tm(str(product), "Product")

                if no_probe_hits == 0:
                    # save matching primers, no matching probe found
                    results.loc[len(results)] = [sequence, contig, df['Name'].iloc[left], df['Name'].iloc[right], "",
                                                 df['Package'].iloc[left], df['Package'].iloc[right], "",
                                                 df['StartPos'].iloc[left], df['EndPos'].iloc[left],
                                                 df['StartPos'].iloc[right], df['EndPos'].iloc[right], "", "",
                                                 df['Tm'].iloc[left], df['Tm'].iloc[right], "", df['dG'].iloc[left], df['dG'].iloc[right], "", 
                                                 int(product_size+1), product_tm, df['MismatchesTotal'].iloc[left],
                                                 df['MismatchesTotal'].iloc[right], "", mismatches1, mismatches2, "",
                                                 "", str(product)]
                else:
                    for index, row in probe_matches.iterrows():
                        tmp_probe_name = row['Package'] + "^" + row['Name']
                        probe_mivi = self._mismatch_visualization(oligos_full_sequences[tmp_probe_name].seq,
                                                                  row['HitSequence'])
                        comment = "More than 1 probe binding to amplicon." if no_probe_hits > 1 else ""
                        # save matching primers and matching probes, multiple probes per primer pair possible
                        results.loc[len(results)] = [sequence, contig, df['Name'].iloc[left], df['Name'].iloc[right],
                                                     row['Name'], df['Package'].iloc[left], df['Package'].iloc[right],
                                                     row['Package'], df['StartPos'].iloc[left], df['EndPos'].iloc[left],
                                                     df['StartPos'].iloc[right], df['EndPos'].iloc[right],
                                                     row['StartPos'], row['EndPos'], 
                                                     df['Tm'].iloc[left], df['Tm'].iloc[right], row['Tm'], df['dG'].iloc[left], df['dG'].iloc[right], row['dG'],
                                                     int(product_size+1), product_tm, 
                                                     df['MismatchesTotal'].iloc[left],
                                                     df['MismatchesTotal'].iloc[right], row['MismatchesTotal'],
                                                     mismatches1, mismatches2, probe_mivi, comment, str(product)]

        return results

    def _parse_results_to_wobbled(self):
        subset = self.results[
            self.results['Primer1'].str.contains("_mut") | self.results['Primer2'].str.contains("_mut")]
        subset_r = subset.replace(['_mut([0-9])+'], [''], regex=True)

        # results without mutations
        unique = self.results.merge(subset, indicator=True, how="outer")
        unique = unique[unique['_merge'] == 'left_only']
        unique = unique.drop("_merge", axis=1)

        results2 = pd.DataFrame(columns=self.pd_col_results)

        for s in subset_r.groupby(
                ['Sequence', 'Contig', 'Primer1', 'Primer2', 'Probe', 'Primer1Package', 'Primer2Package',
                 'ProbePackage',
                 'StartPos1', 'EndPos1', 'StartPos2', 'EndPos2',
                 'StartPos3', 'EndPos3', 'ProductSize', 'Comment',
                 'Product']).groups.items():
            # Fields to be changed: NoMismatchesLeft, NoMismatchesRight, MismatchesLeft, MismatchesRight, NoMismatchesProbe, MismatchesProbe
            sample = subset_r.loc[s[1]]  # get one set
            first_row = sample.iloc[0]

            if len(sample) < 2:
                results2.loc[len(results2)] = first_row

            else:
                primer_left_min, primer_left_max = min(sample['NoMismatchesLeft']), max(sample['NoMismatchesLeft'])
                primer_left = primer_left_min if primer_left_min == primer_left_max else str(primer_left_min) + "-" + str(primer_left_max)
                primer_right_min, primer_right_max = min(sample['NoMismatchesRight']), max(sample['NoMismatchesRight'])
                primer_right = primer_right_min if primer_right_min == primer_right_max else str(primer_right_min) + "-" + str(primer_right_max)
                probe_min, probe_max = min(sample['NoMismatchesProbe']), max(sample['NoMismatchesProbe'])
                probe = probe_min if probe_min == probe_max else str(probe_min) + "-" + str(probe_max)

                primer1tm_min, primer1tm_max = min(sample['Primer1Tm']), max(sample['Primer1Tm'])
                primer1tm = primer1tm_min if primer1tm_min == primer1tm_max else str(primer1tm_min) + "-" + str(primer1tm_max)
                primer2tm_min, primer2tm_max = min(sample['Primer2Tm']), max(sample['Primer2Tm'])
                primer2tm = primer2tm_min if primer2tm_min == primer2tm_max else str(primer2tm_min) + "-" + str(primer2tm_max)
                primer1dg_min, primer1dg_max = min(sample['Primer1dG']), max(sample['Primer1dG'])
                primer1dg = primer1dg_min if primer1dg_min == primer1dg_max else str(primer1dg_min) + "-" + str(primer1dg_max)
                primer2dg_min, primer2dg_max = min(sample['Primer2dG']), max(sample['Primer2dG'])
                primer2dg = primer2dg_min if primer2dg_min == primer2dg_max else str(primer2dg_min) + "-" + str(primer2dg_max)
                producttm_min, producttm_max = min(sample['ProductTm']), max(sample['ProductTm'])
                producttm = producttm_min if producttm_min == producttm_max else str(producttm_min) + "-" + str(producttm_max)
                probetm_min, probetm_max = min(sample['ProbeTm']), max(sample['ProbeTm'])
                probetm = probetm_min if probetm_min == probetm_max else str(probetm_min) + "-" + str(probetm_max)
                probedg_min, probedg_max = min(sample['ProbedG']), max(sample['ProbedG'])
                probedg = probedg_min if probedg_min == probedg_max else str(probedg_min) + "-" + str(probedg_max)

                # Get first row and then replace values of first row (all other fields are identifical)
                results2.loc[len(results2)] = [first_row['Sequence'], first_row['Contig'], first_row['Primer1'],
                                               first_row['Primer2'], first_row['Probe'],
                                               first_row['Primer1Package'], first_row['Primer2Package'],
                                               first_row['ProbePackage'],
                                               first_row['StartPos1'], first_row['EndPos1'],
                                               first_row['StartPos2'], first_row['EndPos2'], first_row['StartPos3'],
                                               first_row['EndPos3'], primer1tm, primer2tm, probetm, primer1dg, primer2dg, probedg, first_row['ProductSize'], producttm, primer_left,
                                               primer_right, probe, '', '',
                                               first_row['MismatchesProbe'], first_row['Comment'], first_row['Product']]
        self.wobbled = unique.append(results2)

    def _save_results(self):
        if len(self.results.index) > 0:
            tmp = self.results['Probe'].str.rsplit("_probe", 1, expand = True)
            self.results['Probe'] = tmp[0]
        self.results.to_csv(self.results_all, index=False, sep=";")

    def _save_wobbled_results(self):
        self.wobbled.to_csv(self.results_wob, index=False, sep=";")

    def _add_oligos_to_automaton(self):
        sequences = SeqIO.parse(open(self.oligo_file), "fasta")
        for fasta in sequences:
            name = str(fasta.id)
            max_mismatches = self.max_probe_mismatches if name.endswith("_probe") else self.max_primer_mismatches
            i = 0
            while i <= max_mismatches:
                for count, elem in enumerate(self._generate_primer_mismatches(fasta.seq, i)):
                    self.aho.add_word(elem.upper(), (name + "__" + str(count), elem.upper()))
                for count, elem in enumerate(self._generate_primer_mismatches(fasta.seq.reverse_complement(), i)):
                    self.aho.add_word(elem.upper(), (name + "__rc__" + str(count), elem.upper()))
                i += 1

    def _generate_primer_mismatches(self, s, d=1):
        N = len(s)
        letters = 'ACGT'
        pool = list(s)

        for indices in combinations(range(N), d):
            for replacements in product(letters, repeat=d):
                skip = False
                for i, a in zip(indices, replacements):
                    if pool[i] == a: skip = True
                if skip: continue

                keys = dict(zip(indices, replacements))
                yield ''.join([pool[i] if i not in indices else keys[i] for i in range(N)])

    def _multiprocess_string_search(self):
        objects = []
        groups = []
        num_seq = sum(1 for x in SeqIO.parse(open(self.output_contigs), "fasta"))
        num_groups = math.ceil(num_seq / self.num_threads)
        # dynamic group allocation does not always work
        # num_groups = 256
        i = 1
        for fasta in SeqIO.parse(open(self.output_contigs), "fasta"):
            if i > num_groups:
                groups.append(objects)
                objects = []
                i = 1
            objects.append([str(fasta.seq).upper(), str(fasta.id)])
            i += 1
        groups.append(objects)

        multiprocessing_results = []
        func = partial(self._string_search_match)
        with Pool(self.num_threads) as p:
            multiprocessing_results = list(tqdm.tqdm(p.imap(func, groups), total=len(groups)))

        self.df_aho = pd.concat(multiprocessing_results, ignore_index=True)
        self.df_aho = self.df_aho.sort_values(['qseqid', 'mismatch'])

    def _string_search_match(self, groups):
        aho = self.aho
        aho.make_automaton()
        oligos_full_sequences = SeqIO.index(self.oligo_file, "fasta")
        df = pd.DataFrame(
            columns=["qseqid", "sseqid", "nident", "qlen", "length", "mismatch", "qstart", "qend", "sstart", "sseq",
                     "sstrand", "send"])

        for object in groups:
            haystack = object[0]
            fasta = object[1]
            for end_pos, (mutated_oligo_name, mutated_oligo_sequence) in aho.iter(haystack):
                print("mutated oligo sequence: " + str(mutated_oligo_sequence))
                start_pos = end_pos - len(mutated_oligo_sequence)
                length = len(mutated_oligo_sequence)
                if "__rc__" in mutated_oligo_name:
                    oligo_name = mutated_oligo_name.rsplit("__rc__", 1)[0]
                    oligo_seq = str(oligos_full_sequences[oligo_name].seq).upper()
                    strand = "minus"
                    start = end_pos + 1
                    end = start_pos + 2
                    sseq = self._get_sequence(self.output_contigs, str(fasta), end, start, -1)
                    mismatches = self._number_of_mismatches(sseq, oligo_seq)
                    nident = length - mismatches
                else:
                    oligo_name = mutated_oligo_name.rsplit("__", 1)[0]
                    oligo_seq = str(oligos_full_sequences[oligo_name].seq).upper()
                    strand = "plus"
                    start = start_pos + 2
                    end = end_pos + 1
                    sseq = self._get_sequence(self.output_contigs, str(fasta), start, end, 1)
                    mismatches = self._number_of_mismatches(sseq, oligo_seq)
                    nident = length - mismatches
                df.loc[len(df)] = [str(oligo_name), str(fasta), str(nident), str(length), str(length), str(mismatches),
                                   str("1"), str(length), str(start), str(sseq), str(strand), str(end)]
        return df

    def _number_of_mismatches(self, s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    def _run_cross_dimer_check(self):
        for package in self.primer_packages:
            oligo_list = []
            for oligo_a in self.primer_packages[package]:
                oligo_list.append((package, oligo_a, oligo_a, self.primer_packages[package][oligo_a], self.primer_packages[package][oligo_a]))
            for oligo_a, oligo_b in combinations(self.primer_packages[package], 2):
                oligo_list.append((package, oligo_a, oligo_b, self.primer_packages[package][oligo_a], self.primer_packages[package][oligo_b]))

            pool = Pool(self.num_threads)
            multiprocessing_results = []
            func = partial(self._cross_dimer_check)
            multiprocessing_results = list(tqdm.tqdm(pool.map(func, oligo_list), total=len(oligo_list)))
            pool.close()
            pool.join()

            df = pd.concat(multiprocessing_results, ignore_index=True)
            df.name = "CD_" + package
            self.cross_dimer_dfs.append(df)

    def _cross_dimer_check(self, seq_tuple):
        package, oligo_a, oligo_b, seq1, seq2 = seq_tuple
        try:
            tr = primer3.calcHeterodimer(seq1, seq2, mv_conc = self.primer_monovalent_cations, dv_conc = self.primer_divalent_cations, dntp_conc = self.primer_dntps, dna_conc = self.primer_annealing_oligo, temp_c = self.primer_annealing_temp)
            deltaG = tr.dg
            Tm = tr.tm
        except:
            deltaG = "NaN"
            Tm = "NaN"
        df = pd.DataFrame(columns=['Package', 'Primer1', 'Primer2', 'dG', 'Tm'])
        df.loc[len(df)] = [str(package), str(oligo_a), str(oligo_b), str(deltaG), str(Tm)]
        return df

    def _calc_thermal_parameters(self, seq1, seq2, seqtype):
        if seqtype == "Primer" or seqtype == "Product":
            mv, dv, dntp, dna, temp = self.primer_monovalent_cations, self.primer_divalent_cations, self.primer_dntps, self.primer_annealing_oligo, self.primer_annealing_temp
        else:
            mv, dv, dntp, dna, temp = self.probe_monovalent_cations, self.probe_divalent_cations, self.probe_dntps, self.probe_annealing_oligo, self.probe_annealing_temp
        try:
            tr = primer3.calcHeterodimer(str(seq1), str(seq2), mv_conc = mv, dv_conc = dv, dntp_conc = dntp, dna_conc = dna, temp_c = temp)
            deltaG = round(tr.dg / 1000, 2)
            Tm = round(tr.tm, 2)
        except:
            deltaG = "NaN"
            Tm = "NaN"
        return Tm, deltaG

    def _calc_Tm(self, seq, seqtype):
        if seqtype == "Primer" or seqtype == "Product":
            mv, dv, dntp, dna = self.primer_monovalent_cations, self.primer_divalent_cations, self.primer_dntps, self.primer_annealing_oligo
        else:
            mv, dv, dntp, dna = self.probe_monovalent_cations, self.probe_divalent_cations, self.probe_dntps, self.probe_annealing_oligo
        try:
            tm = primer3.calcTm(str(seq), mv_conc = mv, dv_conc = dv, dntp_conc = dntp, dna_conc = dna)
            tm = round(tm, 2)
        except:
            tm = "NaN"
        return tm

    def _run_hairpin_check(self):
        for package in self.primer_packages:
            oligo_list = []
            for oligo_name in self.primer_packages[package]:
                oligo_list.append((package, oligo_name, self.primer_packages[package][oligo_name], self.primer_monovalent_cations,
                                   self.primer_divalent_cations, self.primer_dntps, self.primer_annealing_oligo, self.primer_annealing_temp))
            func = partial(self._hairpin_check)
            with Pool(self.num_threads) as p:
                multiprocessing_results = list(tqdm.tqdm(p.imap(func, oligo_list), total=len(oligo_list)))
            df = pd.concat(multiprocessing_results, ignore_index=True)
            df.name = "HP_" + package
            self.hairpin_dfs.append(df)
        for package in self.probe_packages:
            oligo_list = []
            for oligo_name in self.probe_packages[package]:
                oligo_list.append((package, oligo_name, self.probe_packages[package][oligo_name], self.probe_monovalent_cations,
                                   self.probe_divalent_cations, self.probe_dntps, self.probe_annealing_oligo, self.probe_annealing_temp))
            func = partial(self._hairpin_check)
            with Pool(self.num_threads) as p:
                multiprocessing_results = list(tqdm.tqdm(p.imap(func, oligo_list), total=len(oligo_list)))
            df = pd.concat(multiprocessing_results, ignore_index=True)
            df.name = "HP_" + package
            self.hairpin_dfs.append(df)

    def _hairpin_check(self, seq_tuple):
        package, oligo, seq, mv, dv, dntps, annealing_oligo, annealing_temp = seq_tuple
        try:
            tr = primer3.calcHairpin(seq, mv_conc = mv, dv_conc = dv, dntp_conc = dntps, dna_conc = annealing_oligo, temp_c = annealing_temp)
            deltaG = tr.dg
            Tm = tr.tm
        except:
            deltaG = "NaN"
            Tm = "NaN"
        df = pd.DataFrame(columns=['Package', 'Oligo', 'dG', 'Tm'])
        df.loc[len(df)] = [str(package), str(oligo), str(deltaG), str(Tm)]
        return df

    def _analyze_cross_dimer_dataframe(self):
        for count, df in enumerate(self.cross_dimer_dfs):
            # only wobbled primers in wob_r
            wob = df[df['Primer1'].str.contains("_mut") | df['Primer2'].str.contains("_mut")]
            wob_r = wob.replace(['_mut[0-9]'], [''], regex=True)

            # all primeres except wobbled primers in unique
            unique = df.merge(wob, indicator=True, how="outer")
            unique = unique[unique['_merge'] == 'left_only']
            unique = unique.drop("_merge", axis=1)

            # remove unneeded columns for current analysis
            wob_r = wob_r.drop(columns=["Package", "Tm"])
            unique = unique.drop(columns=["Package", "Tm"])

            # empty output dataframe
            new_df = pd.DataFrame(data=None, columns=wob_r.columns)

            # go through wobbled sets
            for s in wob_r.groupby(['Primer1', 'Primer2']).groups.items():
                sample = wob_r.loc[s[1]]
                sample["dG"] = pd.to_numeric(sample["dG"])
                new_df.loc[len(new_df)] = sample.loc[sample['dG'].idxmin()]

            # append unique sets
            new_df = new_df.append(unique)
            new_df["dG"] = pd.to_numeric(new_df["dG"])

            # pivot table output
            pivot = pd.pivot_table(new_df, values = "dG", index=["Primer1"], columns="Primer2").reset_index()

            # reorder pivot table
            pivot['Total'] = pivot.count(axis=1)
            pivot = pivot.sort_values(by='Total', ascending = False)
            pivot = pivot.drop(columns=["Total"])
            pivot.loc['Total']= pivot.count()
            pivot = pivot.sort_values(by='Total', axis = 1, ascending = False)
            pivot = pivot.drop(index=["Total"])

            # replace df in array
            pivot.name = df.name
            self.cross_dimer_dfs_dg.append(pivot)

    def _save_dimer_checks(self):
        writer = pd.ExcelWriter(self.results_dimers, engine='xlsxwriter')
        for df in self.cross_dimer_dfs_dg:
            sheet_name = df.name[:31] if len(df.name) > 31 else df.name
            df.to_excel(writer, sheet_name = sheet_name)
        for df in self.hairpin_dfs:
            sheet_name = df.name[:31] if len(df.name) > 31 else df.name
            df.to_excel(writer, sheet_name = sheet_name)
        writer.save()

    def _check_prebuilt_indexes(self):
        for db in self.prebuilt.split(","):
            if db in self.bowtie_dbs:
                self.bowtie_runs.append(db)

