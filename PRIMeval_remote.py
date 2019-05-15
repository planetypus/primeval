#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from PRIMeval import PRIMeval
import urllib3
import json
import yaml
import os

__author__ = "Rick Conzemius"
__copyright__ = "Copyright 2019, Rick Conzemius"
__license__ = "MIT"
__maintainer__ = "Rick Conzemius"
__version__ = "1.0"
__email__ = "rick.conzemius.fl@ait.ac.at"
__status__ = "Production"

key = "publickey"
update_url   = "http://update_url/"
upload_url   = "http://upload_url/"
run_folder = os.getcwd() + "/runs/"

def primer_run(run_id):
    params_file = run_folder + str(run_id) + "/params.txt"
    with open(params_file, 'r') as stream:
        try:
            params = yaml.load(stream)
        except yaml.YAMLError as exc:
            return False

    primer_mismatches = params['primerMismatchesTotal']
    probe_mismatches = params['probeMismatchesTotal']
    max_product_size = params['productLength']
    cross_check = params['crossCheck']
    dimer_check = params['dimerCheck']
    probes_only = params['checkProbesOnly']
    method = params['method']

    primerMonovalentCations = params['primerMonovalentCations']
    primerDivalentCations = params['primerDivalentCations']
    primerDNTPs = params['primerDNTPs']
    primerConcentration = params['primerConcentration']
    primerAnnealingTemp = params['primerAnnealingTemp']
    probeMonovalentCations = params['probeMonovalentCations']
    probeDivalentCations = params['probeDivalentCations']
    probeDNTPs = params['probeDNTPs']
    probeConcentration = params['probeConcentration']
    probeAnnealingTemp = params['probeAnnealingTemp']

    dbs = params['dbs']

    pe = PRIMevalRun(run_id, primer_mismatches, probe_mismatches, max_product_size, cross_check, probes_only, method, dimer_check,
                     primerMonovalentCations, primerDivalentCations, primerDNTPs, primerConcentration, primerAnnealingTemp,
                     probeMonovalentCations, probeDivalentCations, probeDNTPs, probeConcentration, probeAnnealingTemp, dbs)
    return pe.run_analysis()

def upload_results(runid):
    update_status(runid, "Upload", "Queued", 0, 1)
    update_status(runid, "Upload", "Upload", 1, 0)

    http = urllib3.PoolManager()
    file1 = run_folder + str(runid) + "/output/all_hits.csv"
    file2 = run_folder + str(runid) + "/output/results.csv"
    file3 = run_folder + str(runid) + "/output/results_wobbled.csv"
    file4 = run_folder + str(runid) + "/output/results_dimers.xlsx"
    with open(file1) as fp1:
        file_data1 = fp1.read()
    with open(file2) as fp2:
        file_data2 = fp2.read()
    with open(file3) as fp3:
        file_data3 = fp3.read()
    if os.path.isfile(file4):
        with open(file4, "rb") as fp4:
            file_data4 = fp4.read()
    try:
        if os.path.isfile(file4):
            r = http.request("POST", upload_url + str(runid), fields={'hits': ('all_hits.csv', file_data1),
                                                                      'results_all': ('results.csv', file_data2),
                                                                      'results_wob': ('results_wobbled.csv', file_data3),
                                                                      'results_dimer': ('results_dimers.xlsx', file_data4),
                                                                      'runid': runid, })
        else:
            r = http.request("POST", upload_url + str(runid), fields={'hits': ('all_hits.csv', file_data1),
                                                                      'results_all': ('results.csv', file_data2),
                                                                      'results_wob': ('results_wobbled.csv', file_data3),
                                                                      'runid': runid, })
        success = json.loads(r.data.decode("utf-8"))
        if success['upload_is_valid'] == False:
            update_status(runid, "Upload", "Upload", 0, 0)
            return False
    except:
        update_status(runid, "Upload", "Upload", 0, 0)
        return False
    update_status(runid, "Upload", "Upload", 0, 1)
    return True

def update_status(run_id, step, substep, running, finished, error_msg = ""):
    http = urllib3.PoolManager()
    http.request('GET', update_url + str(run_id) + "/" + str(step) + "/" + str(substep) + "/" + str(running) + "/" + str(finished) + "/" + error_msg, retries=10)
    print("Status updated for step " + str(step) + "/" + str(substep) + " for run " + str(run_id) + ".")

class PRIMevalRun(PRIMeval):

    def run_analysis(self):
        self._prepare_folders()
        update_status(self.run_id, "PRIMEval", "PreparingFolders", 0, 1)
        update_status(self.run_id, "PRIMEval", "ImportingSequences", 1, 0)
        self._import_sequences()
        update_status(self.run_id, "PRIMEval", "ImportingSequences", 0, 1)
        if self.dimer_check == True:
            update_status(self.run_id, "PRIMEval", "DimerCheck", 1, 0)
            self._run_cross_dimer_check()
            self._run_hairpin_check()
            self._analyze_cross_dimer_dataframe()
            self._save_dimer_checks()
            update_status(self.run_id, "PRIMEval", "DimerCheck", 0, 1)
        if self.method == "blast":
            update_status(self.run_id, "PRIMEval", "CreatingBLASTDB", 1, 0)
            self._create_blast_db()
            update_status(self.run_id, "PRIMEval", "CreatingBLASTDB", 0, 1)
            update_status(self.run_id, "PRIMEval", "RunningBLAST", 1, 0)
            self._blast_call()
            update_status(self.run_id, "PRIMEval", "RunningBLAST", 0, 1)
        if self.method == "bowtie":
            update_status(self.run_id, "PRIMEval", "CreatingBowtieIndex", 1, 0)
            self._create_bowtie_index()
            update_status(self.run_id, "PRIMEval", "CreatingBowtieIndex", 0, 1)
            update_status(self.run_id, "PRIMEval", "RunningBowtie", 1, 0)
            self._bowtie_call()
            self._check_prebuilt_indexes()
            self._specificity_calls()
            update_status(self.run_id, "PRIMEval", "RunningBowtie", 0, 1)
            update_status(self.run_id, "PRIMEval", "ConvertingBowtieHits", 1, 0)
            self._multiprocess_convert_bowtie_to_blast()
            update_status(self.run_id, "PRIMEval", "ConvertingBowtieHits", 0, 1)
        if self.method == "aho-corasick":
            update_status(self.run_id, "PRIMEval", "AddOligosToAhoCorasickAutomaton", 1, 0)
            self._add_oligos_to_automaton()
            update_status(self.run_id, "PRIMEval", "AddOligosToAhoCorasickAutomaton", 0, 1)
            update_status(self.run_id, "PRIMEval", "RunningAhoCorasickStringSearch", 1, 0)
            self._multiprocess_string_search()
            update_status(self.run_id, "PRIMEval", "RunningAhoCorasickStringSearch", 0, 1)
        update_status(self.run_id, "PRIMEval", "SplittingOutputFiles", 1, 0)
        self._split_output()
        update_status(self.run_id, "PRIMEval", "SplittingOutputFiles", 0, 1)
        update_status(self.run_id, "PRIMEval", "ProcessingOutputFiles", 1, 0)
        self._multiprocess_split_files()
        update_status(self.run_id, "PRIMEval", "ProcessingOutputFiles", 0, 1)
        if self.probes_only == True:
            update_status(self.run_id, "PRIMEval", "ProcessingProbesOnly", 1, 0)
            self._process_probes_only()
            update_status(self.run_id, "PRIMEval", "ProcessingProbesOnly", 0, 1)
        else:
            update_status(self.run_id, "PRIMEval", "ProcessingHits", 1, 0)
            self._multiprocess_hits()
            update_status(self.run_id, "PRIMEval", "ProcessingHits", 0, 1)
            update_status(self.run_id, "PRIMEval", "SaveAllResults", 1, 0)
            self._save_results()
            update_status(self.run_id, "PRIMEval", "SaveAllResults", 0, 1)
            update_status(self.run_id, "PRIMEval", "ParseWobbledResults", 1, 0)
            self._parse_results_to_wobbled()
            update_status(self.run_id, "PRIMEval", "ParseWobbledResults", 0, 1)
            update_status(self.run_id, "PRIMEval", "SaveWobbledResults", 1, 0)
            self._save_wobbled_results()
            update_status(self.run_id, "PRIMEval", "SaveWobbledResults", 0, 1)
            update_status(self.run_id, "PRIMEval", "CleanUpFolders", 1, 0)
            self._clean_up_folders()
            update_status(self.run_id, "PRIMEval", "CleanUpFolders", 0, 1)
            update_status(self.run_id, "PRIMEval", "Run", 0, 1)

        return True
