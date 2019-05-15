#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
from resource import getrusage as resource_usage, RUSAGE_SELF
from time import time as timestamp

from PRIMeval import PRIMeval

__author__ = "Rick Conzemius"
__copyright__ = "Copyright 2019, Rick Conzemius"
__license__ = "MIT"
__maintainer__ = "Rick Conzemius"
__version__ = "1.0"
__email__ = "rick.conzemius.fl@ait.ac.at"
__status__ = "Production"

class PRIMevalRun(PRIMeval):

    def _unix_runtime(self, function, args=tuple(), kwargs={}):
        start_time, start_resources = timestamp(), resource_usage(RUSAGE_SELF)
        function(*args, **kwargs)
        end_resources, end_time = resource_usage(RUSAGE_SELF), timestamp()

    def run_analysis(self):
        if self.method == "blast":
            self._run_blast_analysis()
        elif self.method == "bowtie":
            self._run_bowtie_analysis()
        elif self.method == "aho-corasick":
            self._run_string_analysis()
        else:
            sys.exit("Invalid method.")

    def _run_blast_analysis(self):
        print("1. Prepare folders")
        self._unix_runtime(self._prepare_folders)
        print("2. Import primers, probes and contigs")
        self._unix_runtime(self._import_sequences)
        print("3. Create BLAST database of contigs")
        self._unix_runtime(self._create_blast_db)
        print("4. Run BLASTn using primers and probes as query")
        self._unix_runtime(self._blast_call)
        print("5. Split BLAST tab output file by query")
        self._unix_runtime(self._split_output)
        print("6. Parse BLAST tab output files for meaningful hits using multiprocessing, then merge the results")
        self._unix_runtime(self._multiprocess_split_files)
        print("7. Compare primers and find matching probes")
        self._unix_runtime(self._multiprocess_hits)
        print("\n8. Generate complete CSV output")
        self._unix_runtime(self._save_results)
        print("9. Process results for wobbled primers")
        self._unix_runtime(self._parse_results_to_wobbled)
        print("10. Generate output with wobbled oligos")
        self._unix_runtime(self._save_wobbled_results)
        print("11. Tidying up")
        self._unix_runtime(self._clean_up_folders)

    def _run_bowtie_analysis(self):
        print("1. Prepare folders")
        self._unix_runtime(self._prepare_folders)
        print("2. Import primers, probes and contigs")
        self._unix_runtime(self._import_sequences)
        print("3. Create Bowtie index with contigs")
        self._unix_runtime(self._create_bowtie_index)
        print("4. Run Bowtie using primers and probes as query")
        self._unix_runtime(self._bowtie_call)
        print("5. Specificity calls")
        self._check_prebuilt_indexes()
        self._specificity_calls()
        print("5. Convert Bowtie output to BLAST outout")
        self._unix_runtime(self._multiprocess_convert_bowtie_to_blast)
        print("6. Split BLAST tab output file by query")
        self._unix_runtime(self._split_output)
        print("7. Parse BLAST tab output files for meaningful hits using multiprocessing, then merge the results")
        self._unix_runtime(self._multiprocess_split_files)
        print("8. Compare primers and find matching probes")
        self._unix_runtime(self._multiprocess_hits)
        print("\n9. Generate complete CSV output")
        self._unix_runtime(self._save_results)
        print("10. Process results for wobbled primers")
        self._unix_runtime(self._parse_results_to_wobbled)
        print("11. Generate output with wobbled oligos")
        self._unix_runtime(self._save_wobbled_results)
        print("12. Tidying up")
        self._unix_runtime(self._clean_up_folders)

    def _run_string_analysis(self):
        print("1. Prepare folders")
        self._unix_runtime(self._prepare_folders)
        print("2. Import primers, probes and contigs")
        self._unix_runtime(self._import_sequences)
        print("3. Add all oligos (including mismatches) to Aho-Corasick Automaton")
        self._unix_runtime(self._add_oligos_to_automaton)
        print("4. Run string search using the Aho-Corasick algorithm")
        self._unix_runtime(self._multiprocess_string_search)
        print("5. Split BLAST tab output file by query")
        self._unix_runtime(self._split_output)
        print("6. Parse BLAST tab output files for meaningful hits using multiprocessing, then merge the results")
        self._unix_runtime(self._multiprocess_split_files)
        print("7. Compare primers and find matching probes")
        self._unix_runtime(self._multiprocess_hits)
        print("8. Generate complete CSV output")
        self._unix_runtime(self._save_results)
        print("9. Process results for wobbled primers")
        self._unix_runtime(self._parse_results_to_wobbled)
        print("10. Generate output with wobbled oligos")
        self._unix_runtime(self._save_wobbled_results)
        print("11. Tidying up")
        self._unix_runtime(self._clean_up_folders)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Align primers and probes to contigs.")
    parser.add_argument("--primermismatches", action="store", dest="primer_mismatches", default=0,
                        help="Number of total mismatches allowed in a primer (default: 0).")
    parser.add_argument("--probemismatches", action="store", dest="probe_mismatches", default=0,
                        help="Number of total mismatches allowed in a probe (default: 0).")
    parser.add_argument("--runid", action="store", dest="run_id", default=999,
                        help="Run ID given by server (default: 999).")
    parser.add_argument("--maxproductsize", action="store", dest="product_size", default=1000,
                        help="Maximum product size allowed (default: 1000)")
    parser.add_argument("--method", action="store", dest="method", default="blast",
                        help="Method. Options: blast, bowtie, aho-corasick (default: blast).")
    parser.add_argument("--crosscheck", action="store", dest="cross_check", default="False",
                        help="Cross-check oligo packages to matches oligos across packages (default: False).")
    parser.add_argument("--probesonly", action="store", dest="probes_only", default="False",
                        help="Only check the probes, ignore primers (default: False).")
    parser.add_argument("--dbs", action="store", dest="dbs", default="",
                        help="Use prebuild Bowtie or BLAST dbs for a specificity check.")
    options = parser.parse_args()

    primer_mismatches = int(options.primer_mismatches)
    probe_mismatches = int(options.probe_mismatches)
    max_product_size = int(options.product_size)
    cross_check = str(options.cross_check)
    probes_only = str(options.probes_only)
    run_id = int(options.run_id)
    method = str(options.method)
    dbs = str(options.dbs)

    pe = PRIMevalRun(run_id, primer_mismatches, probe_mismatches, max_product_size, cross_check, probes_only, method, False, 50, 1.5, 1.75, 300, 25, 10, 0, 0, dbs)

pe.run_analysis()
