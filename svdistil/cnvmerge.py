'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Merge distilled SVs
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
from quicksect import IntervalTree
import networkx as nx
from pathlib import Path
from copy import copy
import re


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_TSV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_OVERLAP = 0.75 
PROGRAM_NAME = "cnvmerge"


try:
    PROGRAM_VERSION = pkg_resources.require(PROGRAM_NAME)[0].version
except pkg_resources.DistributionNotFound:
    PROGRAM_VERSION = "undefined_version"


def exit_with_error(message, exit_status):
    '''Print an error message to stderr, prefixed by the program name and 'ERROR'.
    Then exit program with supplied exit status.

    Arguments:
        message: an error message as a string.
        exit_status: a positive integer representing the exit status of the
            program.
    '''
    logging.error(message)
    print("{} ERROR: {}, exiting".format(PROGRAM_NAME, message), file=sys.stderr)
    sys.exit(exit_status)


def parse_args():
    '''Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    '''
    description = 'Merge distilled SVs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--overlap',
                        metavar='PERCENTAGE',
                        default=DEFAULT_OVERLAP,
                        type=float,
                        help='percentage overlap for CNV equality (default {})'.format(DEFAULT_OVERLAP))
    parser.add_argument('tsv_files',
                        nargs='*',
                        metavar='TSV_FILE',
                        type=str,
                        help='Input TSV files')
    return parser.parse_args()


class CNVIntervals(object):
    def __init__(self):
        self.chroms = {}

    def add(self, chrom, start, end, val):
        if chrom not in self.chroms:
            self.chroms[chrom] = IntervalTree()
        tree = self.chroms[chrom]
        tree.add(start, end, val)

    def lookup(self, chrom, start, end):
        if chrom in self.chroms:
            return self.chroms[chrom].search(start, end)
        else:
            return [] 


# mapping from unique integer (count) to variant record
class Variants(object):
    def __init__(self):
        self.variants = {}
        self.count = 0

    def add(self, variant):
        self.variants[self.count] = variant
        self.count += 1

def get_sample_name(filepath):
    fields = filepath.split('.')
    if len(fields) > 0:
        sample = fields[0]
    else:
        sample = filepath
    return sample

def read_tsv_files(options):
    sample_ids = set()
    variants = Variants()
    for tsv_filename in options.tsv_files:
        logging.info("Processing TSV file from %s...", tsv_filename)
        sample = get_sample_name(tsv_filename)
        sample_ids.add(sample)
        with open(tsv_filename) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                row['sample'] = sample
                variants.add(row)
        logging.info("Processing TSV file from %s: done", tsv_filename)
    return sample_ids, variants


def cnv_intervals(variants):
    logging.info("Computing %i CNV intervals", len(variants))
    intervals = CNVIntervals()
    for idx, (variant_id, variant_info) in enumerate(variants.items()):
        chrom = variant_info['chr']
        start = int(float(variant_info['start']))
        end = int(float(variant_info['end']))
        intervals.add(chrom, start, end, variant_id)
        if (idx + 1) % 100000 == 0:
            logging.info('Computing %i CNV intervals: %i done', len(variants), idx + 1)
    logging.info("Computing %i CNV intervals, done", len(variants))
    return intervals


def is_overlap(start1, end1, start2, end2, min_overlap):
    overlap_start = max(start1, start2)
    overlap_end = min(end1, end2)
    if overlap_start < overlap_end:
        overlap_size = float((overlap_end - overlap_start) + 1)
        cnv1_size = (end1 - start1) + 1
        cnv2_size = (end2 - start2) + 1
        cnv1_overlap = overlap_size / cnv1_size
        cnv2_overlap = overlap_size / cnv2_size
        return cnv1_overlap >= min_overlap and cnv2_overlap >= min_overlap
    return False

def get_intersections(overlap, variants, intervals):
    logging.info("Computing %i CNV intersections...", len(variants))
    overlaps = nx.Graph() 
    for idx, (variant_id, variant_info) in enumerate(variants.items()):
        # make sure all variants are recorded in the graph
        overlaps.add_node(variant_id)
        chrom = variant_info['chr']
        start = int(float(variant_info['start']))
        end = int(float(variant_info['end']))
        this_state = variant_info['state'] 
        intersections = { i for i in intervals.lookup(chrom, start, end) }
        for other_variant in intersections:
            other_variant_id = other_variant.data
            other_variant_info = variants[other_variant_id]
            # don't add self edges
            if variant_id != other_variant_id and \
               is_overlap(start, end, other_variant.start, other_variant.end, overlap) and \
               this_state == other_variant_info['state']:
                overlaps.add_edge(variant_id, other_variant_id)
        if (idx + 1) % 100000 == 0:
          logging.info("Computing %i variant intersections: %i done", len(variants), idx + 1)
    logging.info("Computing %i variant intersections: done", len(variants))
    return overlaps


def list_median(items):
    mid_pos = len(items) // 2
    return sorted(items)[mid_pos] 

def average(items):
    return (sum(items) / len(items))

def build_evidence(variants, samples):
    # evidence: mapping, sample -> set(caller)
    num_positive_samples = 0
    evidence = set()
    for var in variants:
        this_sample = var['sample']
        evidence.add(this_sample)
    results = []
    for sample in samples:
        if sample in evidence:
            num_positive_samples += 1
            results.append(1)
        else:
            results.append(0)
    return num_positive_samples, results


def merge_overlaps(sample_ids, variants, overlaps):
    logging.info("Merging overlapping variants...")
    writer = csv.writer(sys.stdout, delimiter="\t")
    sorted_samples = sorted(sample_ids)
    header = ["chr", "start", "end", "state", "median", "num pos samples"] + sorted_samples 
    writer.writerow(header)
    for component in nx.connected_components(overlaps):
        if len(component) > 0:
            variant_infos = [variants[id] for id in component]
            first_info = variant_infos[0]
            chrom  = first_info['chr']
            state = first_info['state']
            start = min([int(float(info['start'])) for info in variant_infos])
            end = max([int(float(info['end'])) for info in variant_infos])
            avg_median = average([float(info['median']) for info in variant_infos])
            num_positive_samples, evidence = build_evidence(variant_infos, sorted_samples)
            writer.writerow([chrom, start, end, state, avg_median, num_positive_samples] + evidence) 
    logging.info("Merging overlapping variants: done")


def init_logging(log_filename):
    '''If the log_filename is defined, then
    initialise the logging facility, and write log statement
    indicating the program has started, and also write out the
    command line from sys.argv

    Arguments:
        log_filename: either None, if logging is not required, or the
            string name of the log file to write to
    Result:
        None
    '''
    if log_filename is None:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
    else:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
    logging.info('computation started')
    logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    sample_ids, variants = read_tsv_files(options)
    intervals = cnv_intervals(variants.variants)
    overlaps = get_intersections(options.overlap, variants.variants, intervals)
    merge_overlaps(sample_ids, variants.variants, overlaps)
    logging.info("computation ended")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
