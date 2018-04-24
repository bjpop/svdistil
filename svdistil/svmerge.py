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


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_TSV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_BND_WINDOW = 50
PROGRAM_NAME = "svmerge"


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
    parser.add_argument('--window',
                        metavar='WINDOW',
                        default=DEFAULT_BND_WINDOW,
                        type=int,
                        help='window size centered on breakend for determining equality (default {})'.format(DEFAULT_BND_WINDOW))
    parser.add_argument('tsv_files',
                        nargs='*',
                        metavar='VCF_FILE',
                        type=str,
                        help='Input VCF files')
    return parser.parse_args()


class BndIntervals(object):
    def __init__(self):
        self.chroms = {}

    def add(self, chrom, start, end, val):
        if chrom not in self.chroms:
            self.chroms[chrom] = IntervalTree()
        tree = self.chroms[chrom]
        tree.add(start, end, val)

    def lookup(self, chrom, pos):
        if chrom in self.chroms:
            return self.chroms[chrom].search(pos, pos)
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


def get_caller_name(filepath):
    path = Path(filepath)
    filename = path.name
    fields = filename.split('.')
    if len(fields) > 0:
        return fields[0]
    else:
        return filename

def read_tsv_files(options):
    sample_ids = set()
    callers = set()
    variants = Variants()
    for tsv_filename in options.tsv_files:
        logging.info("Processing TSV file from %s...", tsv_filename)
        variant_caller = get_caller_name(tsv_filename)
        callers.add(variant_caller)
        with open(tsv_filename) as file:
            reader = csv.DictReader(file, delimiter="\t")
            for row in reader:
                samples = row['sample'].split(';')
                for sample in samples:
                    new_row = copy(row)
                    new_row['caller'] = variant_caller
                    new_row['sample'] = sample
                    variants.add(new_row)
                    sample_ids.add(sample)
        logging.info("Processing TSV file from %s: done", tsv_filename)
    return callers, sample_ids, variants


def bnd_intervals(window, variants):
    logging.info("Computing %i break end intervals with window size: %i...", len(variants), window)
    window_left = window // 2
    window_right = window - window_left
    intervals_low = BndIntervals()
    intervals_high = BndIntervals()
    for idx, (variant_id, variant_info) in enumerate(variants.items()):
        chrom1 = variant_info['chr1']
        pos1 = int(variant_info['pos1'])
        # make sure interval pos is not negative
        interval1_start = max(0, pos1 - window_left)
        interval1_end = pos1 + window_right
        intervals_low.add(chrom1, interval1_start, interval1_end, variant_id)
        chrom2 = variant_info['chr2']
        pos2 = int(variant_info['pos2'])
        # make sure interval pos is not negative
        interval2_start = max(0, pos2 - window_left)
        interval2_end = pos2 + window_right
        intervals_high.add(chrom2, interval2_start, interval2_end, variant_id)
        if (idx + 1) % 100000 == 0:
            logging.info('Computing %i break end intervals: %i done', len(variants), idx + 1)
    logging.info("Computing %i break end intervals with window size: %i: done", len(variants), window)
    return intervals_low, intervals_high


def get_intersections(variants, intervals_low, intervals_high):
    logging.info("Computing %i variant intersections...", len(variants))
    overlaps = nx.Graph() 
    for idx, (variant_id, variant_info) in enumerate(variants.items()):
        # make sure all variants are recorded in the graph
        overlaps.add_node(variant_id)
        chrom1 = variant_info['chr1']
        pos1 = int(variant_info['pos1'])
        intersections1 = { i.data for i in intervals_low.lookup(chrom1, pos1) }
        chrom2 = variant_info['chr2']
        pos2 = int(variant_info['pos2'])
        intersections2 = { i.data for i in intervals_high.lookup(chrom2, pos2) }
        overlaps_both_ends = intersections1.intersection(intersections2)
        for overlapping_variant_id in overlaps_both_ends:
            # don't add self edges
            if variant_id != overlapping_variant_id:
                overlapping_variant_info = variants[overlapping_variant_id]
                overlaps.add_edge(variant_id, overlapping_variant_id)
        if (idx + 1) % 100000 == 0:
          logging.info("Computing %i variant intersections: %i done", len(variants), idx + 1)
    logging.info("Computing %i variant intersections: done", len(variants))
    return overlaps


def list_median(items):
    mid_pos = len(items) // 2
    return sorted(items)[mid_pos] 

def average(items):
    return (sum(items) / len(items))

def build_evidence(variants, callers, samples):
    # evidence: mapping, sample -> set(caller)
    num_positive_samples = 0
    num_positive_calls = 0
    evidence = {}
    for var in variants:
        this_sample = var['sample']
        this_caller = var['caller']
        if this_sample not in evidence:
            evidence[this_sample] = set()
        evidence[this_sample].add(this_caller)
    results = []
    for sample in samples:
        if sample in evidence:
            num_positive_samples += 1
            positive_callers = evidence[sample]
            num_calls = len(positive_callers)
            num_positive_calls += num_calls
            caller_hits = []
            for caller in callers:
                if caller in positive_callers:
                    caller_hits.append(1)
                else:
                    caller_hits.append(0)
            results.extend([num_calls] + caller_hits)
        else:
            results.extend([0] + [0 for _ in callers])
    return num_positive_samples, num_positive_calls, results


def merge_overlaps(callers, sample_ids, variants, overlaps):
    logging.info("Merging overlapping variants...")
    writer = csv.writer(sys.stdout, delimiter="\t")
    sorted_samples = sorted(sample_ids)
    sorted_callers = sorted(callers)
    sample_headers = [s + " " + f for s in sorted_samples for f in ["count"] + sorted_callers]
    header = ["chr1", "pos1", "chr2", "pos2", "num samples", "avg pos calls"] + sample_headers 
    writer.writerow(header)
    for component in nx.connected_components(overlaps):
        if len(component) > 0:
            variant_infos = [variants[id] for id in component]
            first_info = variant_infos[0]
            chrom1 = first_info['chr1']
            chrom2 = first_info['chr2']
            pos1 = list_median([info['pos1'] for info in variant_infos])
            pos2 = list_median([info['pos2'] for info in variant_infos])
            num_positive_samples, num_positive_calls, evidence = build_evidence(variant_infos, sorted_callers, sorted_samples)
            avg_positive_calls = float(num_positive_calls) / num_positive_samples
            writer.writerow([chrom1, pos1, chrom2, pos2, num_positive_samples, avg_positive_calls] + evidence) 
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
    callers, sample_ids, variants = read_tsv_files(options)
    intervals_low, intervals_high = bnd_intervals(options.window, variants.variants)
    overlaps = get_intersections(variants.variants, intervals_low, intervals_high)
    merge_overlaps(callers, sample_ids, variants.variants, overlaps)
    logging.info("computation ended")


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
