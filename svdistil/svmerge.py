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
from intervaltree import Interval, IntervalTree
import networkx as nx


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
        tree[start:end] = val

    def lookup(self, chrom, pos):
        if chrom in self.chroms:
            return self.chroms[chrom][pos]
        else:
            return set()


# mapping from unique integer (count) to variant record
class Variants(object):
    def __init__(self):
        self.variants = {}
        self.count = 0

    def add(self, variant):
        self.variants[self.count] = variant
        self.count += 1


def read_variants(variants, filename):
    with open(filename) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            variants.add(row)


def read_tsv_files(options):
    '''
    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    variants = Variants()
    for tsv_filename in options.tsv_files:
        logging.info("Processing TSV file from %s", tsv_filename)
        read_variants(variants, tsv_filename)
    return variants


def bnd_intervals(window, variants):
    window_left = window // 2
    window_right = window - window_left
    intervals_low = BndIntervals()
    intervals_high = BndIntervals()
    for variant_id, variant_info in variants.items():
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
    return intervals_low, intervals_high


def get_intersections(variants, intervals_low, intervals_high):
    overlaps = nx.Graph() 
    for variant_id, variant_info in variants.items():
        # make sure all variants are recorded in the graph
        overlaps.add_node(variant_id)
        chrom1 = variant_info['chr1']
        pos1 = int(variant_info['pos1'])
        intersections1 = intervals_low.lookup(chrom1, pos1)
        chrom2 = variant_info['chr2']
        pos2 = int(variant_info['pos2'])
        intersections2 = intervals_high.lookup(chrom2, pos2)
        overlaps_both_ends = intersections1.intersection(intersections2)
        for overlapping_variant in overlaps_both_ends:
            # don't add self edges
            overlapping_variant_id = overlapping_variant.data
            if variant_id != overlapping_variant_id:
                overlapping_variant_info = variants[overlapping_variant_id]
                if (variant_info['sense1'] == overlapping_variant_info['sense1']) and \
                   (variant_info['sense2'] == overlapping_variant_info['sense2']):
                    overlaps.add_edge(variant_id, overlapping_variant_id)
    return overlaps


def list_median(items):
    mid_pos = len(items) // 2
    return sorted(items)[mid_pos] 

def average(items):
    return (sum(items) / len(items))

def merge_overlaps(variants, overlaps):
    writer = csv.writer(sys.stdout, delimiter="\t")
    header = ["chr1", "pos1", "chr2", "pos2", "sense1", "sense2", "qual", "num samples", "samples"]
    writer.writerow(header)
    for component in nx.connected_components(overlaps):
        if len(component) > 0:
            variant_infos = [variants[id] for id in component]
            first_info = variant_infos[0]
            chrom1 = first_info['chr1']
            chrom2 = first_info['chr2']
            sense1 = first_info['sense1']
            sense2 = first_info['sense2']
            pos1 = list_median([info['pos1'] for info in variant_infos])
            pos2 = list_median([info['pos2'] for info in variant_infos])
            qual = average([float(info['qual']) for info in variant_infos]) 
            samples = ";".join([info['sample'] for info in variant_infos])
            writer.writerow([chrom1, pos1, chrom2, pos2, sense1, sense2, qual, len(component), samples]) 


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
    if log_filename is not None:
        logging.basicConfig(filename=log_filename,
                            level=logging.DEBUG,
                            filemode='w',
                            format='%(asctime)s %(levelname)s - %(message)s',
                            datefmt='%m-%d-%Y %H:%M:%S')
        logging.info('program started')
        logging.info('command line: %s', ' '.join(sys.argv))


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    variants = read_tsv_files(options)
    intervals_low, intervals_high = bnd_intervals(options.window, variants.variants)
    overlaps = get_intersections(variants.variants, intervals_low, intervals_high)
    merge_overlaps(variants.variants, overlaps)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
