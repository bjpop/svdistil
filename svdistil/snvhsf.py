'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Annotate SNVs with human splicing finder results 
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
from quicksect import IntervalTree
from pathlib import Path
from copy import copy
import re


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_TSV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
PROGRAM_NAME = "snvhsf"


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
    description = 'Annotate SNVs and indels with gene tiers'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--hsf',
                        metavar='HSF',
                        required=True,
                        type=str,
                        help='human splicing finder TSV file')
    return parser.parse_args()


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
        logging.info('computation started')
        logging.info('command line: %s', ' '.join(sys.argv))


def read_hsf(hsf_filepath):
    classifications = {}
    with open(hsf_filepath) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            effect = row['human splicing finder effect']
            algorithms = row['human splicing finder prediction algorithm']
            gene = row['gene']
            hgvsc = row['hgvsc']
            if effect:
                if gene not in classifications:
                    classifications[gene] = {}
                classifications[gene][hgvsc] = (effect, algorithms)
    return classifications 


def annotate_variants(classifications):
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    header = reader.fieldnames
    output_header = header + ["hsf effect", "hsf algorithms"]
    writer = csv.DictWriter(sys.stdout, output_header, delimiter="\t")
    writer.writeheader()
    for row in reader:
        this_gene = row['gene']
        hsf_effect = ''
        hsf_alg = ''
        if this_gene in classifications:
            gene_classifications = classifications[this_gene]
            this_variant = row['hgvsc']
            hgvsc_fields = this_variant.split(':')
            if len(hgvsc_fields) == 2:
                hgvsc_raw = hgvsc_fields[1]
                if hgvsc_raw in gene_classifications:
                    hsf_effect, hsf_alg = gene_classifications[hgvsc_raw]
        row['hsf effect'] = hsf_effect
        row['hsf algorithms'] = hsf_alg
        writer.writerow(row)



def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    classifications = read_hsf(options.hsf)
    annotate_variants(classifications)

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
