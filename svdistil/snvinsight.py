'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Annotate SNVs with insight classifications 
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
# padding on either side of gene coordinates
DEFAULT_PAD = 2000
PROGRAM_NAME = "snvinsight"


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
    parser.add_argument('--insight',
                        metavar='INSIGHT',
                        required=True,
                        type=str,
                        help='insight TSV file from ClinVar')
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


CLASS_RE = re.compile(r"(?P<significance>.+)\(Last reviewed: (?P<date>[^\)]+)\)")

def parse_classification(classification):
    match = CLASS_RE.match(classification)
    if match is not None:
        significance = match.group('significance')
        date = match.group('date')
        return (significance, date)
    else:
        exit_with_error("Cannot parse classification: {}".format(classification), EXIT_FILE_IO_ERROR)


VAR_RE = re.compile(r"(?P<nm>.+)\((?P<gene>[^\)]+)\):(?P<change>.+)")

def drop_parens(str):
    if str.startswith('('):
        return str[1:-1]
    else:
        return str

def parse_change(change):
    fields = change.split()
    if len(fields) == 2:
        cchange = fields[0]
        pchange = fields[1]
        return (cchange, drop_parens(pchange))
    else:
        return (change, '')

def read_insight(insight_filepath):
    genes = set() 
    classifications = {}
    with open(insight_filepath) as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
            this_variant = row['Name']
            match = VAR_RE.match(this_variant)
            if match is not None:
                nm = match.group('nm').strip()
                gene = match.group('gene').strip()
                change = match.group('change').strip()
                cchange, pchange = parse_change(change)
                #print((this_variant, nm, gene, cchange, pchange))
                genes.add(gene)
                this_classification, this_date = parse_classification(row['Clinical significance (Last reviewed)'])
                if gene not in classifications:
                    classifications[gene] = {}
                classifications[gene][cchange] = (this_classification, this_date)
            else:
                logging.info("Cannot parse variant: {}".format(this_variant))
    return genes, classifications


def annotate_variants(genes, classifications):
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    header = reader.fieldnames
    output_header = header + ["insight class", "insight date"]
    writer = csv.DictWriter(sys.stdout, output_header, delimiter="\t")
    writer.writeheader()
    for row in reader:
        this_gene = row['gene']
        this_variant = row['hgvsc']
        insight_class = ''
        insight_date = ''
        if this_gene and this_variant:
            this_change = this_variant.split(':')[1]
            if this_gene in genes:
                insight_class, insight_date = classifications[this_gene].get(this_change, ('', ''))
        row['insight class'] = insight_class
        row['insight date'] = insight_date
        writer.writerow(row)



def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    genes, classifications = read_insight(options.insight)
    annotate_variants(genes, classifications)

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
