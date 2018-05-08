'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Test number of variants per gene across samples
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
import numpy as np
from scipy import stats
import plotly
import plotly.plotly as py
import plotly.graph_objs as go


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_TSV_FILE_ERROR = 3
DEFAULT_VERBOSE = False
# padding on either side of gene coordinates
DEFAULT_PAD = 2000
PROGRAM_NAME = "snvgene"


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


samples = ["0131313009", "0131326001", "0151032037", "0151052001",  "0151078001", "0151081002", "0151095001", "0350124001", "0350355001", "0636307001", "0636307011", "0656045001", "0656050001", "0757017010", "0757045001", "0757045010", "9930087001", "C4055110001", "E60176"]

def process_variants():
    counts = {}
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    for row in reader:
        gene = row['gene']
        num_carriers = int(row['num carriers'])
        if gene and num_carriers == 1:
            if gene not in counts:
                counts[gene] = {}
            for sample in samples:
                if "1" in row[sample]:
                    if sample not in counts[gene]:
                        counts[gene][sample] = 0
                    counts[gene][sample] += 1
                    break
    sample_headers = [[sample + " count", sample + " z score", sample + " p value"] for sample in samples]
    flat_sample_headers = [x for sublist in sample_headers for x in sublist]
    headers = "\t".join(["gene"] + flat_sample_headers)
    print(headers)
    for gene in sorted(counts):
        sample_counts = []
        for sample in samples:
            if sample in counts[gene]:
                sample_counts.append(counts[gene][sample]) 
            else:
                sample_counts.append(0)
        zscores = stats.zscore(np.array(sample_counts))
        pvalues = stats.norm.sf(abs(zscores))*2
        scores = [[c, z, p] for (c, z, p) in zip(sample_counts, zscores, pvalues)]
        row_data = [gene] + [str(x) for sublist in scores for x in sublist] 
        print("\t".join(row_data))
        if gene == 'MSH2':
            plot(sample_counts)


def plot(counts):
    counts = np.array(counts)
    trace1 = go.Histogram(
        x=counts,
        histnorm='count',
        name='MSH2 sample-unique variants',
        xbins=dict(
            start=counts.min(),
            end=counts.max(),
            size=5
        ),
    )

    data = [trace1]

    layout = go.Layout(
        title='MSH2 sample-unique variants',
        xaxis=dict(
            title='Number of unique variants'
        ),
        yaxis=dict(
            title='Number of samples'
        ),
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename='MSH2.histogram.html')
    '''
    histogram_data = [{
        'x' : counts,
        'type': 'histogram'
    }]
    histogram_layout = {'title': 'histogram of unique variants in MSH2',
        'xaxis': {'title': 'number of unique variants'}, 'yaxis': {'title': 'count'}}
    histogram_fig = {'data': histogram_data, 'layout': histogram_layout}
    py.offline.plot(histogram_fig, filename="MSH2.histogram")
    '''


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    process_variants()

# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
