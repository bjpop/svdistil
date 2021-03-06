'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Filter variants based on quality score threshold
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import pathlib
from cyvcf2 import VCF


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_VCF_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_MIN_QUAL_THRESHOLD = 0.0
PROGRAM_NAME = "svqualfilter"


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
    description = 'Filter variants based on quality score'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--qual',
                        metavar='MIN_QUAL_THRESHOLD',
                        type=float,
                        default=DEFAULT_MIN_QUAL_THRESHOLD,
                        help='minimum QUAL threshold, variants below this will be discarded')
    parser.add_argument('vcf_files',
                        nargs='*',
                        metavar='VCF_FILE',
                        type=str,
                        help='Input VCF files')
    return parser.parse_args()


def process_variants(outfile, qual_thresh, vcf):
    for var in vcf:
        qual = var.QUAL
        info = var.INFO
        assembly1 = info['AS']
        assembly2 = info['RAS']
        if qual >= qual_thresh and assembly1 > 0 and assembly2 > 0:
            outfile.write(str(var))


def process_files(options):
    '''
    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    for vcf_filename in options.vcf_files:
        logging.info("Processing VCF file from %s", vcf_filename)
        vcf = VCF(vcf_filename)
        outfilepath = pathlib.Path(vcf_filename).with_suffix(".filtered.vcf")
        with outfilepath.open(mode='w') as outfile:
            process_variants(outfile, options.qual, vcf)


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
    process_files(options)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
