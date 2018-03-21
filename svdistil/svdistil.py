'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Convert DNA structural variants in VCF files into BED format.
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv
import re
from cyvcf2 import VCF


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_VCF_FILE_ERROR = 3
DEFAULT_VERBOSE = False
DEFAULT_MIN_QUAL_THRESHOLD = 0.0
PROGRAM_NAME = "svdistil"


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
    description = 'Convert DNA structural variants in VCF files into BED format'
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

'''
According to VCF 4.2 spec, section 5.4

There are 4 possible ways to create the ALT in a SVTYPE=BND. In each of the 4 cases,
the assertion is that s (the REF) is replaced with t, and then some piece starting at
position p is joined to t. The cases are:

s t[p[ piece extending to the right of p is joined after t
s t]p] reverse comp piece extending left of p is joined after t
s ]p]t piece extending to the left of p is joined before t
s [p[t reverse comp piece extending right of p is joined before t
'''

# Match the 
INFO_ALT_REGEX = re.compile(r"\w*(\[|\])(?P<chrom>\w+)\:(?P<pos>\d+)(\[|\]).*")

def parse_bnd(info_alt):
    if len(info_alt) == 1:
        first_alt = info_alt[0]
    else:
        exit_with_error("BND ALT field without exactly one entry: {}".format(info_alt),
            EXIT_VCF_FILE_ERROR) 
    match = INFO_ALT_REGEX.match(first_alt)
    if match is not None:
        return match.group('chrom'), match.group('pos')
    else:
        exit_with_error("Cannot parse coordinate from BND variant ALT field: {}".format(first_alt),
            EXIT_VCF_FILE_ERROR) 

# Normalise the coordinates of an SV 
class NormSV(object):
    def __init__(self, var):
        info = var.INFO
        sv_type = info.get("SVTYPE", ".")
        # VCF parser returns zero-based coordinate 
        pos1 = var.start + 1
        if sv_type == 'BND':
            self.chrom1 = var.CHROM
            self.pos1 = pos1 
            self.chrom2, self.pos2 = parse_bnd(var.ALT)
        else:
            self.chrom1 = self.chrom2 = var.CHROM
            self.pos1 = pos1 
            self.pos2 = var.end
        if self.chrom1 != self.chrom2:
            self.type = "ITX" 
        else:
            self.type = "BND"



# XXX check this
def get_samples_with_variant(samples, genotypes):
    return [sample for (sample, gt) in zip(samples, genotypes) if gt != 0]


def process_variants(writer, qual_thresh, samples, vcf):
    for var in vcf:
        qual = var.QUAL
        if qual >= qual_thresh:
            samples_with_variant = get_samples_with_variant(samples, var.gt_types)
            qual_str = "{:.2f}".format(qual)
            norm = NormSV(var)
            row = [norm.chrom1, norm.pos1, norm.chrom2, norm.pos2, qual_str, norm.type] + samples_with_variant 
            writer.writerow(row)


def process_files(options):
    '''
    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    writer = csv.writer(sys.stdout, delimiter="\t")
    for vcf_filename in options.vcf_files:
        logging.info("Processing VCF file from %s", vcf_filename)
        vcf = VCF(vcf_filename)
        samples = vcf.samples
        process_variants(writer, options.qual, samples, vcf)


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
