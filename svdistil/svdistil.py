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
from functools import total_ordering
from cyvcf2 import VCF


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_VCF_FILE_ERROR = 3
DEFAULT_VERBOSE = False
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
                        help='minimum QUAL threshold, variants below this will be discarded')
    parser.add_argument('--ispass',
                        action='store_true',
                        help='only keep variants whose filter field is PASS')
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

# t[p[, bp1 is right of pos1, bp2 is left of pos2 
INFO_ALT_REGEX_1 = re.compile(r"(?P<replacement>\w+)\[(?P<chrom>\w+)\:(?P<pos>\d+)\[")
# t]p], bp1 is right of pos1, bp2 is right of pos2
INFO_ALT_REGEX_2 = re.compile(r"(?P<replacement>\w+)\](?P<chrom>\w+)\:(?P<pos>\d+)\]")
# ]p]t, bp1 is left of pos1, bp2 is right of pos2
INFO_ALT_REGEX_3 = re.compile(r"\](?P<chrom>\w+)\:(?P<pos>\d+)\](?P<replacement>\w+)")
# [p[t, bp1 is left of pos1, bp2 is left of pos2
INFO_ALT_REGEX_4 = re.compile(r"\[(?P<chrom>\w+)\:(?P<pos>\d+)\[(?P<replacement>\w+)")

# XXX handle multiple ALTs
def parse_bnd(info_alt):
    if len(info_alt) == 1:
        first_alt = info_alt[0]
    else:
        exit_with_error("BND ALT field without exactly one entry: {}".format(info_alt),
            EXIT_VCF_FILE_ERROR) 
    match1 = INFO_ALT_REGEX_1.match(first_alt)
    match2 = INFO_ALT_REGEX_2.match(first_alt)
    match3 = INFO_ALT_REGEX_3.match(first_alt)
    match4 = INFO_ALT_REGEX_4.match(first_alt)
    if match1 is not None:
        return match1.group('chrom'), int(match1.group('pos')), match1.group('replacement'), "R", "L"
    elif match2 is not None:
        return match2.group('chrom'), int(match2.group('pos')), match2.group('replacement'), "R", "R"
    elif match3 is not None:
        return match3.group('chrom'), int(match3.group('pos')), match3.group('replacement'), "L", "R"
    elif match4 is not None:
        return match4.group('chrom'), int(match4.group('pos')), match4.group('replacement'), "L", "L"
    else:
        exit_with_error("Cannot parse coordinate from BND variant ALT field: {}".format(first_alt),
            EXIT_VCF_FILE_ERROR) 

@total_ordering
class Chrom(object):
    def __init__(self, name):
        if name.startswith('chr'):
            self.name = name[3:]
        else:
            self.name = name
        if len(self.name) == 0:
            exit_with_error("Empty chromosome name", EXIT_VCF_FILE_ERROR)
    def __eq__(self, other):
        return self.name == other.name
    def __lt__(self, other):
        return self.name < other.name
    def __str__(self):
        return "chr" + self.name 
    def __hash__(self):
        return hash(self.name)


# breakside indicates the side L|R of the position in which the breakpoint occurs,
# when using the + strand orientation.
#
#       DNA
#       --------X   breakpoint occurs on the right of X
#
#                   DNA
#               X-------- breakpoint occurs on the left of X
#
# we drop the chr from the start of chrom names

@total_ordering
class BreakEnd(object):
    def __init__(self, chrom, pos, breakside):
        self.chrom = Chrom(chrom)
        self.pos = pos
        self.breakside = breakside
    def __eq__(self, other):
        return (self.chrom, self.pos, self.breakside) == \
               (other.chrom, other.pos, other.breakside) 
    def __lt__(self, other):
        return (self.chrom, self.pos, self.breakside) < \
               (other.chrom, other.pos, other.breakside) 

# Normalise the coordinates of an SV 
#
# bnd_low: is the "lowest" of the 2 breakends.
#    - if they are on the same chrom, then it is the one with the least position
#    - if they are on different chroms, then it is the one with the lowest chrom
# bnd_high: is correspondingly the "highest" of the 2 breakends. 
#
# replacement: is an inserted sequence, if present
class NormSV(object):
    def __init__(self, var):
        info = var.INFO
        sv_type = info.get("SVTYPE", ".")
        # VCF parser returns zero-based coordinate 
        pos1 = var.start + 1
        if sv_type == 'BND':
            bnd2_chrom, bnd2_pos, replacement, breakside1, breakside2 = parse_bnd(var.ALT)
            bnd1 = BreakEnd(var.CHROM, pos1, breakside1)
            bnd2 = BreakEnd(bnd2_chrom, bnd2_pos, breakside2)
            self.bnd_low = min(bnd1, bnd2)
            self.bnd_high = max(bnd1, bnd2)
            self.replacement = replacement
        # the following are all on the same chrom
        elif sv_type in ['DEL', 'INV', 'DUP', 'INS'] :
            bnd1 = BreakEnd(var.CHROM, var.start, '.')
            end = int(info["END"])
            bnd2 = BreakEnd(var.CHROM, end, '.')
            self.bnd_low = min(bnd1, bnd2)
            self.bnd_high = max(bnd1, bnd2)
            self.replacement = ''
        else:
            exit_with_error("Unsupported SVTYPE: {}".format(sv_type), EXIT_VCF_FILE_ERROR)


# XXX check this
def get_samples_with_variant(samples, genotypes):
    return [sample for (sample, gt) in zip(samples, genotypes) if gt != 0]


def keep_variant(qual_thresh, filter_pass, var):
    passes_qual_thresh = False
    passes_filter = False

    if qual_thresh is None:
        passes_qual_thresh = True
    elif var.QUAL is None:
        passes_qual_thresh = qual_thresh == 0 
    else:
        passes_qual_thresh = var.QUAL >= qual_thresh

    if filter_pass:
        passes_filter = var.FILTER is None
    else:
        passes_filter = True

    return passes_qual_thresh and passes_filter


def process_variants(qual_thresh, filter_pass, samples, vcf):
    results = set()
    for var in vcf:
        qual = var.QUAL
        samples_with_variant = get_samples_with_variant(samples, var.gt_types)
        if keep_variant(qual_thresh, filter_pass, var) and len(samples_with_variant) > 0:
            if qual is not None:
                qual_str = "{:.2f}".format(qual)
            else:
                qual_str = "."
            norm = NormSV(var)
            bnd_low = norm.bnd_low
            bnd_high = norm.bnd_high
            samples_str = ";".join(samples_with_variant)
            new_row = (bnd_low.chrom, bnd_low.pos, bnd_high.chrom, bnd_high.pos,
                       bnd_low.breakside, bnd_high.breakside, len(norm.replacement), qual_str, samples_str)
            results.add(new_row)
    return sorted(set(results))


def process_files(options):
    '''
    Arguments:
       options: the command line options of the program
    Result:
       None
    '''
    writer = csv.writer(sys.stdout, delimiter="\t")
    header = ["chr1", "pos1", "chr2", "pos2", "sense1", "sense2", "insertlen", "qual", "sample"]
    writer.writerow(header)
    for vcf_filename in options.vcf_files:
        logging.info("Processing VCF file from %s", vcf_filename)
        vcf = VCF(vcf_filename)
        samples = vcf.samples
        results = process_variants(options.qual, options.ispass, samples, vcf)
        for row in results:
            writer.writerow(row)



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
