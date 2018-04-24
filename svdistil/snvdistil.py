'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Convert SNVs and INDELS variants in VCF files into BED format.
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
PROGRAM_NAME = "snvdistil"


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


# XXX check this
def get_samples_with_variant(samples, genotypes):
    return [(sample, gt) for (sample, gt) in zip(samples, genotypes) if gt != 0]


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


def parse_gt_bases(ref, alts, base):
    seq_map = { seq: str(pos) for (pos, seq) in enumerate(alts, 1) }
    seq_map[ref] = "0"
    seq_map["."] = "." 
    alleles = base.split('/')
    return "/".join([seq_map[a] for a in alleles])

'''
Consequence|IMPACT|Codons|Amino_acids|Gene|SYMBOL|Feature|EXON|PolyPhen|SIFT|Protein_position|BIOTYPE|HGVSc|HGVSp|cDNA_position|CDS_position|HGVSc|HGVSp|cDNA_position|CDS_position|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MaxEntScan_alt|MaxEntScan_diff|MaxEntScan_ref|PICK
'''
def parse_vep(consequence_info):
    fields = consequence_info.split("|")
    consequence = fields[0]
    #impact = fields[1]
    gene = fields[5]
    gnomad_af = fields[20]
    return [consequence, gene, gnomad_af]


def process_variants(writer, qual_thresh, filter_pass, samples, vcf):
    results = set()
    for var in vcf:
        qual = var.QUAL
        samples_with_variant = get_samples_with_variant(samples, var.gt_types)
        if keep_variant(qual_thresh, filter_pass, var) and len(samples_with_variant) > 0:
            if qual is not None:
                qual_str = "{:.2f}".format(qual)
            else:
                qual_str = "."
            alt_str = ";".join(var.ALT)
            consequences = var.INFO['CSQ'].split(',')
            annotations = parse_vep(consequences[0]) 
            genotypes = [parse_gt_bases(var.REF, var.ALT, base) for base in var.gt_bases]
            new_row = tuple([var.CHROM, var.POS, var.REF, alt_str, var.var_type] + genotypes)
            writer.writerow(new_row)


def process_files(options):
    writer = csv.writer(sys.stdout, delimiter="\t")
    header = ["chr", "pos", "qual", "sample"]
    writer.writerow(header)
    for vcf_filename in options.vcf_files:
        logging.info("Processing VCF file from %s", vcf_filename)
        vcf = VCF(vcf_filename)
        samples = vcf.samples
        process_variants(writer, options.qual, options.ispass, samples, vcf)



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
