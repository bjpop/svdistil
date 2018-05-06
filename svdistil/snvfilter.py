'''
Module      : Main
Description : The main entry point for the program.
Copyright   : (c) Bernie Pope, 2018
License     : BSD-2-Clause 
Maintainer  : bjpope@unimelb.edu.au
Portability : POSIX

Filter SNVs
'''

from argparse import ArgumentParser
import sys
import logging
import pkg_resources
import csv


EXIT_FILE_IO_ERROR = 1
EXIT_COMMAND_LINE_ERROR = 2
EXIT_TSV_FILE_ERROR = 3
EXIT_FILTER_ERROR = 4
DEFAULT_VERBOSE = False
# padding on either side of gene coordinates
PROGRAM_NAME = "snvfilter"


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
    description = 'Filter SNVs'
    parser = ArgumentParser(description=description)
    parser.add_argument('--version',
                        action='version',
                        version='%(prog)s ' + PROGRAM_VERSION)
    parser.add_argument('--log',
                        metavar='LOG_FILE',
                        type=str,
                        help='record program progress in LOG_FILE')
    parser.add_argument('--filter',
                        metavar='FILTER_PYTHON_FILE',
                        required=True,
                        type=str,
                        help='file containing python filter function')
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


def get_filter(filter_filepath):
    # default filter is the identity function,
    # will include all rows in the output
    # WARNING: this code execs arbitrary Python code. Do not use this in
    # an untrusted environment, such as a web application!
    locals = {'filter': None}
    try:
        with open(filter_filepath) as filter_file:
            contents = filter_file.read()
            exec(contents, None, locals)
    except Exception as e:
        print(e)
        exit(EXIT_FILTER_ERROR)
    filter = locals['filter']
    if filter is not None:
        return filter
    else:
        exit_with_error("Could not load filter function", EXIT_FILTER_ERROR)


def filter_variants(row_filter):
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    fieldnames = reader.fieldnames
    writer = csv.DictWriter(sys.stdout, fieldnames, delimiter="\t")
    writer.writeheader()
    for row in reader:
        for filtered_row in row_filter(row):
            writer.writerow(filtered_row)


def main():
    "Orchestrate the execution of the program"
    options = parse_args()
    init_logging(options.log)
    row_filter = get_filter(options.filter)
    filter_variants(row_filter)


# If this script is run from the command line then call the main function.
if __name__ == '__main__':
    main()
