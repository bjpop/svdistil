[![travis](https://travis-ci.org/bjpop/svdistil.svg?branch=master)](https://travis-ci.org/bjpop/svdistil)

# Overview 

Convert DNA structural variants in VCF files into BED format.

# Licence

This program is released as open source software under the terms of [BSD-2-Clause License](https://raw.githubusercontent.com/bjpop/svdistil/master/LICENSE).

# Installing

Svdistil can be installed using `pip` in a variety of ways (`%` indicates the command line prompt):

1. Inside a virtual environment:
```
% python3 -m venv svdistil_dev
% source svdistil_dev/bin/activate
% pip install -U /path/to/svdistil
```
2. Into the global package database for all users:
```
% pip install -U /path/to/svdistil
```
3. Into the user package database (for the current user only):
```
% pip install -U --user /path/to/svdistil
```


# General behaviour


# Usage 


```
% svdistil -h
```


## Logging

If the ``--log FILE`` command line argument is specified, svdistil will output a log file containing information about program progress. The log file includes the command line used to execute the program, and a note indicating which files have been processes so far. Events in the log file are annotated with their date and time of occurrence. 

```
% svdistil --log 
# normal svdistil output appears here
# contents of log file displayed below
```
```
% cat bt.log
```

# Exit status values

Svdistil returns the following exit status values:

* 0: The program completed successfully.
* 1: File I/O error. This can occur if at least one of the input FASTA files cannot be opened for reading. This can occur because the file does not exist at the specified path, or svdistil does not have permission to read from the file. 
* 2: A command line error occurred. This can happen if the user specifies an incorrect command line argument. In this circumstance svdistil will also print a usage message to the standard error device (stderr).
* 3: Input FASTA file is invalid. This can occur if svdistil can read an input file but the file format is invalid. 


# Error handling

# Testing

## Unit tests

```
% cd svdistil/python/svdistil
% python -m unittest -v svdistil_test
```

# Bug reporting and feature requests

Please submit bug reports and feature requests to the issue tracker on GitHub:

[svdistil issue tracker](https://github.com/bjpop/svdistil/issues)
