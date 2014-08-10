#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3

###!/usr/bin/python3

# imports
import sys
import argparse

from Bio import Entrez, SeqIO
Entrez.email = 'timiat@yahoo.com' # good manners start with an introduction

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "

# command line options
def parseArgs():
    parser = argparse.ArgumentParser(description='Take a list of gene names\
            and an assembly, provide coordinates for those genes in the \
            given assembly. Currently only supports salmonella SL1344 \
            but could easily be extended if I ever use this script again.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    # TODO want to make this store either a string or a boolean depending on
    # user input
    parser.add_argument('-b', help='Output in BED format', metavar='BED\
           format', default='out.bed')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')
    parser.add_argument('-d', help='debug (print a lot of shite, and all \
            reports)', action='store_true')
    parser.add_argument('-t', help='test mode with toy input', \
            action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = vars(parser.parse_args())

    if not args['f']:
        if not args['t']:
            argError("Please specify a filename with -f.", parser)

    return vars(parser.parse_args())

def argError(message, parser):
    parser.print_help()
    print("\n" + message)
    sys.exit(1)

#####

# take a list of names, return a list of Seq objects
def getSequences(names):
    sequences = []
    for name in names:
        print(name)
        handle = Entrez.efetch(db='nucleotide', \
                id="186972394", \
                rettpye='gb', \
                retmode='text')
        print(handle.read())

#         record = Entrez.read(handle)
#         print(record['IdList'])
#         record = seqio.read(handle, 'fasta')
#         handle.close()
#         print(record)

# this works
#     handle = Entrez.einfo()
#     print(handle.read())

    return sequences


##### 

if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False

    if args['t']:
        names = ['hns', 'stpA']
    else: 
        names = [line.rstrip() for line in open(args['f'])]
    print(names)
    getSequences(names)




