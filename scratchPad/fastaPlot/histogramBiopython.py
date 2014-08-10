#!/usr/bin/python3

# imports
import sys
import argparse

from Bio import SeqIO
import pylab

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "

def parseArgs():
    parser = argparse.ArgumentParser(description='Description goes here.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')
    parser.add_argument('-d', help='debug (print a lot of shite, and all reports)', action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #return vars(args)
    return vars(parser.parse_args())

def plotSizes(sizes):
    pylab.hist(sizes, bins=10)
    pylab.title("%i Contigs\nLengths %i to %i" \
    % (len(sizes),min(sizes),max(sizes)))
    pylab.xlabel("Contig length (bp)")
    pylab.ylabel("Count")
    pylab.show()


if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False

    sizes = [len(rec) for rec in SeqIO.parse(args['f'], "fasta")]
    plotSizes(sizes)



