#!/usr/bin/python

# imports
import sys
import argparse

from Bio import SeqIO

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "

def parseArgs():
    parser = argparse.ArgumentParser(description='Take a FASTA file and \
            find the longest n sequences, write to a separate file.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    parser.add_argument('-n', help='Split file into files of n sequences', metavar='int', type=int)
    parser.add_argument('-d', help='debug (print a lot of shite, and all reports)', action='store_true')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())

    if not args['f']:
        parser.print_help()
        print("\nPlease specify a filename with -f.")
        sys.exit(1)
    if not args['n']:
        parser.print_help()
        print("\nPlease specify how many sequences per file you want with -n .")
        sys.exit(1)

    return args

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

    

if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False

    record_iter = SeqIO.parse(open(args['f']),"fasta")
    for i, batch in enumerate(batch_iterator(record_iter, args['n'])) :
        filename = "group_%i.fasta" % (i+1)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, "fasta")
        handle.close()
        print("Wrote %i records to %s" % (count, filename))




