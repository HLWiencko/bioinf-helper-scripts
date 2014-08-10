#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3

###!/usr/bin/python3

# imports
import sys
import argparse
from collections import defaultdict

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

def reverseCompliment(inputString):
    translate = {'A':'T', 'T':'A','G':'C','C':'G'}
    rev = inputString[::-1]
    outputString = ''
    for char in rev:
        outputString = outputString + translate[char]
    return outputString

def getSeqs(lines):
    # take a fasta file as a list of lines and return a dictionary of name -> sequence
    label = ''
    sequence = ''
    seqs = {}
    for line in lines:
        if '>' in line:
            seqs[label] = sequence
            sequence = ''
            label = line.lstrip('> ').rstrip('\n')
        else: 
            sequence = sequence + line
    seqs[label] = sequence
    return seqs

def writeSeq(data, outPrefix=''):
    # expects a dictionary name -> sequence, writes a fasta file
    outfile = outPrefix + 'out.fasta'
    with open(outfile, 'w') as out:
        for label in data:
            for line in makeFasta(label, data[label]):
                out.write(line + '\n')
    print("Reverse compliment of %s written to %s" % (infile, outfile))

def makeFasta(label, sequence):
    # takes a name and a sequence and returns a list [>label, seq1, seq2, seq3]
    fasta = []
    label = '>' + label
    fasta.append(label)
    buf = ''
    for l in sequence:
        if len(buf) < 60:
            buf += l
        else:
            fasta.append(buf)
            buf = ''
    fasta.append(buf)
    return fasta
    

if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False

