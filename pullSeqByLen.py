#!/Library/Frameworks/Python.framework/Versions/3.4/bin/python3

###!/usr/bin/python3

# imports
import sys
import argparse

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "

def parseArgs():
    parser = argparse.ArgumentParser(description='Take a FASTA file and \
            find sequences over a given length, write to a separate file.\
            Optional: write to individual files.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    parser.add_argument('-n', help='Minimum sequence length (bp)', metavar='int', type=int)
    parser.add_argument('-s', help='Write to individual files', action='store_true')
    parser.add_argument('-d', help='debug (print a lot of shite, and all reports)', action='store_true')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())

    if not args['f']:
        argError("Please specify a filename with -f.", parser)
    if not args['n']:
        argError("Please specify a minimum sequence length with -n.", parser)

    return args

def argError(message, parser):
    parser.print_help()
    print("\n" + message)
    sys.exit(1)


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

def chooseSeqs(data, n, outPrefix=''):
    # expects a dictionary name -> sequence, returns a filtered dictionary
    filtered = {}
    for label in data:
        if len(data[label]) >= n:
            filtered[label] = data[label]
    return filtered


def writeSeqs(data, n, outPrefix=''):
    # expects a dictionary name -> sequence, writes a fasta file
    outfile = outPrefix + 'out.fasta'
    with open(outfile, 'w') as out:
        for label in data:
            for line in makeFasta(label, data[label]):
                out.write(line + '\n')
    print("Sequences over %i written to %s" % (n, outfile))
def writeSplitFiles(data, n, outPrefix=''):
    i = 1
    numSeq = len(data)
    for label in data:
        # pad with leading zeroes so the files list nicely
        prefix = outPrefix + 'seq' + str(i).zfill(len(str(numSeq))) + '-'
        writeSeqs({label: data[label]}, n, prefix)
        i += 1
    

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

    with open(args['f'], 'r') as f:
        lines = f.readlines()
    seqs = getSeqs(lines)
    chosenSeqs = chooseSeqs(seqs, args['n'])
    # This gives some weird output in stdout, but it does the job.
    if args['s']:
        writeSplitFiles(chosenSeqs, args['n'])
    else:
        writeSeqs(chosenSeqs, args['n'])




