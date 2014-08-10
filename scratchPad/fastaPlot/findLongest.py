#!/usr/bin/python3

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
            find the longest n sequences, write to a separate file.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    parser.add_argument('-best', help='Return the longest n contigs', metavar='int', type=int)
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


    return args

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
    return seqs

def countSeqs(sequences):
    # takes a dictionary of sequences and returns a dictionary of number -> sequences with that length
    totals = {}
    for seq in sequences.values():
        length = len(seq)
        if length in totals:
            totals[length] += 1
        else:
            totals[length] = 1
    return totals

def writeBest(data, n, outPrefix=''):
    outfile = outPrefix + 'best' + str(n) + '.fasta'
    best = sorted(data, key=lambda x: len(data[x]), reverse=True)[0:n]
    with open(outfile, 'w') as out:
        for label in best:
            for line in makeFasta(label, data[label]):
                out.write(line + '\n')
    print("Best %i contigs written to %s" % (n, outfile))

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
    totals = countSeqs(seqs)
    if debug:
        print(lineno(), totals)
    if args['best']:
        writeBest(seqs, args['best'])




