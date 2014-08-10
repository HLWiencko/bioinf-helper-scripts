#!/usr/bin/python3

# imports
import sys
import argparse
import matplotlib.pyplot as plt
from random import expovariate

## for me
import inspect
def lineno():
    """Returns the current line number in our program."""
    return "line " + str(inspect.currentframe().f_back.f_lineno) + ": "

def parseArgs():
    parser = argparse.ArgumentParser(description='Take a FASTA file and \
            produce a plot with the distribution of sequence lengths.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    parser.add_argument('-plot', help='Show the plot', action='store_true')
    parser.add_argument('-best', help='Return the longest n contigs', metavar='int', type=int)
    parser.add_argument('-d', help='debug (print a lot of shite, and all reports)', action='store_true')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')

    #if len(sys.argv) == 1:
    #    parser.print_help()
    #    sys.exit(1)

    #return vars(args)
    return vars(parser.parse_args())

# this kind of doesn't work. Use histogramBiopython.py instead
def plotDict(data, show=True):
    # expects data in a dictionary
    # number -> total elements with that number

    # create a Figure object
    fig = plt.figure(figsize=(5, 4))
    # create an Axes object
    ax = fig.add_subplot(1,1,1) # one row, one column, first plot
    # plot the data
    ax.scatter(list(data.keys()), list(data.values()), color='red', marker=',')
    ax.set_title('Distribution')
    ax.set_xlabel('length')
    ax.set_ylabel('# sequences')

    # add lines for median and percentiles

    # show the plot
    if show:
        plt.show()

def findPercentile(data, percentile):
    # expects data in a dictionary (doesn't have to be sorted)
    # number -> total elements with that number

    if 100 > percentile < 0:
        print(lineno(), "Percentile should be between 0 and 100.")
        return false

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

    if not args['f']:
        x_data = [x for x in range(0, 100)]
        y_data = [int(expovariate(1/50)) for x in x_data]
        totals = dict(zip(x_data, y_data))
        if debug:
            print(lineno(), totals)
    else: 
        with open(args['f'], 'r') as f:
            lines = f.readlines()
        seqs = getSeqs(lines)
        totals = countSeqs(seqs)
        if debug:
            print(lineno(), totals)
        if args['best']:
            writeBest(seqs, args['best'])

    if args['plot']:
        plotDict(totals)
        




