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
    parser = argparse.ArgumentParser(description='Take a file with a command \
            in it with a placeholder and iterate over a range to create a \
            file with the same command over and over with only the iterator \
            as the difference. Will pad the left side of the number with \
            zeroes. \
            output: commandList.txt \
            requires a filename, a range in the format high-low, and the \
            placeholder.')

    parser.add_argument('-f', help='Filename (first line only will be used)', metavar='Filename')
    parser.add_argument('-r', help='iteratorRange', metavar='numRange')
    parser.add_argument('-p', help='placeholder to substiute', metavar='Placeholder')
    parser.add_argument('-pad', help='number of places to pad with zeroes', metavar='padding', default=0, type=int)
    parser.add_argument('-d', help='debug (print a lot of shite, and all reports)', action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = vars(parser.parse_args())

    if not args['f']:
        argError("Please specify a filename with -f.", parser)
    if not args['r']:
        argError("Please specify an iterator range with -r.", parser)
    if not args['p']:
        argError("Please specify what to replace with -p.", parser)

    return args

def argError(message, parser):
    parser.print_help()
    print("\n" + message)
    sys.exit(1)

    

def makeCommands(template, placeholder, numRange, prefix='', pad=0):
    # Takes a string, a number range, and a number to replace and 
    # returns a list of strings with that number replaced iteratively
    newCommands = []
    low, high = numRange.split('-', 1)

    #for i in range(int(low), (int(high)-int(low))):
    for i in range(int(low), (int(high) + 1)):
        newCommands.append(template.replace(placeholder, str(i).zfill(pad)))

    return newCommands

if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False

    with open(args['f'], 'r') as f:
        line = f.readline()
    commands = makeCommands(line, args['p'], args['r'], '', args['pad'])
    with open('./commandList.txt', 'w') as out:
            for line in commands:
                out.write(line)

                

