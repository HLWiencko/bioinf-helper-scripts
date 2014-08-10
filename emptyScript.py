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
    parser = argparse.ArgumentParser(description='Description goes here.')

    parser.add_argument('-f', help='Filename', metavar='filename')
    #parser.add_argument('-p', help='Prefix for output', metavar='prefix')
    parser.add_argument('-d', help='debug (print a lot of shite, and all \
            reports)', action='store_true')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    #return vars(args)

    args = vars(parser.parse_args())

    if not args['f']:
        argError("Please specify a filename with -f.", parser)

    return vars(parser.parse_args())

def argError(message, parser):
    parser.print_help()
    print("\n" + message)
    sys.exit(1)



if __name__ == '__main__':
    args = parseArgs()
    if args['d']: 
        debug = True
    else:
        debug = False




