#!/usr/bin/env python

import argparse
import sys

# -------------------------------------

def get_arg():
    """
    Get Arguments
    :rtype: object
    """
    # parse arguments

    prog_description = 'preprocess'
    parser = argparse.ArgumentParser(description=prog_description)
    parser.add_argument('-d', '--scripts', help='the git repository directory')
    parser.add_argument('-o', '--outputdir', default='aggregate_preprocess', help='the output directory')
    args = parser.parse_args()

    # need this to get local modules
    sys.path.append(args.scripts)
    global hp
    from helpers import helpers as hp

    # add key-value pairs to the args dict
    vars(args)['step'] = 'preprocess'
    vars(args)['olog'] = args.outputdir + '/../' + 'log.out'
    vars(args)['elog'] = args.outputdir + '/../' + 'log.err'

    return args

# -------------------------------------

def main():
    """Main function"""

    # get arguments
    args = get_arg()

    hp.echostep(args.step)

    # print args
    print(args)
    print

    # mkdir -p
    hp.mkdirp(args.outputdir)

    # end of step
    hp.echostep(args.step, start=0)

# -------------------------------------

if __name__ == '__main__':

    main()
