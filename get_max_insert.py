#!/usr/bin/env python


import sys
import re
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--std-devs', dest='std_devs', default=5)
    parser.add_argument('infile', nargs='?', type=argparse.FileType('r'),
            default=sys.stdin)
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'),
            default=sys.stdout)
    args = parser.parse_args()

    std_devs = int(args.std_devs)
    infile = args.infile
    outfile = args.outfile

    median_exp = re.compile(r'Actual FR median insert size:\s*(\d+)')
    mad_exp = re.compile(r'Actual FR median absolute deviation:\s*(\d+)')

    for line in infile:
        if median_exp.match(line):
            median = int(median_exp.match(line).group(1))
            mad = int(mad_exp.match(infile.next()).group(1))
            break

    max_insert = median + (std_devs * mad)

    outfile.write("%d" % max_insert)
