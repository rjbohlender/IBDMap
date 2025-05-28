"""
Ryan Bohlender 2019

Generate the genome-wide average difference between case-case and case-control or control-control pairs.
"""

import sys
import glob
import argparse as ap
import gzip
import numpy as np
from collections import OrderedDict
from pathlib import Path


class GeneticMap(object):
    """Container for the genetic map"""
    def __init__(self, dir_path):
        self.data = OrderedDict()
        files = glob.glob(dir_path + '*.gz')
        for path in files:
            with gzip.open(path, 'rt') as f:
                for l in f:
                    # Order = pos chr cM
                    if l.startswith('pos'):  # Skip header
                        continue
                    l = l.strip().split()
                    if l[1] not in self.data:
                        self.data[l[1]] = OrderedDict()
                    self.data[l[1]][int(l[0])] = float(l[2])

    def find_nearest(self, chrom, pos):
        """Find the nearest genetic position and return a pair of values
           to be interpolated.
        """
        nearest = (abs(np.array(list(self.data[chrom].keys())) - pos)).argmin()
        keys = list(self.data[chrom].keys())
        near_pos = keys[nearest]
        if near_pos > pos:  # to the right
            if nearest > 0:
                other_pos = keys[nearest - 1]
                return other_pos, near_pos
            else:
                other_pos = keys[nearest + 1]
                return near_pos, other_pos
        elif nearest == len(keys) - 1:
            other_pos = keys[nearest - 1]
            return other_pos, near_pos
        else:
            other_pos = keys[nearest + 1]
            return near_pos, other_pos


def parse_pheno(ifile):
    """Parse the input phenotype file"""
    data = {}
    with ifile.open('r') as f:
        i = 0
        for l in f:
            l = l.strip().split()
            if i == 0:
                i += 1
                continue
            if l[1] == "NA":
                continue
            data[l[0]] = int(l[1])

    return data


def parse_file(ifiles, gmap):
    """Parse the ibd file"""
    pairmap = {}
    breakpoints = 0
    total_dis = 0
    for ifile in ifiles:
        with open(ifile, 'r') as f:
            predis = 0
            for l in f:
                if l.startswith('#'):
                    continue
                breakpoints += 1

                l = l.strip().split()
                chrom = l[0]
                pos = int(l[1])

                if pos in gmap.data[chrom]:
                    dis = gmap.data[chrom][pos]
                else:
                    pair = gmap.find_nearest(chrom, pos)
                    dis = gmap.data[chrom][pair[0]] + (pos - pair[0]) / (pair[1] - pair[0]) * (gmap.data[chrom][pair[1]] - gmap.data[chrom][pair[0]])
                curlen = dis - predis
                predis = dis
                total_dis += dis
                if len(l) > 2:
                    samples = [x.split(':')[1].split('-') for x in l[2:]]
                    try:
                        for x1, x2 in samples:
                            try:
                                if x1 in pairmap:
                                    if x2 in pairmap[x1]:
                                        pairmap[x1][x2] += 1 * curlen
                                    else:
                                        pairmap[x1][x2] = 1 * curlen
                                else:
                                    if x2 in pairmap:
                                        if x1 in pairmap[x2]:
                                            pairmap[x2][x1] += 1 * curlen
                                        else:
                                            pairmap[x2][x1] = 1 * curlen
                                    else:
                                        pairmap[x1] = {}
                                        pairmap[x1][x2] = 1 * curlen
                            except KeyError:
                                continue
                    except ValueError:
                        print(l, file=sys.stderr)
                        print(samples, file=sys.stderr)
                        sys.exit(1)

    for k1 in pairmap.keys():
        for k2 in pairmap[k1].keys():
            pairmap[k1][k2] /= total_dis

    return pairmap


def main():
    """Entrypoint"""
    parser = ap.ArgumentParser()
    parser.add_argument('data_dir', help="The directory containing the data.")
    parser.add_argument('file_glob', help="The glob pattern for the files.")
    parser.add_argument('genetic_map', help="The genetic map directory.")
    args = parser.parse_args()

    gmap = GeneticMap(args.genetic_map)

    pattern = args.data_dir + '/' + args.file_glob + "*.txt"
    files = glob.glob(pattern)
    print(pattern, file=sys.stderr)
    print(files, file=sys.stderr)
    print('\n'.join(files), file=sys.stderr)

    pairmap = parse_file(files, gmap)

    print('id1\tid2\tavg')
    for k1 in pairmap.keys():
        for k2 in pairmap[k1].keys():
            v = pairmap[k1][k2]
            print(f'{k1}\t{k2}\t{v}')


if __name__ == "__main__":
    main()
