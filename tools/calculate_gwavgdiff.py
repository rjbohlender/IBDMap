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
                    l = l.strip().split()
                    if l[0] not in self.data:
                        self.data[l[0]] = OrderedDict()
                    self.data[l[0]][int(l[1])] = float(l[2])

    def find_nearest(self, chrom, pos):
        """Find the nearest genetic map value and return a pair of values
           to be interpolated.
        """
        if pos in self.data[chrom]:
            return self.data[chrom][pos], self.data[chrom][pos]
        else:
            nearest = (abs(np.array(list(self.data[chrom].keys())) - breakpoint)).argmin()
            near_pos = self.data[chrom].keys()[nearest]
            if near_pos > pos:  # to the right
                if nearest > 0:
                    other_pos = self.data[chrom].keys()[nearest - 1]
                    return self.data[chrom][other_pos], self.data[chrom][near_pos]
                else:
                    other_pos = self.data[chrom].keys()[nearest + 1]
                    return self.data[chrom][near_pos], self.data[chrom][other_pos]
            elif nearest == len(self.data[chrom].keys()) - 1:
                other_pos = self.data[chrom].keys()[nearest - 1]
                return self.data[chrom][other_pos], self.data[chrom][near_pos]
            else:
                other_pos = self.data[chrom].keys()[nearest + 1]
                return self.data[chrom][near_pos], self.data[chrom][other_pos]


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


def parse_file(pheno, ifile, counts):
    """Parse the ibd file"""
    casecaserate = 0
    casecontrolrate = 0
    controlcontrolrate = 0
    with ifile.open('r') as f:
        breakpoints = 0
        for l in f:
            if l.startswith('#'):
                continue
            breakpoints += 1

            l = l.strip().split()
            if len(l) > 2:
                casecase = 0
                casecontrol = 0
                controlcontrol = 0
                samples = [x.split(':')[1].split('-') for x in l[2:]]
                try:
                    for x1, x2 in samples:
                        try:
                            if pheno[x1] == 0 and pheno[x2] == 0:
                                controlcontrol += 1
                            elif pheno[x1] == 1 and pheno[x2] == 0:
                                casecontrol += 1
                            elif pheno[x1] == 0 and pheno[x2] == 1:
                                casecontrol += 1
                            elif pheno[x1] == 1 and pheno[x2] == 1:
                                casecase += 1
                        except KeyError:
                            continue
                except ValueError:
                    print(l, file=sys.stderr)
                    print(samples, file=sys.stderr)
                    sys.exit(1)
                casecaserate += casecase / counts[0]
                casecontrolrate += casecontrol / counts[1]
                controlcontrolrate += controlcontrol / counts[2]

    return casecaserate / breakpoints, casecontrolrate / breakpoints, controlcontrolrate / breakpoints


def main():
    """Entrypoint"""
    parser = ap.ArgumentParser()
    parser.add_argument('data_dir', help="The directory containing the data.")
    parser.add_argument('file_glob', help="The glob pattern for the files.")
    parser.add_argument('phen_file', help="The phenotype file.")
    parser.add_argument('genetic_map', help="The genetic map directory.")
    args = parser.parse_args()

    gmap = GeneticMap(args.genetic_map)

    pheno = parse_pheno(Path(args.phen_file))
    ncase = sum(x for x in pheno.values())
    ncontrol = sum(x == 0 for x in pheno.values())

    ncasecase = ncase * (ncase - 1) / 2.
    ncasecontrol = ncase * ncontrol
    ncontrolcontrol = ncontrol * (ncontrol - 1) / 2.

    casecaserate = 0
    casecontrolrate = 0
    controlcontrolrate = 0

    pattern = args.data_dir + '/' + args.file_glob + "*.txt"
    files = glob.glob(pattern)
    print(pattern, file=sys.stderr)
    print(files, file=sys.stderr)
    print('\n'.join(files), file=sys.stderr)

    for f in files:
        casecase, casecontrol, controlcontrol = parse_file(pheno, Path(f), (ncasecase, ncasecontrol, ncontrolcontrol))
        casecaserate += casecase
        casecontrolrate += casecontrol
        controlcontrolrate += controlcontrol
    casecaserate /= len(files)
    casecontrolrate /= len(files)
    controlcontrolrate /= len(files)

    print("casecaserate: {}".format(casecaserate))
    print("casecontrolrate: {}".format(casecontrolrate))
    print("controlcontrolrate: {}".format(controlcontrolrate))


if __name__ == "__main__":
    main()
