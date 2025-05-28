"""Combine multiple runs of IBDMap.

Requires Python >= 3.6
"""

import glob
import argparse as ap
import sys
import gzip
from typing import Dict, List, Tuple
from pathlib import Path
import numpy as np
import scipy.stats as stats
import ibdlib
from datetime import datetime


def is_gzipped(fpath: Path) -> bool:
    """Check if the provided file is gzipped or not.

    Args:
        fpath: The path to the file to check.

    Returns:
        A bool indicating if the file is gzipped (true) or not (false).
    """

    with fpath.open('rb') as f:
        magic = f.read(2)

    return magic == b'\x1f\x8b'


def check_and_open(fpath: str):
    """Check if a file is gzipped and return an open file handle.

    Args:
        fpath: The path to the file.
    Returns:
        A file handle pointing to the file.
    """

    file_ = Path(fpath)

    if is_gzipped(file_):
        f = gzip.open(str(file_), 'rt')
    else:
        f = file_.open('r')
    return f


def get_pdist(stat_dist: np.ndarray, method: str) -> np.ndarray:
    """Convert statistics to percentiles.

    Args:
        stat_dist: Distribution of statistics.
        method: The method to be used, passed to stats.rankdata. One of average, min, or max.

    Returns:
        The distribution of values converted to percentiles.
    """
    return stats.rankdata(stat_dist, method=method) / (len(stat_dist) + 1)


def main():
    """Entrypoint."""
    global args
    parser = ap.ArgumentParser()
    parser.add_argument('prefix',
                        type=str,
                        help="File prefix, to which suffix will be appended to form a complete file name")
    parser.add_argument('suffix', type=str,
                        help="File suffix. Use {i} for the chromosome, {j} for the permutation set.")
    parser.add_argument('nperm', type=int, help="Total number of permutations in each run.")
    parser.add_argument('nruns', type=int, help="Total number of runs.")
    parser.add_argument('gmap', type=str, help="Path to genetic map files.")
    parser.add_argument(
        '--method',
        type=str,
        default='average',
        help="Method for calculating p-value. One of {average, min, max}, default is average.")
    parser.add_argument('--at', default=0, type=int,
                        help="The starting index of the output files. Starting value of j.")
    parser.add_argument('--single', default=None, type=int, help="Run only a single chromosome.")
    parser.add_argument('--null', default=None, type=int,
                        help="Run an alternate single chromosome for the null distribution.")
    parser.add_argument('--no_avg', default=False, action='store_true', help="Don't calculate the genomewide average.")
    parser.add_argument('--two_sided', default=False, action='store_true', help="Calculate for a two-sided test.")
    parser.add_argument('--phenotypes', default=1, type=int, help="Number of phenotypes we're parsing.")
    parser.add_argument('--output', required=True, help="Output path. will be appended with number.")
    parser.add_argument('--separation', type=int, help="Amount of space between printed markers. Units are basepairs.")
    parser.add_argument('--new', default=False, action='store_true',
                        help="Does the output have the proportion of pairs in 3 columns?")
    parser.add_argument('--fdr', default=False, action='store_true',
                        help="Control FDR instead of FWE.")
    args = parser.parse_args()

    cppargs = ibdlib.IBDRArgs()
    cppargs.at = args.at
    cppargs.phenotypes = args.phenotypes
    cppargs.nperm = args.nperm
    cppargs.nruns = args.nruns
    cppargs.new_ = args.new
    cppargs.fdr = args.fdr
    cppargs.prefix = args.prefix
    cppargs.suffix = args.suffix
    if args.single is not None:
        cppargs.is_single = True
        cppargs.single = args.single
    else:
        cppargs.is_single = False
    if args.null is not None:
        cppargs.is_null = True
        cppargs.null = args.null
    else:
        cppargs.is_null = False

    t1 = datetime.now()
    gmaps = glob.glob(args.gmap + '/*.gmap.gz')
    gmap = ibdlib.GeneticMap(gmaps)
    t2 = datetime.now()
    print("GMAP time v2: {}".format(t2 - t1), file=sys.stderr)

    t1 = datetime.now()
    original, evd, deltas, ibdlen, fdr = ibdlib.run_stages(cppargs, gmap)  # Runtime J * M * N / K * (5N + 5)
    # ibdlib.run_stages(cppargs, gmap)  # Runtime J * M * N / K * (5N + 5)
    t2 = datetime.now()
    print("run_stages time v2: {}".format(t2 - t1), file=sys.stderr)

    if args.single:
        chroms = [args.single]
    else:
        chroms = list(range(1, 23))

    if args.fdr:
        for phen in range(args.phenotypes):
            fdr[phen] = stats.rankdata(fdr[phen], method='max', axis=1)
            fdr[phen] = np.divide(fdr[phen], args.nruns * args.nperm + 1)

    # Collapse averages, need dot(avgs, ibdlen)
    total_perms = args.nperm * args.nruns
    opath = Path('/'.join(args.output.split('/')[0:len(args.output.split('/')) - 1]))
    opath.mkdir(parents=True, exist_ok=True)

    t1 = datetime.now()
    for phen in range(args.phenotypes):
        opf = open(f'{args.output}.{phen}', 'w')
        print('# {}'.format(' '.join(sys.argv)), file=opf)
        print("CHROM\tPOS\tcM\tPVal\tPValCI\tPAdj\tPAdjCutoff\tSuccess\tPermutation\tDelta", file=opf)

        # Significance level for p-value in permutation correcting for multiple tests
        p_adjust_cutoff = np.percentile(evd[phen], 5.)
        memoize = {}  # See if we can save time by memoizing all our Rstar lookups
        for i in chroms:
            for k, res in original[phen][i].items():
                p = res[0]
                succ = res[1]
                d = deltas[phen][i][k]
                # Poisson CI
                upper = stats.chi2.ppf(0.975, 2 * succ + 2) / 2.
                if succ > 0:
                    lower = stats.chi2.ppf(0.025, 2 * succ) / 2.
                else:
                    lower = 0.

                # Adjusted p-value
                if args.fdr:
                    if p in memoize:
                        Rstar = memoize[p]
                    else:
                        Rstar = np.sum(fdr[phen] <= p, axis=1)
                        memoize[p] = Rstar
                    rp = sum(res[0] <= p for _, res in original[phen][i].items())
                    pm = p * fdr[phen].shape[0]
                    rb = np.percentile(Rstar, 95)
                    if rp - rb >= pm:
                        p_adjust = np.mean(Rstar / (Rstar + rp - pm))
                    else:
                        p_adjust = stats.percentileofscore(evd[phen], p, kind='weak') / 100.
                else:
                    p_adjust = stats.percentileofscore(evd[phen], p, kind='weak') / 100.
                print(f"{i}\t{k}\t{ibdlen[i][k]}\t{p}\t{lower}," +
                      f"{upper}\t{p_adjust}\t{p_adjust_cutoff}\t{succ}\t{total_perms}\t{d}", file=opf)
    t2 = datetime.now()
    print("output time v2: {}".format(t2 - t1), file=sys.stderr)


if __name__ == "__main__":
    main()
