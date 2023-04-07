import glob
import argparse as ap
import sys
from itertools import starmap
import multiprocessing as mp
import gzip
import zstandard as zstd
from typing import Dict, List, Tuple
from geneticmap import GeneticMap
from pathlib import Path
from datetime import datetime
import numpy as np
import scipy.stats as stats


def is_compressed(fpath: Path) -> str:
    """Check if the provided file is compressed with gzip or zstd or not.

    Args:
        fpath: The path to the file to check.

    Returns:
        A str indicating if the file is compressed (zstd, gzip) or not (no)..
    """

    gzip = b'\x1f\x8b'
    zstd = b'\x28\xb5\x2f\xfd'

    with fpath.open('rb') as f:
        magic = f.read(4)

    if magic[:2] == gzip:
        return 'gzip'
    elif magic == zstd:
        return 'zstd'
    else:
        return 'no'


def check_and_open(fpath: str):
    """Check if a file is gzipped and return an open file handle.

    Args:
        fpath: The path to the file.
    Returns:
        A file handle pointing to the file.
    """

    file_ = Path(fpath)

    compressed = is_compressed(file_)
    if compressed == 'gzip':
        f = gzip.open(str(file_), 'rt')
    elif compressed == 'zstd':
        f = zstd.open(str(file_), 'rt')
    else:
        f = file_.open('r')
    return f


def parse_line(line: str) -> Tuple[int, int, np.ndarray]:
    """Line parser.

    Example line pair:
    22	16504399	0.00034038	8.31453e-07	-2.27348e-07

    Args:
        line: The line to be split and processed.

    Returns:
        Tuple of useful values from the line.
    """
    line = line.strip().split()
    assert len(line) >= 2

    line_start = 5
    if line[0].startswith('chr'):
        chrom = int(line[0][3:])
    else:
        chrom = int(line[0])
    pos = int(line[1])
    vals = np.array([float(x) for x in line[line_start:]], dtype=np.float64)
    return chrom, pos, vals


def ibdlen(fpath: str, gmap: GeneticMap) -> Tuple[int, float]:
    """Calculate the IBD length for a single chromosome and permutation set.
    :param fpath: The file path to the IBDmap file.
    :param gmap: The genetic map object.
    :return: The number of breakpoints and the total length of covered physical distance.
    """
    breakpoints = 0
    total = 0
    predis = 0
    with check_and_open(fpath) as f:
        for line in f:
            chrom, pos, vals = parse_line(line)
            breakpoints += 1

            if pos in gmap.gmap[chrom]:
                dis = gmap.gmap[chrom][pos]
            else:
                close = gmap.find_nearest(chrom, pos)
                dis = gmap.gmap[chrom][close[0]] + ((int(pos) - close[0]) / (close[1] - close[0])) * (
                        gmap.gmap[chrom][close[1]] - gmap.gmap[chrom][close[0]])
            total += float(dis) - float(predis)
            predis = dis
    return breakpoints, total


def main():
    """Main function."""
    parser = ap.ArgumentParser(description='Summarize IBDmap files.')
    parser.add_argument('--prefix',
                        type=str,
                        required=True,
                        help="File prefix, to which suffix will be appended to form a complete file name")
    parser.add_argument('--suffix', type=str, required=True,
                        help="File suffix. Use {i} for the chromosome, {j} for the permutation set.")
    parser.add_argument('--nperm', type=int, required=True, help="Total number of permutations in each run.")
    parser.add_argument('--nruns', type=int, required=True, help="Total number of runs.")
    parser.add_argument('--gmap', type=str, required=True, help="Path to genetic map files.")
    parser.add_argument('--at', type=int, default=0, help="The permutation set to start at.")
    parser.add_argument('--unweighted', action='store_true', help='Use unweighted average.')
    parser.add_argument('--output', required=True, help="Output path.")
    parser.add_argument('--print_evd', default=False, action='store_true', help="Print the EVD to stdout.")
    args = parser.parse_args()

    ttotal1 = datetime.now()
    t1 = datetime.now()
    gmaps = glob.glob(args.gmap + '/*.gmap.gz')
    gmap = GeneticMap(gmaps)
    t2 = datetime.now()
    print("GMAP time: {}".format(t2 - t1), file=sys.stderr)

    # Loop over chromosomes to determine the number of IBD segments and the total length of IBD segments.
    t1 = datetime.now()
    breakpoints = 0
    total = 0
    map_args = [(args.prefix + args.suffix.format(i=i, j=args.at), gmap) for i in range(1, 23)]
    with mp.Pool() as pool:
        results = pool.starmap(ibdlen, map_args)
    for result in results:
        breakpoints += result[0]
        total += result[1]
    t2 = datetime.now()
    print("IBD Length Time: {}".format(t2 - t1), file=sys.stderr)
    print("Total Breakpoints: {}".format(breakpoints), file=sys.stderr)
    print("Total Length: {}".format(total), file=sys.stderr)

    # Create numpy arrays to store the data.
    ibdfrac = np.zeros(breakpoints, dtype=np.float64)
    data = np.zeros((breakpoints, args.nperm + 1), dtype=np.float64)
    bp_ids = []

    # Loop over chromosomes again to fill the arrays.
    idx = 0
    for i in range(1, 23):
        start_idx = idx
        for j in range(args.at, args.at + args.nruns):
            idx = start_idx
            fpath = args.prefix + args.suffix.format(i=i, j=j)
            offset = j - args.at
            predis = 0
            with check_and_open(fpath) as f:
                for line in f:
                    chrom, pos, vals = parse_line(line)
                    if offset == 0:
                        bp_ids.append((chrom, pos))
                        if pos in gmap.gmap[chrom]:
                            dis = gmap.gmap[chrom][pos]
                        else:
                            close = gmap.find_nearest(chrom, pos)
                            dis = gmap.gmap[chrom][close[0]] + ((int(pos) - close[0]) / (close[1] - close[0])) * (
                                    gmap.gmap[chrom][close[1]] - gmap.gmap[chrom][close[0]])
                        ibdfrac[idx] = (float(dis) - float(predis)) / total
                        predis = dis
                        data[idx, offset:] = vals
                    else:
                        data[idx, (offset * args.nperm):] = vals[1:]
                    idx += 1

    # Calculate the average and subtract it from the data.
    if args.unweighted:
        avgs = np.mean(data, axis=0)
    else:
        avgs = np.matmul(data.T, ibdfrac)

    def subtract(a):
        return a - avgs

    data = np.apply_along_axis(subtract, 1, data)
    delta = data[:, 0]

    # Calculate the empirical p-value.
    succ = np.zeros(breakpoints)
    for i in range(1, args.nperm * args.nruns + 1):
        succ += data[:, 0] >= data[:, i]
    empp = (succ + 1.) / (args.nperm * args.nruns + 1)

    # Generate the evd
    data = data[:, 1:]
    data = stats.rankdata(data, axis=1) / (args.nperm * args.nruns + 1)
    evd = np.min(data, axis=0)

    if args.print_evd:
        print("EVD:\n{}".format("\t".join(str(x) for x in evd)))

    # Calculate the adjusted p-value.
    adjp = np.zeros(breakpoints)
    cutoff = np.percentile(evd, 5)
    for i in range(breakpoints):
        adjp[i] = stats.percentileofscore(evd, empp[i], kind='weak') / 100.

    # Calculate the confidence interval.
    upper = stats.chi2.ppf(0.975, 2 * succ + 2) / 2.
    lower = stats.chi2.ppf(0.025, 2 * succ) / 2.
    lower[np.isnan(lower)] = 0

    # Generate header information and write results to file.
    opf = open(f'{args.output}', 'w')
    print('# {}'.format(' '.join(sys.argv)), file=opf)
    print('# Genome-wide Average: {}'.format(avgs[0]), file=opf)
    print('# Total breakpoints: {}'.format(breakpoints), file=opf)
    print("CHROM\tPOS\tcM\tPVal\tPValCI\tPAdj\tPAdjCutoff\tSuccess\tPermutation\tDelta", file=opf)

    idx = np.argsort(empp)

    for i in idx:
        print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(
            bp_ids[i][0],
            bp_ids[i][1],
            ibdfrac[i],
            empp[i],
            '{},{}'.format(lower[i], upper[i]),
            adjp[i],
            cutoff,
            succ[i],
            args.nperm * args.nruns,
            delta[i]
        ), file=opf)
    ttotal2 = datetime.now()
    print('Total runtime: {}'.format(ttotal2 - ttotal1), file=sys.stderr)


if __name__ == '__main__':
    main()

