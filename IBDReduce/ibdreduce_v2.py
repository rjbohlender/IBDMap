"""Combine multiple runs of carvaIBD.

Requires Python >= 3.6

First line is the value of the statistic at the locus without the average
across permutations. The second line is the value of the average difference at
the position across permutations, scaled to the total number of breakpoints.
"""

import glob
import argparse as ap
import sys
import concurrent.futures
import gzip
from typing import Dict, List, Tuple
from pathlib import Path
import numpy as np
import scipy.stats as stats
import ibdlib
from datetime import datetime
import multiprocessing as mp

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


def parse_avg(i: int, args: ap.Namespace, ibd_frac: Dict[int, dict]) \
        -> Tuple[List[Dict[int, Dict[int, float]]], List[Dict[int, np.ndarray]]]:
    """Second pass of parsing for the genome wide average of the statistic across markers.

    The calculation is done for the original statistic and each permutation. Average is calculated across markers.

    Runtime:

    Args:
        i: The chromosome number.
        args: The program arguments.
        ibd_frac: The distribution of genetic distances between markers.

    Returns:
        The average of original values by chromsome for each phenotype, and the average of each permutation by
        chromosome.
    """
    arr_size = args.nruns * args.nperm
    permuted = [{i: np.zeros(arr_size)} for _ in range(args.phenotypes)]
    original = [{i: {}} for _ in range(args.phenotypes)]


    #    start = min(list(range(args.at, args.at + args.nruns)))
    start = args.at
    for j in range(args.at, args.at + args.nruns):
        file_ = Path(args.prefix + args.suffix.format(i=i, j=j))

        if is_gzipped(file_):
            f = gzip.open(str(file_), 'rt')
        else:
            f = file_.open('r')

        for lineno, l in enumerate(f):
            chrom, orig, pos, vals = parse_line(l, args.new)

            if not args.no_avg:
                orig *= ibd_frac[chrom][pos]
                vals *= ibd_frac[chrom][pos]

            if np.all(vals == 0):
                continue
            if j == args.at:
                if args.two_sided:
                    original[lineno % args.phenotypes][chrom][pos] = np.abs(orig)
                    permuted[lineno % args.phenotypes][chrom][0:args.nperm] += np.abs(vals)
                else:
                    original[lineno % args.phenotypes][chrom][pos] = orig
                    permuted[lineno % args.phenotypes][chrom][0:args.nperm] += vals
            else:
                try:
                    if args.two_sided:
                        permuted[lineno % args.phenotypes][chrom][
                        args.nperm * (j - args.at):args.nperm * (j - args.at + 1)] += np.abs(vals)
                    else:
                        permuted[lineno % args.phenotypes][chrom][
                        (args.nperm * (j - args.at)):(args.nperm * (j - args.at + 1))] += vals
                except ValueError as err:
                    print(
                        'permuted: {} vals: {}'.format(permuted[lineno % args.phenotypes][chrom].shape, vals.shape),
                        file=sys.stderr)
                    raise err
                except KeyError as err:
                    print(permuted, file=sys.stderr)
                    print(f'chrom: {chrom} pos: {pos} j: {j}', file=sys.stderr)
                    raise err
    return original, permuted


def parse_line(line: str, new: bool = False) -> Tuple[int, float, int, np.ndarray]:
    """Line parser.

    Example line pair:
    22	16504399	0.00034038	8.31453e-07	-2.27348e-07

    Args:
        line: The line to be split and processed.
        new: Whether we have the new format or not.

    Returns:
        Tuple of useful values from the line.
    """
    line = line.strip().split()
    assert len(line) >= 2

    orig_idx = 2
    line_start = 3
    if new:
        orig_idx = 5
        line_start = 6
    if line[0].startswith('chr'):
        chrom = int(line[0][3:])
    else:
        chrom = int(line[0])
    pos = int(line[1])
    orig = float(line[orig_idx])
    vals = np.array([float(x) for x in line[line_start:]])
    return chrom, orig, pos, vals


def get_pdist(stat_dist: np.ndarray, method: str) -> np.ndarray:
    """Convert statistics to percentiles.

    Args:
        stat_dist: Distribution of statistics.
        method: The method to be used, passed to stats.rankdata. One of average, min, or max.

    Returns:
        The distribution of values converted to percentiles.
    """
    return stats.rankdata(stat_dist, method=method) / (len(stat_dist) + 1)


def parse(i: int, orig_avg: List[float], permuted_avg: List[np.ndarray], breakpoints: Dict[int, int],
          args: ap.Namespace) \
        -> Tuple[List[Dict[int, dict]],
                 List[Dict[int, np.ndarray]],
                 List[Dict[int, np.ndarray]],
                 List[Dict[int, Dict[int, float]]]]:
    """Final pass parser for the extreme value distribution.

    Args:
        i: The chromosome number.
        orig_avg: The average of the observed statistics.
        permuted_avg: The average of the permuted statistics.
        breakpoints: The number of breakpoints for each chromosome.
        args: The program arguments.

    Returns:
        The original statistics for each marker, converted to empirical p-value, and the extreme value distribution of
        statistics for the current chromosome, as well as the full distribution of resampled statistics if FDR argument
        is passed..
    """
    original = []
    extreme_val_dist = []
    deltas = []
    arr_size = args.nruns * args.nperm
    fdr = []
    for k in range(args.phenotypes):
        original.append({i: {}})
        extreme_val_dist.append({i: np.zeros(arr_size)})
        extreme_val_dist[k][i].fill(float('Inf'))
        deltas.append({i: {}})
        if args.fdr:
            # Corresponds to R* in Yekutieli and Benjamini (1999)
            fdr.append({i: np.zeros((breakpoints[i], arr_size))})
        else:
            fdr.append({i: np.array([])})

    orig_buffer = [{i: {}} for _ in range(args.phenotypes)]

    files = []
    for j in range(args.at, args.at + args.nruns):
        files.append(check_and_open(args.prefix + args.suffix.format(i=i, j=j)))
    lineno = -1
    breakpointno = -1
    evd_buffer = np.zeros(arr_size)
    line = None
    while True:
        lineno += 1
        if lineno % args.phenotypes == 0:
            breakpointno += 1
        for j, f in enumerate(files):
            line = f.readline()
            if line == '':
                continue
            chrom, orig, pos, vals = parse_line(line, args.new)
            if all(vals == 0):
                continue
            if not args.no_avg:
                orig -= orig_avg[lineno % args.phenotypes]
                vals -= permuted_avg[lineno % args.phenotypes][args.nperm * j:args.nperm * (j + 1)]

            deltas[lineno % args.phenotypes][chrom][pos] = orig
            if pos not in orig_buffer[lineno % args.phenotypes][chrom]:
                orig_buffer[lineno % args.phenotypes][chrom][pos] = np.array([0, 0])

            orig_buffer[lineno % args.phenotypes][chrom][pos] += np.array([np.sum(vals >= orig), len(vals)])
            evd_buffer[args.nperm * j:args.nperm * (j + 1)] = vals
            if args.fdr:
                fdr[lineno % args.phenotypes][i][breakpointno, args.nperm * j:args.nperm * (j + 1)] = vals
        if line == '':
            break
        evd_buffer = get_pdist(evd_buffer, method=args.method)
        extreme_val_dist[lineno % args.phenotypes][i] = \
            np.minimum(extreme_val_dist[lineno % args.phenotypes][i], evd_buffer)
    for lineno, phen in enumerate(orig_buffer):
        for chrom, chrom_data in phen.items():
            for pos, data in chrom_data.items():
                original[lineno][chrom][pos] = ((data[0] + 1) / (data[1] + 1), data[0])
    return original, extreme_val_dist, fdr, deltas


def run_avg(args: ap.Namespace, ibd_frac: Dict[int, Dict[int, float]]) \
        -> Tuple[List[float], List[np.ndarray]]:
    """Encapsulate running of the averaging process.

    Observed and permuted statistics are averaged along columns, with each variant weighted by the fraction of total
    genetic distance between the variant and the prior variant.

    Args:
        args: Program arguments.
        ibd_frac: The distribution of genetic distances as weights for averaging DeltaIBD statistics.

    Returns:
        The average of the observed statistics, and the average of the permuted statistics, weighted by distance
        between markers.
    """
    arr_size = args.nruns * args.nperm
    permuted_avg = [{} for _ in range(args.phenotypes)]
    original_avg = [{} for _ in range(args.phenotypes)]
    if args.single:
        if args.null:
            chroms = [args.null]
        else:
            chroms = [args.single]
    else:
        chroms = list(range(1, 23))

    future = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in chroms:
            future.append(executor.submit(parse_avg, i, args, ibd_frac))
        for fut in future:
            res = fut.result()
            for j in range(args.phenotypes):
                original_avg[j].update(res[0][j])
                permuted_avg[j].update(res[1][j])
    original_avg = [
        sum(sum(v for _, v in original_avg[j][i].items()) for i in chroms) for j in range(args.phenotypes)]
    p_avg = [np.zeros(arr_size) for _ in range(args.phenotypes)]
    for i, phen_set in enumerate(permuted_avg):
        for _, v in phen_set.items():
            p_avg[i] += v
    return original_avg, p_avg


def run_parse(args: ap.Namespace, original_avg: List[Dict[int, Dict[int, float]]], p_avg: List[np.ndarray],
              breakpoints: Dict[int, int]) \
        -> Tuple[List[dict], List[np.ndarray], List[np.ndarray], List[Dict[int, Dict[int, float]]]]:
    """Encapsulate running of the final parsing pass and calculation of the extreme value distribution.

    Args:
        args: The program arguments.
        original_avg: The average of the observed statistics.
        p_avg: The averages of the permuted statistics.
        breakpoints: The number of breakpoints on each chromosome.

    Returns:
        The observed statistics and the extreme value distribution.
    """
    arr_size = args.nruns * args.nperm
    original = [{} for _ in range(args.phenotypes)]
    extreme_val_dist = [{} for _ in range(args.phenotypes)]
    future = []
    deltas = [{} for _ in range(args.phenotypes)]
    if args.single:
        chroms = [args.single]
    else:
        chroms = list(range(1, 23))
    fdr = [{i: [] for i in chroms} for _ in range(args.phenotypes)]
    total_breakpoints = 0
    for i in chroms:
        total_breakpoints += breakpoints[i]
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in chroms:
            future.append(executor.submit(parse, i, original_avg, p_avg, breakpoints, args))
        for i, fut in enumerate(future):
            res = fut.result()
            for j in range(args.phenotypes):
                original[j].update(res[0][j])
                extreme_val_dist[j].update(res[1][j])
                fdr[j].update(res[2][j])
                deltas[j].update(res[3][j])

    # Flatten dict across chromosomes
    evd = []
    full_dist = [np.zeros((total_breakpoints, arr_size)) for _ in range(args.phenotypes)]
    for phen in range(args.phenotypes):
        last = 0
        evd.append(np.zeros(arr_size))
        evd[phen].fill(float('Inf'))
        for chrom, dist in extreme_val_dist[phen].items():
            evd[phen] = np.minimum(dist, evd[phen])
        if args.fdr:
            for chrom in chroms:
                dist = fdr[phen][chrom]
                full_dist[phen][last:(last + dist.shape[0]), :] = dist
                last += dist.shape[0]

    return original, evd, full_dist, deltas

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
    pool = mp.Pool(mp.cpu_count())

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
