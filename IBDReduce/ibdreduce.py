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
from geneticmap import GeneticMap
from pathlib import Path
import numpy as np
import scipy.stats as stats


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


def ibdlen_parse(i: int, gmap: GeneticMap, args: ap.Namespace) -> Dict[int, dict]:
    """First pass of parsing for the distribution of distances between markers and length of the genome.

    Args:
        i: The chromosome number.
        gmap: The genetic map used to calculate distances.
        args: The program arguments.

    Returns:
        A dictionary of genetic distances between markers for each chromosome.
    """
    ibdlen = {i: {}}
    file_ = Path(args.prefix + args.suffix.format(i=i, j=args.at))

    if is_gzipped(file_):
        f = gzip.open(str(file_), 'rt')
    else:
        f = file_.open('r')

    predis = 0
    prechr = 0
    prepos = 0
    for lineno, line in enumerate(f):
        chrom, orig, pos, vals = parse_line(line, args.new)

        if pos in gmap.gmap[chrom]:
            dis = gmap.gmap[chrom][pos]
        else:
            close = gmap.find_nearest(chrom, pos)
            dis = gmap.gmap[chrom][close[0]] + ((int(pos) - close[0]) / (close[1] - close[0])) * (
                    gmap.gmap[chrom][close[1]] - gmap.gmap[chrom][close[0]])
        if chrom == prechr:
            ibdlen[prechr][prepos] = float(dis) - float(predis)
        elif prechr != 0:
            ibdlen[prechr][prepos] = float(predis) + 1  # for final variants
        prechr = chrom
        prepos = pos
        predis = dis
    ibdlen[prechr][prepos] = 1  # for final variants
    return ibdlen


def parse_avg(i: int, args: ap.Namespace, ibd_frac: Dict[int, dict]) \
        -> Tuple[List[Dict[int, dict]], List[Dict[int, np.ndarray]]]:
    """Second pass of parsing for the genome wide average of the statistic across markers.

    The calculation is done for the original statistic and each permutation. Average is calculated across markers.

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
            if j == 1:
                if args.two_sided:
                    original[lineno % args.phenotypes][chrom][pos] = np.abs(orig)
                    permuted[lineno % args.phenotypes][chrom][0:args.nperm] += np.abs(vals)
                else:
                    original[lineno % args.phenotypes][chrom][pos] = orig
                    permuted[lineno % args.phenotypes][chrom][0:args.nperm] += vals
            else:
                try:
                    if args.two_sided:
                        permuted[
                            lineno % args.phenotypes][chrom][
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


def parse(i: int, orig_avg: List[float], permuted_avg: List[np.ndarray], args: ap.Namespace) \
        -> Tuple[List[Dict[int, dict]], List[Dict[int, np.ndarray]]]:
    """Final pass parser for the extreme value distribution.

    Args:
        i: The chromosome number.
        orig_avg: The average of the observed statistics.
        permuted_avg: The average of the permuted statistics.
        args: The program arguments.

    Returns:
        The original statistics for each marker, converted to empirical p-value, and the extreme value distribution of
        statistics for the current chromosome.
    """
    original = []
    extreme_val_dist = []
    arr_size = args.nruns * args.nperm
    for k in range(args.phenotypes):
        original.append({i: {}})
        extreme_val_dist.append({i: np.zeros(arr_size)})
        extreme_val_dist[k][i].fill(float('Inf'))

    files = []
    for j in range(args.at, args.at + args.nruns):
        files.append(check_and_open(args.prefix + args.suffix.format(i=i, j=j)))
    lineno = -1
    evd_buffer = np.zeros(arr_size)
    line = None
    while True:
        lineno += 1
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

            succ = np.sum(vals >= orig)
            try:
                original[lineno % args.phenotypes][chrom][pos] = ((succ + 1) / (len(vals) + 1), succ)
            except Exception as e:
                print(f'i: {i}, chrom: {chrom}, pos: {pos}')
                raise e
            evd_buffer[args.nperm * j:args.nperm * (j + 1)] = vals
        if line == '':
            break
        evd_buffer = get_pdist(evd_buffer, method=args.method)
        extreme_val_dist[lineno % args.phenotypes][i] = \
            np.minimum(extreme_val_dist[lineno % args.phenotypes][i], evd_buffer)
    return original, extreme_val_dist


def run_ibdlen(args: ap.Namespace, gmap: GeneticMap) -> Tuple[Dict[int, Dict[int, float]], Dict[int, Dict[int, float]]]:
    """Encapsulate running of IBD length parsing and calculations.

    Args:
        args: Program arguments.
        gmap: The recombination map.

    Returns:
        Mapping of chromosomes and positions against genetic distances between positions, and fractions of total
        distance between variants.
    """
    ibdlen = {}
    future = []
    if args.single:
        chroms = [args.single]
        if args.null:
            chroms.append(args.null)
    else:
        chroms = list(range(1, 23))
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in chroms:
            future.append(executor.submit(ibdlen_parse, i, gmap, args))
        for fut in future:
            res = fut.result()
            ibdlen.update(res)

    sum_ibdlen = 0
    for i in chroms:
        for k, v in ibdlen[i].items():
            sum_ibdlen += ibdlen[i][k]

    ibd_frac = {}
    for i in chroms:
        ibd_frac[i] = {}
        for k, v in ibdlen[i].items():
            ibd_frac[i][k] = v / sum_ibdlen

    return ibdlen, ibd_frac


def run_avg(args: ap.Namespace, ibd_frac: Dict[int, Dict[int, float]]) \
        -> Tuple[List[Dict[int, Dict[int, float]]], List[np.ndarray]]:
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


def run_parse(args: ap.Namespace, original_avg: List[Dict[int, Dict[int, float]]], p_avg: List[np.ndarray]) \
        -> Tuple[List[dict], List[np.ndarray]]:
    """Encapsulate running of the final parsing pass and calculation of the extreme value distribution.

    Args:
        args: The program arguments.
        original_avg: The average of the observed statistics.
        p_avg: The averages of the permuted statistics.

    Returns:
        The observed statistics and the extreme value distribution.
    """
    arr_size = args.nruns * args.nperm
    original = [{} for _ in range(args.phenotypes)]
    extreme_val_dist = [{} for _ in range(args.phenotypes)]
    future = []
    if args.single:
        chroms = [args.single]
    else:
        chroms = list(range(1, 23))
    with concurrent.futures.ThreadPoolExecutor() as executor:
        for i in chroms:
            future.append(executor.submit(parse, i, original_avg, p_avg, args))
        for fut in future:
            res = fut.result()
            for j in range(args.phenotypes):
                original[j].update(res[0][j])
                extreme_val_dist[j].update(res[1][j])

    # Flatten dict across chromosomes
    evd = []
    for phen in range(args.phenotypes):
        evd.append(np.zeros(arr_size))
        evd[phen].fill(float('Inf'))
        for chrom, dist in extreme_val_dist[phen].items():
            evd[phen] = np.minimum(dist, evd[phen])
    return original, evd


def main():
    """Entrypoint."""
    parser = ap.ArgumentParser()
    parser.add_argument('prefix', type=str, help="File prefix, to which 'chr{i}.txt.{j}.results' will be appended")
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
    args = parser.parse_args()

    gmaps = glob.glob(args.gmap + '/*.gmap.gz')
    gmap = GeneticMap(gmaps)

    ibdlen, ibd_frac = run_ibdlen(args, gmap)
    original_avg, p_avg = run_avg(args, ibd_frac)
    original, evd = run_parse(args, original_avg, p_avg)

    if args.single:
        chroms = [args.single]
    else:
        chroms = list(range(1, 23))
    # Collapse averages, need dot(avgs, ibdlen)
    total_perms = args.nperm * args.nruns
    for phen in range(args.phenotypes):
        opf = open(f'{args.output}.{phen}', 'w')
        print("CHROM\tPOS\tcM\tPVal\tPValCI\tPAdj\tPAdjCutoff\tSuccess\tPermutation", file=opf)

        # Significance level for p-value in permutation correcting for multiple tests
        p_adjust_cutoff = np.percentile(evd[phen], 5.)
        for i in chroms:
            for k, res in original[phen][i].items():
                p = res[0]
                succ = res[1]
                # Poisson CI
                upper = stats.chi2.ppf(0.975, 2 * succ + 2) / 2.
                if succ > 0:
                    lower = stats.chi2.ppf(0.025, 2 * succ) / 2.
                else:
                    lower = 0.

                # Adjusted p-value
                p_adjust = stats.percentileofscore(evd[phen], p, kind='weak') / 100.
                print(f"{i}\t{k}\t{ibdlen[i][k]}\t{p}\t{lower}," +
                      f"{upper}\t{p_adjust}\t{p_adjust_cutoff}\t{succ}\t{total_perms}", file=opf)


if __name__ == "__main__":
    main()
