import glob
import argparse as ap
import sys
import multiprocessing as mp
import gzip
import zstandard as zstd
from typing import Tuple
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
        A str indicating if the file is compressed (zstd, gzip) or not (no).
    """

    gzip_magic = b"\x1f\x8b"
    zstd_magic = b"\x28\xb5\x2f\xfd"

    with fpath.open("rb") as f:
        magic = f.read(4)

    if magic[:2] == gzip_magic:
        return "gzip"
    elif magic == zstd_magic:
        return "zstd"
    else:
        return "no"


def check_and_open(fpath: str):
    """Check if a file is gzipped and return an open file handle.

    Args:
        fpath: The path to the file.
    Returns:
        A file handle pointing to the file.
    """

    file_ = Path(fpath)

    compressed = is_compressed(file_)
    if compressed == "gzip":
        f = gzip.open(str(file_), "rt")
    elif compressed == "zstd":
        f = zstd.open(str(file_), "rt")
    else:
        f = file_.open("r")
    return f


def parse_line(line: str) -> Tuple[int, int, np.ndarray, np.ndarray]:
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
    if line[0].startswith("chr"):
        chrom = int(line[0][3:])
    else:
        chrom = int(line[0])
    pos = int(line[1])
    obs_parts = np.array([float(line[2]), float(line[3]), float(line[4])])
    vals = np.array([float(x) for x in line[line_start:]], dtype=np.float64)
    return chrom, pos, obs_parts, vals


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
            chrom, pos, obs, vals = parse_line(line)
            breakpoints += 1

            if pos in gmap.gmap[chrom]:
                dis = gmap.gmap[chrom][pos]
            else:
                close = gmap.find_nearest(chrom, pos)
                dis = gmap.gmap[chrom][close[0]] + (
                    (int(pos) - close[0]) / (close[1] - close[0])
                ) * (gmap.gmap[chrom][close[1]] - gmap.gmap[chrom][close[0]])
            total += float(dis) - float(predis)
            predis = dis
    return breakpoints, total


def parse_pheno(ifile):
    """

    :param ifile:
    :return:
    """
    data = {'0': 0, '1': 0, 'NA': 0}
    with open(ifile, 'r') as f:
        for i, l in enumerate(f):
            if i == 0:
                continue
            l = l.strip().split()
            if l[1] == 'NA':
                continue
            data[str(int(float(l[1])))] += 1
    return data['1'], data['0']


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
    parser.add_argument('--single', default=None, type=int, help="Run only a single chromosome.")
    parser.add_argument('--no_avg', default=False, action='store_true',
                        help="Don't calculate the genomewide average.")
    parser.add_argument('--null', default=None, type=int,
                        help="Run an alternate chromosome for the null dist. Requires --single. Implies --unweighted.")
    parser.add_argument('--fdr', default=False, action='store_true',
                        help="Control FDR instead of FWE.")
    parser.add_argument('--two_sided', default=False, action='store_true',
                        help="Conduct a two_sided test.")
    parser.add_argument("--lrt", default=None, help="Calculate likelihood ratio at each breakpoint. Requires phenotype file as argument.")
    args = parser.parse_args()

    if args.null and not args.single:
        print("--null requires --single.", file=sys.stderr)
        sys.exit(1)

    if args.no_avg and args.two_sided:
        print("Cannot use --no_avg and --two_sided at the same time.", file=sys.stderr)
        sys.exit(1)

    ttotal1 = datetime.now()
    t1 = datetime.now()
    gmaps = glob.glob(args.gmap + "/*.gmap.gz")
    gmap = GeneticMap(gmaps)
    t2 = datetime.now()
    print("GMAP time: {}".format(t2 - t1), file=sys.stderr)

    if args.single:
        chrom = [args.single]
        if args.null:
            chrom.append(args.null)
    else:
        chrom = list(range(1, 23))

    # Loop over chromosomes to determine the number of IBD segments and the total length of IBD segments.
    t1 = datetime.now()
    breakpoints = 0
    total = 0
    null_idx = 0
    null_breakpoints = 0
    map_args = [(args.prefix + args.suffix.format(i=i, j=args.at), gmap) for i in chrom]
    if args.null:
        null_idx = chrom.index(args.null)
    with mp.Pool(processes=min(mp.cpu_count(), len(map_args))) as pool:
        results = pool.starmap(ibdlen, map_args)
    for i, result in enumerate(results):
        breakpoints += result[0]
        if args.null and i == null_idx:
            null_breakpoints = result[0]
        total += result[1]
    t2 = datetime.now()
    print("IBD Length Time: {}".format(t2 - t1), file=sys.stderr)
    print("Total Breakpoints: {}".format(breakpoints), file=sys.stderr)
    print("Total Length: {}".format(total), file=sys.stderr)

    # Create numpy arrays to store the data.
    ibdfrac = np.zeros(breakpoints, dtype=np.float64)
    data = np.zeros((breakpoints, args.nperm * args.nruns + 1), dtype=np.float64)
    bp_ids = []

    # Loop over chromosomes again to fill the arrays.
    idx = 0
    obs_rates = np.zeros((breakpoints, 3), dtype=np.float64)
    for i in chrom:
        start_idx = idx
        if args.null and i == args.null:
            null_idx = idx
        for j in range(args.at, args.at + args.nruns):
            idx = start_idx
            fpath = args.prefix + args.suffix.format(i=i, j=j)
            offset = j - args.at
            predis = 0
            with check_and_open(fpath) as f:
                for line in f:
                    chrom, pos, obs, vals = parse_line(line)
                    try:
                        if offset == 0:
                            obs_rates[idx,:] = obs
                            bp_ids.append((chrom, pos))
                            if pos in gmap.gmap[chrom]:
                                dis = gmap.gmap[chrom][pos]
                            else:
                                close = gmap.find_nearest(chrom, pos)
                                dis = gmap.gmap[chrom][close[0]] + (
                                    (int(pos) - close[0]) / (close[1] - close[0])
                                ) * (
                                    gmap.gmap[chrom][close[1]] - gmap.gmap[chrom][close[0]]
                                )
                            ibdfrac[idx] = (float(dis) - float(predis)) / total
                            predis = dis
                            data[idx, offset:(args.nperm + 1)] = vals
                        else:
                            data[
                                idx, (offset * args.nperm):((offset + 1) * args.nperm)
                            ] = vals[1:]
                    except ValueError as e:
                        print(f"Error at {i} {j} {idx} {offset} {fpath}", file=sys.stderr)
                        raise e
                    idx += 1

    # Calculate the average and subtract it from the data.
    if args.no_avg:
        avgs = np.zeros(args.nperm * args.nruns + 1)
    elif args.unweighted:
        avgs = np.mean(data, axis=0)
    else:
        if args.null:
            avgs = np.mean(data[null_idx:(null_idx + null_breakpoints), :], axis=0)
        else:
            avgs = np.matmul(data.T, ibdfrac)

    # Likelihood ratio for localization
    llik_ratio = np.zeros(breakpoints)
    if args.lrt:
        ncase, ncontrol = parse_pheno(args.lrt)

        cscs = ncase * (ncase - 1.) / 2.
        cscn = ncase * ncontrol

        obs_counts = obs_rates[:, 0:2]

        obs_counts[:, 0] *= cscs
        obs_counts[:, 1] *= cscn
        obs_counts = obs_counts.astype(np.int64)

        weights = np.array([cscs / (cscs + cscn), cscn / (cscs + cscn)])

        mean_counts = np.sum(obs_counts * weights, 1)
        mean_counts = mean_counts.astype(np.int64)
        mean_rates = np.sum(obs_rates[:, 0:2] * weights, 1)

        for i in range(len(mean_rates)):
            null = stats.poisson.logpmf(obs_counts[i, 0], weights[0] * mean_counts[i])
            null += stats.poisson.logpmf(obs_counts[i, 1], weights[1] * mean_counts[i])
            alt = stats.poisson.logpmf(obs_counts[i, 0], obs_counts[i, 0])
            alt += stats.poisson.logpmf(obs_counts[i, 1], obs_counts[i, 1])

            llik_ratio[i] = -2 * (null - alt)

    def subtract(a):
        return a - avgs

    data = np.apply_along_axis(subtract, 1, data)
    if args.two_sided:
        data = np.abs(data)
    delta = data[:, 0]

    # Calculate the empirical p-value.
    succ = np.zeros(breakpoints)
    for i in range(1, args.nperm * args.nruns + 1):
        succ += data[:, 0] <= data[:, i]
    empp = (succ + 1.0) / (args.nperm * args.nruns + 1)

    # Generate the evd
    data = data[:, 1:]
    data = stats.rankdata(-data, axis=1) / (args.nperm * args.nruns + 1)
    evd = np.min(data, axis=0)

    if args.print_evd:
        print("EVD:\n{}".format("\t".join(str(x) for x in evd)))

    # Calculate the adjusted p-value.
    adjp = np.zeros(breakpoints)
    cutoff = np.percentile(evd, 5)
    for i in range(breakpoints):
        adjp[i] = stats.percentileofscore(evd, empp[i], kind="weak") / 100.0

    # Calculate the confidence interval.
    upper = stats.chi2.ppf(0.975, 2 * succ + 2) / 2.0
    lower = stats.chi2.ppf(0.025, 2 * succ) / 2.0
    lower[np.isnan(lower)] = 0

    upper /= data.shape[1]
    lower /= data.shape[1]

    if args.fdr:
        # Permutation p-values are likely to be repeated, so we save on computation at the cost of memory
        memoize = {}
        for i, p in enumerate(empp):
            # Rstar is the distribution of rejected null hypotheses over permutations
            if p in memoize:
                Rstar = memoize[p]
            else:
                Rstar = np.sum(data <= p, axis=0)
                memoize[p] = Rstar
            # rp is the number of rejected null hypotheses in the observed set
            rp = np.sum(data[:, 0] <= p)
            # pm is the expected number of rejected null hypotheses at a given alpha = p
            pm = p * data.shape[0]
            # rb is the 95th percentile of the Rstar distribution
            rb = np.percentile(Rstar, 95)
            # Condition for applying the fdr adjustment
            if rp - rb >= pm:
                adjp[i] = np.mean(Rstar / (Rstar + rp - pm))
            else:
                adjp[i] = sum(Rstar >= 1) / len(Rstar)

    # Generate header information and write results to file.
    with open(args.output, "w", encoding="utf-8") as output_file:
        output_file.write(f"#  {' '.join(sys.argv)}\n")
        output_file.write(f"# Genome-wide Average: {avgs[0]}\n")
        output_file.write(f"# Total breakpoints: {breakpoints}\n")
        if args.lrt:
            output_file.write(
                "CHROM\tPOS\tcM\tPVal\tPValCI\tPAdj\tPAdjCutoff\tSuccess\tPermutation\tDelta\tLLik\n"
            )
        else:
            output_file.write(
                "CHROM\tPOS\tcM\tPVal\tPValCI\tPAdj\tPAdjCutoff\tSuccess\tPermutation\tDelta\n"
            )

        idx = np.argsort(empp)

        for i in idx:
            if args.lrt:
                output_file.write(
                    f"{bp_ids[i][0]}\t{bp_ids[i][1]}\t{ibdfrac[i]}\t{empp[i]}\t{lower[i]},{upper[i]}\t{adjp[i]}\t{cutoff}\t{succ[i]}\t{args.nperm * args.nruns}\t{delta[i]}\t{llik_ratio[i]}\n"
                )
            else:
                output_file.write(
                    f"{bp_ids[i][0]}\t{bp_ids[i][1]}\t{ibdfrac[i]}\t{empp[i]}\t{lower[i]},{upper[i]}\t{adjp[i]}\t{cutoff}\t{succ[i]}\t{args.nperm * args.nruns}\t{delta[i]}\n"
                )

    ttotal2 = datetime.now()
    print(f"Total runtime: {ttotal2 - ttotal1}")


if __name__ == "__main__":
    main()
