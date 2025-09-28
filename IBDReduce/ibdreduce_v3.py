import glob
import argparse as ap
import sys
import multiprocessing as mp
import gzip
import zstandard as zstd
from typing import Tuple, Optional
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


def write_permutation_values(
    data: np.ndarray, idx: int, offset: int, vals: np.ndarray, nperm: int
) -> None:
    """Insert observed and permutation statistics into the shared matrix.

    Args:
        data: The matrix tracking observed and permutation statistics for all
            breakpoints.
        idx: Row index of the breakpoint being updated.
        offset: Zero-based permutation run offset.
        vals: Array containing the observed value followed by permutation
            statistics for the current run.
        nperm: Number of permutations produced per run.
    """

    if offset == 0:
        data[idx, 0] = vals[0]
        data[idx, 1 : (nperm + 1)] = vals[1:]
    else:
        start = 1 + offset * nperm
        end = 1 + (offset + 1) * nperm
        data[idx, start:end] = vals[1:]


def compute_fdr_adjusted_pvalues(empp: np.ndarray, permutation_pvalues: np.ndarray) -> np.ndarray:
    """Calculate FDR-adjusted p-values from permutation p-values.

    Args:
        empp: Empirical p-values derived from the observed statistics.
        permutation_pvalues: Ranked permutation-derived p-values for each breakpoint.

    Returns:
        Array of FDR-adjusted p-values corresponding to ``empp``.
    """

    memoize = {}
    adjp = np.zeros_like(empp, dtype=float)
    num_breakpoints = permutation_pvalues.shape[0]

    for i, p in enumerate(empp):
        if p in memoize:
            Rstar = memoize[p]
        else:
            Rstar = np.sum(permutation_pvalues <= p, axis=0)
            memoize[p] = Rstar

        rp = np.sum(empp <= p)
        pm = p * num_breakpoints
        rb = np.percentile(Rstar, 95)

        if rp - rb >= pm:
            adjp[i] = np.mean(Rstar / (Rstar + rp - pm))
        else:
            adjp[i] = np.sum(Rstar >= 1) / len(Rstar)

    return adjp


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


def parse_pheno(ifile: str, phenotype: Optional[str] = None) -> Tuple[int, int]:
    """

    :param ifile:
    :return:
    """
    data = {'0': 0, '1': 0, 'NA': 0}
    with open(ifile, 'r') as f:
        for i, l in enumerate(f):
            # we need to have specific logic for if the user passes a phenotype
            if phenotype and i == 0: # This line only occurs if there is a phenotype string and we are at the header
                header_list = l.strip().split() # We want to pass the header line and the determine where in the list the phenotype column is
                phenotype_indx = header_list.index(phenotype)
                continue
            else: # If no phenotype string then we skip the header
                phenotype_indx = 1 # If no phenotype name is passed then we are going to treat this like a 2 column file where the first column are the ids and the second column is the phenotype status
                if i == 0:
                    continue

            l = l.strip().split()
            try:
                data[l[phenotype_indx]] += 1
            except KeyError: 
                if l[phenotype_indx] in ["N/A", "-1"]: # This case will be hit if the phenotype value is not 0, 1, or NA.
                    data["NA"] += 1 # We can assume the user is meaning for these values to be excluded so just increase the count
                else: # This case account for if the user has float values as the phenotype. Ex: 1.0 or 0.0. This could error if the user gave a different string than NA or -1
                    try:
                        data[str(int(float(l[phenotype_indx])))] += 1
                    except ValueError: # This except allows the propgram to close if it encounters anything else
                        print(f"unrecognized value of {l[phenotype_indx]} at line {i+1}. Acceptable values are 1/0 for cases and controls and NA, N/A, or -1 for exclusions. Please format your data correctly and rerun the program.")
                        sys.exit(1)
    print(f"Using the provided phenotype: {phenotype}, {data['1']} cases were identified, {data['0']} were controls were identified, and {data['NA']} exclusions were identified.")
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
    parser.add_argument("--phenotype-name", default=None, help="Name of the phenotype column in the phenotype file that the user wishes to use in ibdreduce. This option only needs to be used if the user wishes to use a phenotype matrix with multiple phenotypes, otherwise te program assumes that the phenotype is in the second columns and the ids are in the first column") 

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
                        write_permutation_values(data, idx, offset, vals, args.nperm)
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
        ncase, ncontrol = parse_pheno(args.lrt, args.phenotype_name)

        cscs = ncase * (ncase - 1.) / 2.
        cscn = ncase * ncontrol

        obs_counts = obs_rates[:, 0:2]

        obs_counts[:, 0] *= cscs
        obs_counts[:, 1] *= cscn
        obs_counts = obs_counts.astype(np.int64)

        weights = np.array([cscs / (cscs + cscn), cscn / (cscs + cscn)])

        total_counts = np.sum(obs_counts * weights, 1)
        total_counts = total_counts.astype(np.int64)
        mean_rates = np.sum(obs_rates[:, 0:2] * weights, 1)

        for i in range(len(mean_rates)):
            # Null is Pois(obs_cscs, Expected_counts[cscs]) * Pois(obs_cscn, Expected_counts[cscn])
            null = stats.poisson.logpmf(obs_counts[i, 0], weights[0] * total_counts[i])
            null += stats.poisson.logpmf(obs_counts[i, 1], weights[1] * total_counts[i])
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
    # The total number of permutations is nperm * nruns, where +1 is from the observed value.
    for i in range(1, args.nperm * args.nruns + 1):
        # Add each success across breakpoints.
        succ += data[:, 0] <= data[:, i]
    empp = (succ + 1.0) / (args.nperm * args.nruns + 1)

    # Generate the evd
    # Drop the first column, which contains the observed values.
    permutation_values = data[:, 1:]
    # Rank the data in descending order, so that the smallest values are ranked highest.
    permutation_pvalues = stats.rankdata(-permutation_values, axis=1) / (args.nperm * args.nruns + 1)
    # Calculate the extreme value distribution (EVD) by taking the minimum across permutations.
    evd = np.min(permutation_pvalues, axis=0)

    if args.print_evd:
        print("EVD:\n{}".format("\t".join(str(x) for x in evd)))

    # Calculate the adjusted p-value.
    adjp = np.zeros(breakpoints)
    # FWE p-value cutoff is the 5th percentile of the EVD.
    cutoff = np.percentile(evd, 5)
    for i in range(breakpoints):
        adjp[i] = stats.percentileofscore(evd, empp[i], kind="weak") / 100.0

    # Calculate the confidence interval.
    upper = stats.chi2.ppf(0.975, 2 * succ + 2) / 2.0
    lower = stats.chi2.ppf(0.025, 2 * succ) / 2.0
    lower[np.isnan(lower)] = 0

    upper /= permutation_pvalues.shape[1]
    lower /= permutation_pvalues.shape[1]

    if args.fdr:
        adjp = compute_fdr_adjusted_pvalues(empp, permutation_pvalues)

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
