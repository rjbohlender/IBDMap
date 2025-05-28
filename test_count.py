import argparse as ap
import gzip
import sys


def main():
    parser = ap.ArgumentParser()
    parser.add_argument('phenotypes')
    parser.add_argument('ibd')
    args = parser.parse_args()

    pheno = {}
    with open(args.phenotypes, 'r') as f:
        for l in f:
            l = l.strip().split()
            if l[1] == 'NA':
                continue
            pheno[l[0]] = int(l[1])

    cscs = 0
    cscn = 0
    cncn = 0
    seen = set()
    duplicates = 0
    autozygous = 0
    with gzip.open(args.ibd, 'rt') as f:
        for l in f:
            if l.startswith('chr'):
                continue
            chrom, pos, segs, pairs, add, dele = l.strip().split('\t')
            for p in add.split():
                try:
                    a, b = p.split(':')[1].split('-')
                except:
                    print(p)
                    sys.exit()
                if a == b:
                    autozygous += 1
                    continue
                if f'{a},{b}' in seen or f'{b},{a}' in seen:
                    duplicates += 1
                    continue
                else:
                    if a <= b:
                        seen.add(f'{a},{b}')
                    else:
                        seen.add(f'{b},{a}')

                if a in pheno and b in pheno:
                    x = pheno[a]
                    y = pheno[b]
                    if x == 1 and y == 1:
                        cscs += 1
                    elif x == 1 and y == 0:
                        cscn += 1
                    elif x == 0 and y == 1:
                        cscn += 1
                    else:
                        cncn += 1
            break
    print(f'cscs: {cscs} cscn: {cscn} cncn: {cncn}')
    print(f'duplicates: {duplicates}')
    print(f'autozygous: {autozygous}')
    print(f'total: {cscs + cscn + cncn}')


if __name__ == "__main__":
    main()