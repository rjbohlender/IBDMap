
import os
import sys
import gzip
import multiprocessing as mp
import argparse as ap
from functools import partial
from bisect import bisect_left


'''
Ryan's script that converts ilash & hap-ibd output to IBDMap input

iLASH input file location: 

'''


class Indices:
    """Index container class. Provides indices for relevant fields."""
    def __init__(self, f: str):
        ## set input format
        self.id1_indx = 0
        self.id2_indx = 2
        self.chr_indx = 4
        self.str_indx = 5
        self.end_indx = 6

        if f.lower() == 'germline':
            self.cM_indx = 10
            self.unit = 11
        elif f.lower() == 'ilash':
            self.cM_indx = 9
        elif f.lower() in ['hap-ibd', 'hapibd']:
            self.cM_indx = 7
        elif f.lower() == 'rapid':
            self.chr_indx = 0
            self.id1_indx = 1
            self.id2_indx = 2
            self.cM_indx = 7
        else:
            if f.lower().split(':')[0] in ['other','others']:
                indx = in_format.lower().split(':')[1].split(';')
                self.id1_indx = int(indx[0]) - 1
                self.id2_indx = int(indx[1]) - 1
                self.chr_indx = int(indx[2]) - 1
                self.str_indx = int(indx[3]) - 1
                self.end_indx = int(indx[4]) - 1
                self.cM_indx = int(indx[5]) - 1
            else:
                sys.exit('unrecognized or incorrect format: GREMLINE/iLASH/RaPID/hap-ibd/other:id1;id2;chr;str;start bp;end bp;cM')


def find_nearest(pos_int, pos_list):
    """Return the index of the nearest position."""
    nearest_pos = min(pos_list, key=lambda x:abs(x - pos_int))
    if nearest_pos < pos_int:
        return pos_list.index(nearest_pos) + 1
    else:
        return pos_list.index(nearest_pos)


class newPOS:
    __slots__ = 'add', 'rem'
    def __init__(self, add, rem):
        self.add = add
        self.rem = rem


def is_gzipped(ipf):
    """Check if given file is gzipped."""
    with open(ipf, 'rb') as f:
        magic = f.read(2)
        return magic == b'\x1f\x8b'


def get_phenotypes(ipf):
    """Read phenotype and build possibele ID pairs"""
    uniqID = {}
    dupID = []
    pheno_file = open(ipf, 'r')
    IDnum = 0
    for line in pheno_file:
        line = line.strip()
        line = line.split('\t')
        if line[1] == 'NA':
            continue
        elif str(line[0]) in uniqID:
            dupID.append(line[0])
        else:
            uniqID[str(line[0])] = IDnum
            IDnum += 1

    print('identified '+str(len(uniqID))+' unique IDs')
    pheno_file.close()
    return uniqID, dupID

def write_temp(args, fh, pos, pair, add):
    """Open the temporary file, write our the addition or removal, then close."""
    if pos not in fh:
        fh[pos] = open(args.temp + '/{}'.format(pos), 'w')
    print('{}\t{}'.format('add' if add else 'rm', pair), file=fh[pos])

# def write_temp(args, fh, pos, pair, add):
#     if pos not in fh:
#         fh[pos] = open(args.temp + '/{}'.format(pos), 'w')
#     with fh[pos]:
#         print('{}\t{}'.format('add' if add else 'rm', pair), file=fh[pos])

def read_temp(args, pos):
    """Parse the position and update"""
    data = {'add': [], 'rm': []}
    with open(args.temp + '/{}'.format(pos), 'r') as f:
        for l in f:
            l = l.strip().split()
            data[l[0]].append(l[1])
    for k in data.keys():
        if len(data[k]) == 0:
            data[k].append('NA')
    return data



def IBDsumm(args, idx, uniqID, dupID):
    """Read and format IBD data."""
    IBDdata = {}
    IBDindex = {}
    allpos = []
    fh = {}

    ipf = args.input
    if is_gzipped(ipf):
        data = gzip.open(ipf, 'rt')
    else:
        data = open(ipf, 'rt')

    for line in data:
        line = line.strip().split()
        id1 = str(line[idx.id1_indx])
        id2 = str(line[idx.id2_indx])
        chrom = str(line[idx.chr_indx])
        start = min(int(line[idx.str_indx]), int(line[idx.end_indx]))
        end = max(int(line[idx.str_indx]), int(line[idx.end_indx]))
        cM = str(line[idx.cM_indx])
        if id1 not in uniqID or id2 not in uniqID or float(cM) < args.min or ( 'unit' in vars() and str(line[unit]) != 'cM'):
            continue
        else:
            if uniqID[id1] < uniqID[id2]:
                pair =  '{0}:{1}-{2}'.format(cM, id1, id2)
            else:
                pair =  '{0}:{1}-{2}'.format(cM, id2, id1)

            # Maintain sorted and unique allpos
            start_pos = bisect_left(allpos, start)
            if start_pos == len(allpos):
                allpos.insert(start_pos, start)
            elif allpos[start_pos] != start:
                allpos.insert(start_pos, start)

            end_pos = bisect_left(allpos, end)
            if end_pos == len(allpos):
                allpos.insert(end_pos, end)
            elif allpos[end_pos] != end:
                allpos.insert(end_pos, end)

            # Add the pair
            write_temp(args, fh, start, pair, True)
            # Remove the pair
            write_temp(args, fh, end, pair, False)

    for k, v in fh.items():
        v.close()

    print('Identified {} breakpoints on chr{}'.format(len(allpos), args.chrom))
    out = open('{}.chr{}.small.txt'.format(args.output, args.chrom), 'w')
    # Header
    out.write('chr\tpos\tadd\tdel\n')
    for pos in allpos:
        data = read_temp(args, pos)
        out.write('{}\t{}\t{}\t{}\n'.format(args.chrom, str(pos), ' '.join(data['add']), ' '.join(data['rm'])))
    out.close()
    os.system('gzip '+ args.output + '.chr' + str(args.chrom) + '.small.txt')


def main():
    """Entrypoint."""
    parser = ap.ArgumentParser()
    parser.add_argument('-i', '--input', type=str, help='Input IBD segment file.')
    parser.add_argument('-p', '--pheno', type=str, help='Input phenotype file.')
    parser.add_argument('-o', '--output', type=str, help='Prefix to the output file.')
    parser.add_argument('-f', '--format', type=str, help='Which program is the input derived from? May be GERMLINE, iLASH, RaPID, hap-ibd, or other.')
    parser.add_argument('-m', '--min', type=float, help='Minimum centimorgan length of segments to be considered.')
    parser.add_argument('-c', '--chrom', type=str, help='The current chromosome.')
    parser.add_argument('-T', '--temp', type=str, help='Path to the temporary file storage directory.')
    args = parser.parse_args()

    idx = Indices(args.format)
    print('Input: {}'.format(args.input))
    print('Phenotype file: {}'.format(args.pheno))
    print('Output: {}'.format(args.output))
    print('Input file format: {}'.format(args.format))
    print('Min output IBD length: {}'.format(args.min))

    uniqID, dupID = get_phenotypes(args.pheno)

    IBDsumm(args, idx, uniqID, dupID)


if __name__ == "__main__":
    main()