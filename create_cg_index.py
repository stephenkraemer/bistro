# TODO: update docstrings
# TODO: put on backlog: make index creation available on command line
"""
Retrieve CG positions with adjacent CHH positions from FASTA

FASTA is expected to be available in per chromosome files. Output is one
index across all chromosomes, in tab-separated BED6+1 format. Every line
represents one CG, the interval boundaries are given on the + strand.
The fourth column contains a JSON serialized list of tuples indicating
CHH positions adjacent to the CpG position and their respective strand.
The window size for including CHHs can be defined by a parameter which is
150 bp by default. Chromosomes are printed out in lexical order and are
prefixed with --prefix (empty string allowed, 'chr' by default).


Input:
Chromosome fasta filenames should be 'chr*.fa.gz' - this is the default
when the files are retrieved from the UCSC genome browser

Exemplary output:
    chrom  start  end  name            score  strand  chh_positions
    1      100    102  chr1_100_102    .      +       [('+', 100), ...]

Example:
python3 create_cg_index.py \
--cpg_content 0.05 \
--gc_content 0.5 \
--max_read_length 150 \
--out_bed /home/stephen/projects/mqc/mqc/test/data/chr11.cg.bed.gz \
--prefix chr \
/home/stephen/projects/mqc/mqc/test/data/chr11.fa.gz

# fastas must be given as absolute paths
(or several fasta files as in: )
/home/stephen/temp/chr*.fa.gz
"""

import gzip
import itertools
import numpy as np
import os
import re
import subprocess


def main():
    args = get_args()
    remove_output_file_if_present(args.out_bed_gz)

    for curr_chrom_fasta_abspath in sorted(args.chrom_fastas):
        # This sorting order means that chromosomes are printed out in lexical order
        chh_pos_strand_arr, cpg_pos_arr = initialize_pos_arrays(curr_chrom_fasta_abspath, args)
        chrom_str, chrom_seq = load_fasta(curr_chrom_fasta_abspath, args.prefix)
        print('Working on chromosome {}'.format(chrom_str))
        (chh_pos_strand_arr, cpg_pos_arr,
         ind_last_filled_row_chh_arr, ind_last_filled_cpg_arr) = get_motif_positions(
            chrom_seq, chh_pos_strand_arr, cpg_pos_arr)

        with gzip.open(args.out_bed_gz, 'at') as cg_index_fobj:
            index_smallest_chh_pos_to_consider = 0
            for curr_cg_pos_index in range(0, ind_last_filled_cpg_arr + 1):
                """
                Using iteration over index values because direct iteration over numpy array may not
                guarantee order (see nditer documentation, - didn't really check this in detail)
                """
                curr_cpg_pos = cpg_pos_arr[curr_cg_pos_index]
                adjacent_chh_pos_per_strand, index_smallest_chh_pos_to_consider = \
                    find_adjacent_chh_pos(curr_cpg_pos, chh_pos_strand_arr,
                                          index_smallest_chh_pos_to_consider,
                                          ind_last_filled_row_chh_arr,
                                          args.max_read_length)
                write_line(chrom_str, curr_cpg_pos, adjacent_chh_pos_per_strand, cg_index_fobj)


def find_motif(curr_base, up_1, up_2):
    """ Given three bases, determine presence of sequence motif

    This function is looking for CG, CHH and DDG motifs on the Watson
    strand (DDG on Watson is CHH on Crick)
    """

    if curr_base == 'C':
        if up_1 == 'G':
            return 'CG'
        elif up_1 != 'G' and up_2 != 'G':
            return 'CHH'
    elif up_1 != 'C' and up_2 == 'G':
        return 'DDG'
    else:
        return None


def write_line(chrom, curr_cpg_pos, adjacent_chh_pos_per_strand, cg_index_fobj):
    """ Write BED6+1 line for CG position with adjacent CHH positions

    The CHH positions are given as list of (strand, position) tuples,
    serialized as JSON string
    """

    name = 'CG'
    score = '.'
    strand = '+'
    csv_chh_pos_plus = ','.join(adjacent_chh_pos_per_strand['+'])
    if not csv_chh_pos_plus:
        csv_chh_pos_plus = '.'
    csv_chh_pos_minus = ','.join(adjacent_chh_pos_per_strand['-'])
    if not csv_chh_pos_minus:
        csv_chh_pos_minus = '.'
    print(chrom, curr_cpg_pos, curr_cpg_pos + 2, name, score, strand,
          csv_chh_pos_plus, csv_chh_pos_minus,
          sep='\t', end='\n', file=cg_index_fobj)


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('chrom_fastas', nargs='+')
    parser.add_argument('--cpg_content', default=0.05, type=float)
    parser.add_argument('--gc_content', default=0.50, type=float)
    parser.add_argument('--max_read_length', default=150, type=int)
    parser.add_argument('--out_bed_gz', required=True)
    parser.add_argument('--prefix', default='chr')
    return parser.parse_args()


def remove_output_file_if_present(output_path):
    """
    Chromosomes are parsed iteratively and the cg positions are appended
    to the output file, so the output file has to be deleted before we
    begin.
    """
    try:
        os.remove(output_path)
    except FileNotFoundError:
        pass


def initialize_pos_arrays(curr_chrom_fasta_abspath, args):
    """ Preallocate numpy arrays for chh and cpg positions

    This is not performance-criticial since this script works quite
    fast. Therefore, estimation of the required array size is very
    rough: we determine the number of bases in a chromosome and then
    take the fraction of CpGs and of G or C bases to determine the size
    of the CpG and CHH arrays respectively.
    """
    chrom_size = int(subprocess.check_output(
        'zcat {} | wc -c'.format(curr_chrom_fasta_abspath), shell=True))

    chh_pos_strand_arr_prealloc_size = np.round(args.gc_content * chrom_size).astype(int)
    chh_pos_strand_arr = np.empty([chh_pos_strand_arr_prealloc_size, 2], dtype=int)

    cg_pos_prealloc_size = np.round(args.cpg_content * chrom_size).astype(int)
    cpg_pos_arr = np.empty([cg_pos_prealloc_size], dtype=int)

    return chh_pos_strand_arr, cpg_pos_arr


def load_fasta(curr_chrom_fasta_abspath, prefix):
    """ Return chrom name (str) and chromosome sequence (str iterator)

    Chromosome name is prefixed with args.prefix
    """

    chrom_str = (prefix +
                 re.search('chr(.+?)\.fa\.gz$', curr_chrom_fasta_abspath).group(1))
    with gzip.open(curr_chrom_fasta_abspath, 'rt') as chr_fobj:
        _header = next(chr_fobj)
        """
        Newline removal from FASTA string is important for correct determination
        of the chromosome position
        """
        chrom_seq = iter(chr_fobj.read().replace('\n', ''))
    return chrom_str, chrom_seq


def get_motif_positions(chrom_seq, chh_pos_strand_arr, cpg_pos_arr):
    """ Fill array of CHH and CpG positions on a chromsome

    The arrays are preallocated and may contain unused positions at the end,
    which are still initialized to 0. Therefore, this function also returns
    the indices of the last array elements which were filled with motif
    positions
    """

    ind_last_filled_row_chh_arr, ind_last_filled_cpg_arr = (-1, -1)
    up_1 = next(chrom_seq)
    up_2 = next(chrom_seq)
    for chrom_pos in itertools.count(start=0, step=1):
        curr_base = up_1
        up_1 = up_2
        try:
            up_2 = next(chrom_seq)
        except StopIteration:
            break
        motif = find_motif(curr_base, up_1, up_2)
        if motif == 'CG':
            cpg_pos_arr[ind_last_filled_cpg_arr + 1] = chrom_pos
            ind_last_filled_cpg_arr += 1
        elif motif == 'CHH':
            chh_pos_strand_arr[ind_last_filled_row_chh_arr + 1, :] = (1, chrom_pos)
            ind_last_filled_row_chh_arr += 1
        elif motif == 'DDG':
            chh_pos_strand_arr[ind_last_filled_row_chh_arr + 1, :] = (0, chrom_pos + 2)
            ind_last_filled_row_chh_arr += 1

    return chh_pos_strand_arr, cpg_pos_arr, ind_last_filled_row_chh_arr, ind_last_filled_cpg_arr


def find_adjacent_chh_pos(curr_cpg_pos, chh_pos_strand_arr,
                          index_smallest_chh_pos_to_consider,
                          ind_last_filled_row_chh_arr,
                          max_read_length):
    """
    Find CHH positions surrounding a CpG position

    This function is intended for use within a loop across CpG
    positions. It finds all CHH positions within
    [CpG_position - max_read_length, CpG_position + max _read_length].
    In addition in returns the index of the smallest CHH position (in
    the CHH position array) which needs to be considered for the next
    CpG position in the CpG-loop (which is expected to be bigger than
    the current CpG position).

    The CpG position list is expected to contain all CpGs from a single
    chromosome.
    """
    adjacent_chh_pos_per_strand = {'+': [],
                                   '-': []}
    found_adjacent_chh = False
    current_chh_index = index_smallest_chh_pos_to_consider
    while True:
        strand_binary, curr_chh_pos = chh_pos_strand_arr[current_chh_index, :]
        strand_str = '+' if strand_binary == 1 else '-'
        if curr_chh_pos < curr_cpg_pos - max_read_length:
            current_chh_index += 1
            if current_chh_index > ind_last_filled_row_chh_arr:
                break
        elif curr_chh_pos < curr_cpg_pos + max_read_length:
            if not found_adjacent_chh:
                index_smallest_chh_pos_to_consider = current_chh_index
                found_adjacent_chh = True
            adjacent_chh_pos_per_strand[strand_str].append(
                (curr_chh_pos - curr_cpg_pos).astype(str))
            current_chh_index += 1
            if current_chh_index > ind_last_filled_row_chh_arr:
                break
        else:
            if not found_adjacent_chh:
                index_smallest_chh_pos_to_consider = current_chh_index
            break
    return adjacent_chh_pos_per_strand, index_smallest_chh_pos_to_consider


if __name__ == '__main__':
    main()
