import gzip
import itertools
import re
from collections import OrderedDict
from os.path import dirname, isabs, join
import numpy as np
import os
import subprocess
from typing import List, Iterable, Dict, Any
from joblib import Parallel, delayed
import mqc.utils


def make_index():
    pass

class IndexFile:
    def __init__(self, bed_abspath):
        self.index_fobj = gzip.open(bed_abspath, 'rt')
        next(self.index_fobj)  # skip header

    def __next__(self):
        next_index_line = next(self.index_fobj)
        return IndexPosition(next_index_line)

    def __iter__(self):
        return self


class IndexPosition:
    def __init__(self, index_line: str):
        fields = index_line.rstrip().split('\t')
        self.chrom = fields[0].replace('chr', '')
        self.start = int(fields[1])
        self.end = int(fields[2])
        # TODO: use name for motif? -> name as motif on the respective strand?
        self.name = fields[3]
        self.score = fields[4]
        self.strand = fields[5]  # '+' or '-'
        # TODO: backlog -> only process these fields if required (perhaps implement through properties?)

        def parse_control_pos_str(s):
            if s == '.':
                return []
            else:
                return [int(curr_pos_str) for curr_pos_str in s.split(',')]

        self.watson_conv_control_cyts = parse_control_pos_str(fields[6])
        self.crick_conv_control_cyts = parse_control_pos_str(fields[7])

        self.original_motif = self.name
        # TODO: check this!
        self.watson_motif = (self.name
                             if self.strand == '+'
                             else reverse_complement_seq(self.name))

    def __str__(self):
        return '{}:{}-{}, motif: {} (on strand {})'.format(
            self.chrom, self.start, self.end,
            self.original_motif, self.strand)


def reverse_complement_seq(seq):
    return seq.translate(str.maketrans('CGHD', 'GCDH')).reverse()


def fasta_to_index(fasta_fp: str, output_fp: str,
                   chr_name, motifs: List[str], annotations: Dict[str, Any]):

    fasta_seq = read_fasta(fasta_fp)
    bed_lines = bed_lines_generator(fasta_seq, chr_name,
                                    motifs=motifs, annotations=annotations)
    with gzip.open(output_fp, 'wt') as out_fobj:

        base_header = ['chr', 'start', 'end', 'motif', 'score', 'strand']
        anno_cols = [k for k, v in annotations.items() if v]
        header_line = '\t'.join(base_header + anno_cols) + '\n'

        out_fobj.write(header_line)
        for curr_line in bed_lines:
            out_fobj.write('\t'.join(curr_line) + '\n')


def read_fasta(fasta_fp) -> str:
    """ Return chrom sequence from fasta file (plain or gzip) as string"""

    with mqc.utils.open_gzip_or_plain_file(fasta_fp) as chr_fobj:
        _header = next(chr_fobj)
        chrom_seq = (chr_fobj.read()
                     .replace('\n', '')
                     .upper())
    return chrom_seq

def bed_lines_generator(fasta_seq: str, chr_name: str, motifs: List[str],
                        annotations: Dict[str, Any]) -> Iterable[List[str]]:

    for i in range(0, len(fasta_seq)):
        curr_base = fasta_seq[i]
        if curr_base in ['C', 'G']:

            bed_line = [chr_name, str(i), str(i+1)]

            # Note that len(triplet_seq) < 3 at the ends of fasta_seq
            if curr_base == 'C':
                triplet_seq = fasta_seq[i:i+3]
                strand = '+'
            else:
                start = 0 if (i - 2 < 0) else (i - 2)
                triplet_seq = reverse_complement_seq(fasta_seq[start:i+1])
                strand = '-'

            motif = classify_motif(triplet_seq)
            if motif == '.':
                # Motifs at the boundaries of fasta_seq may not be callable
                continue

            bed_line += [motif, '.', strand]

            for k, v in annotations.items():
                if v:
                    field = ANNOTATORS[k](fasta_seq, i, strand, triplet_seq, v)
                    bed_line.append(field)

            yield bed_line


def classify_motif(seq):

    if len(seq) == 1:
        return '.'

    if seq[1] == 'G':
        return 'CG'
    if seq[1] == 'N':
        return 'CN'

    if len(seq) < 3:
        return '.'

    if seq[2] == 'G':
        return 'CHG'

    return 'CHH'

def reverse_complement_seq(seq):
    return seq.translate(str.maketrans('CGAT',
                                       'GCTA'))[-1::-1]


def add_triplet_seq(fasta_seq, i, strand, triplet_seq, value):
    if len(triplet_seq) == 1:
        triplet_seq += "$$"
    elif len(triplet_seq) == 2:
        triplet_seq += "$"
    return triplet_seq

def add_seq_context(fasta_seq, i, strand, triplet_seq, value):

    n_bp_context = value

    try:
        start = (i - n_bp_context)
        end = (i + n_bp_context + 1)
        seq_context = fasta_seq[start:end]
    except IndexError:  # we are at the start of fasta_seq
        seq_context =  '.'

    if len(seq_context) < 2 * n_bp_context + 1:  # we are at the end
        seq_context =  '.'

    if strand == '-':
        seq_context = reverse_complement_seq(seq_context)

    return seq_context

ANNOTATORS = {'triplet_seq': add_triplet_seq,
              'seq_context': add_seq_context}



def parallel_index_generation(fasta_path_template: str, output_path_template: str,
                              motifs: List[str], annotations: OrderedDict,
                              cores: int, chr_prefix='chr',
                              fasta_to_index_fn=fasta_to_index):
    """Determine file paths and start parallel index generation

    Parameters
    ----------
    chr_prefix:
        string prepended to chromosome index
    fasta_path_template:
        template path specifying where to find the fasta files. Must have
        one field, named '{chr}'. This field will be expanded to match the different chromosome
        indices of the individual fasta files. Must be absolute path.
    output_path_template:
        template path specifying where to save the generated index files. Must
        have one field, named '{chr}'. This field will be replaced with the chromosome
        index of the corresponding fasta file
    motifs:
        the motifs to be considered for the index generation (CG, CHH, CHG/CWG)
    annotations:
        strings specifying the desired annotation columns
    cores:
        index generation can be run in parallel with *cores* processes
    fasta_to_index_fn:
        dependency injection used for unit test
    """

    # Note: Need to pass fasta_to_index_fn to this function
    # to facilitate unit test

    if not (isabs(fasta_path_template) and isabs(output_path_template)):
        raise ValueError('Template paths must be specified as absolute paths')

    fasta_dir = dirname(fasta_path_template)
    fasta_path_regex = fasta_path_template.replace(
        '{chr}', r'(\d{1,2}|[a-zA-Z]{1,2})')
    matches_in_fasta_dir = [re.match(fasta_path_regex, join(fasta_dir, fp))
                            for fp in os.listdir(fasta_dir)]
    job_params = [
        {'fasta_fp': m.group(0),
         'output_fp': output_path_template.replace('{chr}', m.group(1)),
         'chr_name': chr_prefix + m.group(1)}
        for m in matches_in_fasta_dir if m]

    if not job_params:
        raise(ValueError('Could not find any FASTA files!'))

    # Saving return values from fasta_to_index_fn helps with patching
    # the function in unit test. Would not be necessary otherwise.
    out_fps = Parallel(n_jobs=cores)(delayed(fasta_to_index_fn)
                                     (fasta_fp=curr_params['fasta_fp'],
                                      output_fp=curr_params['output_fp'],
                                      chr_name=curr_params['chr_name'],
                                      motifs=motifs,
                                      annotations=annotations)
                                     for curr_params in job_params)

    print('Completed index generation')

    return out_fps
