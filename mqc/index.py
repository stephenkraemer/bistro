import gzip
import itertools
import re
from collections import OrderedDict
from os.path import dirname, isabs, join
import numpy as np
import os
import os.path as op
import subprocess
from typing import List, Iterable, Dict, Any
from joblib import Parallel, delayed
from mqc.utils import open_gzip_or_plain_file


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
    # TODO: support additional index fields
    def __init__(self, index_line: str):
        fields = index_line.rstrip().split('\t')
        # TODO: why?
        self.chrom = fields[0].replace('chrom', '')
        self.start = int(fields[1])
        self.end = int(fields[2])
        self.motif = fields[3]
        self.score = fields[4]
        self.strand = fields[5]  # '+' or '-'
        self.watson_base = 'C' if (self.strand == '+') else 'G'

    def __str__(self):
        return '{}:{}-{}, motif: {} (on strand {})'.format(
            self.chrom, self.start, self.end,
            self.original_motif, self.strand)


def start_parallel_index_generation(genome_fasta: str, index_output_dir: str,
                                    motifs: List[str], annotations: OrderedDict,
                                    cores: int):
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

    chroms = find_chroms(genome_fasta)
    genome_name = re.sub(r'.fa(.gz)?$', '', op.basename(genome_fasta))
    os.makedirs(index_output_dir, exist_ok=True)

    motifs_str = '-'.join(motifs)
    # NOTE: the hardcoded template is used by PileupRun._single_run_over_index_file
    #       to find which chromosome it is working on
    output_file_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
    output_files = [output_file_template.format(chrom=c) for c in chroms]

    Parallel(n_jobs=cores)(delayed(fasta_to_index)
                           (genome_fasta,
                            chrom=curr_chrom,
                            output_fp=curr_output_fp,
                            motifs=motifs,
                            annotations=annotations)
                           for curr_chrom, curr_output_fp in zip(chroms, output_files))

    print('Completed index generation')


def find_chroms(genome_fasta):
    chroms = []
    with open_gzip_or_plain_file(genome_fasta) as fobj:
        for line in fobj:
            if line.startswith('>'):
                chrom = line.split()[0].replace('>', '')
                chroms.append(chrom)
    return chroms


def fasta_to_index(genome_fasta, chrom:str, output_fp: str,
                   motifs: List[str], annotations: Dict[str, Any]):

    print(f"Working on {chrom}")
    chrom_seq = read_fasta(genome_fasta, chrom)
    bed_lines = bed_lines_generator(chrom_seq, chrom,
                                    motifs=motifs, annotations=annotations)

    with gzip.open(output_fp, 'wt') as out_fobj:
        base_header = ['#chrom', 'start', 'end', 'motif', 'score', 'strand']
        anno_cols = [k for k, v in annotations.items() if v]
        header_line = '\t'.join(base_header + anno_cols) + '\n'

        out_fobj.write(header_line)
        for curr_line in bed_lines:
            out_fobj.write('\t'.join(curr_line) + '\n')


def read_fasta(genome_fasta, chrom) -> str:

    def read_lines(fobj):
        lines = []
        for curr_line in fobj:
            if curr_line.startswith('>'):
                break
            lines.append(curr_line.strip().upper())
        chrom_seq = ''.join(lines)
        return chrom_seq

    chrom_seq = ''
    with open_gzip_or_plain_file(genome_fasta) as fobj:
        for curr_line in fobj:
            if curr_line.startswith('>'):
                curr_chrom = curr_line.split()[0].replace('>', '')
                if curr_chrom == chrom:
                    chrom_seq = read_lines(fobj)
                    break

    if not chrom_seq:
        raise ValueError(f"Could not find chromosome {chrom} in\n{genome_fasta}")

    return chrom_seq



def bed_lines_generator(fasta_seq: str, chr_name: str, motifs: List[str],
                        annotations: Dict[str, Any]) -> Iterable[List[str]]:

    for i in range(0, len(fasta_seq)):

        curr_base = fasta_seq[i]
        if curr_base in ['C', 'G']:

            # Note that len(triplet_seq) < 3 at the ends of fasta_seq
            if curr_base == 'C':
                triplet_seq = fasta_seq[i:i+3]
                strand = '+'
            else:
                start = 0 if (i - 2 < 0) else (i - 2)
                triplet_seq = reverse_complement_seq(fasta_seq[start:i+1])
                strand = '-'

            motif = classify_motif(triplet_seq)
            # TODO: document that this discard CN motifs, which are not choosable currently (?)
            if motif not in motifs:  # this also discards motif='.' and motif='CN'
                # Motifs at the boundaries of fasta_seq may not be callable
                continue

            bed_line = [chr_name, str(i), str(i+1), motif, '.', strand]

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

    #TODO: case seq[2] == N

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



