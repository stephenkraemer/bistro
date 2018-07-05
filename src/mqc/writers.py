"""Write methylation stats to compressed BED files

Compressed as gz, currently non-optional.

Supports different formats
- Bed
- StratifiedBed (global methylation stats and stratified by BS-Seq strand
  and mate)
- Bismark methylation extractor output

Planned
- VCF (in combination with SNP calling module)
- biscuit like epiread format
"""

import gzip
import os
import os.path as op
from abc import ABCMeta
from itertools import chain
from pathlib import Path
from typing import List, Dict, Any, IO

import numpy as np

from mqc.flag_and_index_values import (
    methylation_status_flags as mflags,
    meth_status_indices as mstat_ids,
    strat_call_indices as scall_ids,
)
from mqc.pileup.pileup import MotifPileup
from mqc.visitors import Visitor


class BedWriter(Visitor):
    """Standard BED6 + [beta_value n_meth n_total] format"""
    # TODO: subclass McallWriterABC

    def __init__(self, config, chrom: str) -> None:
        # TODO-important: avoid hard coding
        meth_calls_path_template = (
                config['paths']['meth_calls_basepath']
                + f"_{config['sample']['name']}_{{motif}}_{chrom}.bed.gz")

        os.makedirs(op.dirname(meth_calls_path_template), exist_ok=True, mode=0o770)

        self.meth_calls_fobj_dict = {
            motif: gzip.open(meth_calls_path_template.format(motif=motif), 'wt')
            for motif in config['run']['motifs']}

        header_no_endl = '\t'.join(['#chrom', 'start', 'end',
                                    'motif', 'score', 'strand',
                                    'beta_value', 'n_meth', 'n_total'])
        for fobj in self.meth_calls_fobj_dict.values():
            fobj.write(header_no_endl + '\n')

    def process(self, motif_pileup: MotifPileup) -> None:
        motif = motif_pileup.idx_pos.motif
        line_no_endl = '\t'.join([
            motif_pileup.idx_pos.chrom,
            str(motif_pileup.idx_pos.start),
            str(motif_pileup.idx_pos.end),
            motif,
            '.',
            motif_pileup.idx_pos.strand,
            f"{motif_pileup.beta_value:.8f}",
            str(motif_pileup.n_meth),
            str(motif_pileup.n_total),
        ])
        self.meth_calls_fobj_dict[motif].write(line_no_endl + '\n')


class McallWriterABC(Visitor, metaclass=ABCMeta):
    """ABC for classes writing meth. calling statistics to text files

    Args:
        calls_by_chrom_motif_fp: must contain [motif] and [chrom] fields
        motifs: list of motifs covered in the run
        header_no_newline: optional header line for all output files.
            Starts with "#"
        chrom: chrom ID
    """

    def __init__(self, calls_by_chrom_motif_fp: str,
                 motifs: List[str], header_no_newline: str,
                 chrom: str) -> None:
        self.calls_by_chrom_motif_fp = str(calls_by_chrom_motif_fp)
        self.motifs = motifs
        self.meth_calls_fobj_dict: Dict[str, IO[Any]] = {}
        if header_no_newline == '':
            header_line_is_ok = True
        elif ('\t' in header_no_newline
              and header_no_newline.startswith('#')
              and not header_no_newline.endswith('\n')):
            header_line_is_ok = True
        else:
            header_line_is_ok = False
        assert header_line_is_ok, \
            f'Header line {header_no_newline} does not have correct format'
        self.header_no_newline = header_no_newline
        self.chrom = chrom

    def setup(self) -> 'McallWriterABC':
        """Open file objects for each motif covered by the run, add header"""

        self._open_mcall_files()
        if self.header_no_newline:
            self._write_header_line()
        return self

    def _open_mcall_files(self) -> None:
        for motif in self.motifs:
            out_fp = (self.calls_by_chrom_motif_fp
                      .replace('[motif]', motif)
                      .replace('[chrom]', self.chrom)
                      )
            # Note that there may be one output folder per motif,
            # so parent creation must come after pattern expansion
            Path(out_fp).parent.mkdir(parents=True, exist_ok=True, mode=0o770)
            self.meth_calls_fobj_dict[motif] = gzip.open(out_fp, 'wt')

    def _write_header_line(self) -> None:
        for fobj in self.meth_calls_fobj_dict.values():
            fobj.write(self.header_no_newline + '\n')


class BismarkWriter(McallWriterABC):
    """Generate output as specified by Bismark methylation extractor

    Args:
        calls_by_chrom_motif_fp: must contain
            [motif] and [chrom] fields
        motifs: list of motifs covered in the run
        chrom: chrom ID
    """

    def __init__(self, calls_by_chrom_motif_fp: str, motifs: List[str],
                 chrom: str) -> None:
        super().__init__(calls_by_chrom_motif_fp=calls_by_chrom_motif_fp,
                         motifs=motifs,
                         header_no_newline='',
                         chrom=chrom)
        self.symbol_table = {
            'CG':  {
                mflags.is_methylated:   ('+', 'Z'),
                mflags.is_unmethylated: ('-', 'z')
            },
            'CHG': {
                mflags.is_methylated:   ('+', 'X'),
                mflags.is_unmethylated: ('-', 'x')
            },
            'CHH': {
                mflags.is_methylated:   ('+', 'H'),
                mflags.is_unmethylated: ('-', 'h')
            },
        }

    def process(self, motif_pileup: MotifPileup) -> None:
        """Extract Bismark output lines from MotifPileup

        Ignores reads with discard-flags (qc fail, trimming or overlap)

        Only considers reads with 'methylated' or 'unmethylated' status,
        not SNPs, NA, Ref...
        """
        pos_str = str(motif_pileup.idx_pos.start)
        curr_fout = self.meth_calls_fobj_dict[motif_pileup.idx_pos.motif]
        curr_symbol_table = self.symbol_table[motif_pileup.idx_pos.motif]
        for read in motif_pileup.reads:
            if (read.qc_fail_flag
                    or read.trimm_flag
                    or read.overlap_flag):
                continue
            try:
                strand_symbol, meth_symbol = curr_symbol_table[
                    read.meth_status_flag]
            except KeyError:
                # read has no meth status (SNP, Ref...)
                continue
            curr_fout.write('\t'.join(
                (read.alignment.query_name, strand_symbol,
                 self.chrom, pos_str, meth_symbol)) + '\n')


class StratifiedBedWriter(McallWriterABC):
    """Stratified methylation calls (BED6 + methylation stats)

    The first 9 columns are identical to the columns of the standard
    BED output format (BED6 + beta_value n_meth n_total).
    The following columns indicate the meth stats per
    BS-Seq strand, and then per mate. Floats are written as %.8f

    Subclasses McallWriterABC and implements the process method. File
    handling behavior is inherited.

    Args:
        calls_by_chrom_motif_fp: must contain [motif] and [chrom] fields
        motifs: list of motifs covered in the run
        chrom: chrom ID
    """

    header_no_newline = '\t'.join(
        ['#chrom', 'start', 'end', 'motif', 'score', 'strand',
         'beta_value', 'n_meth', 'n_total',
         'c_bc_beta_value', 'c_bc_n_meth', 'c_bc_n_total',
         'c_bc_rv_beta_value', 'c_bc_rv_n_meth', 'c_bc_rv_n_total',
         'w_bc_beta_value', 'w_bc_n_meth', 'w_bc_n_total',
         'w_bc_rv_beta_value', 'w_bc_rv_n_meth', 'w_bc_rv_n_total',
         'mate1_beta_value', 'mate1_n_meth', 'mate1_n_total',
         'mate2_beta_value', 'mate2_n_meth', 'mate2_n_total',
         ])

    def __init__(self, calls_by_chrom_motif_fp: str, motifs: List[str], chrom: str) -> None:
        super().__init__(calls_by_chrom_motif_fp=calls_by_chrom_motif_fp,
                         motifs=motifs,
                         header_no_newline=self.header_no_newline,
                         chrom=chrom)

    def process(self, motif_pileup: MotifPileup) -> None:
        """Create stratified meth. calls BED line

        Use MotifPileup attributes filled by StratifiedMethCaller.
        """

        strat_calls_slice = slice(0, scall_ids.all)

        line = '\t'.join(chain(
            [motif_pileup.idx_pos.chrom, str(motif_pileup.idx_pos.start),
             str(motif_pileup.idx_pos.end), motif_pileup.idx_pos.motif,
             '.', motif_pileup.idx_pos.strand,
             f"{motif_pileup.beta_value:.8f}", str(motif_pileup.n_meth), str(motif_pileup.n_total),
             ],
            *zip(np.char.mod('%.8f', motif_pileup.strat_beta_arr[strat_calls_slice]),
                 motif_pileup.meth_counts_arr[strat_calls_slice, mstat_ids.n_meth].astype(str),
                 motif_pileup.meth_counts_arr[strat_calls_slice, mstat_ids.n_total].astype(str)),
            ('\n',)
        ))
        # astype(str) is sufficient because meth_counts_array dtype is int

        self.meth_calls_fobj_dict[motif_pileup.idx_pos.motif].write(line)
