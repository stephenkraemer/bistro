import gzip
import os
import os.path as op
from abc import ABCMeta
from pathlib import Path
from typing import List, Optional

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup
from mqc.flag_and_index_values import methylation_status_flags as mflags


class BedWriter(Visitor):

    def __init__(self, config, chrom: str):

        # TODO: avoid hard coding
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

    def process(self, motif_pileup: MotifPileup):
        motif = motif_pileup.idx_pos.motif
        line_no_endl = '\t'.join([
            motif_pileup.idx_pos.chrom,
            str(motif_pileup.idx_pos.start),
            str(motif_pileup.idx_pos.end),
            motif,
            '.',
            motif_pileup.idx_pos.strand,
            f"{motif_pileup.beta_value:.6f}",
            str(motif_pileup.n_meth),
            str(motif_pileup.n_total),
        ])
        self.meth_calls_fobj_dict[motif].write(line_no_endl + '\n')


class McallWriter(Visitor, metaclass=ABCMeta):
    def __init__(self, calls_by_chrom_motif_fp: str,
                 motifs: List[str], header_no_newline: Optional[str],
                 chrom: str):
        """Write MotifPileup info to mcall file

        Args:
            calls_by_chrom_motif_fp: must contain
                [motif] and [chrom] fields
            motifs: list of motifs covered in the run
            header_no_newline: optional header line for all output files
            chrom: chrom ID
        """
        self.calls_by_chrom_motif_fp = str(calls_by_chrom_motif_fp)
        self.motifs = motifs
        self.meth_calls_fobj_dict = {}
        self.header_no_newline = header_no_newline
        self.chrom = chrom

    def setup(self):
        """Open file objects for each motif covered by the run, add header"""

        self._open_mcall_files()
        if self.header_no_newline:
            self._write_header_line()
        return self

    def _open_mcall_files(self):
        for motif in self.motifs:
            out_fp = (self.calls_by_chrom_motif_fp
                      .replace('[motif]', motif)
                      .replace('[chrom]', self.chrom)
                      )
            # Note that there may be one output folder per motif,
            # so parent creation must come after pattern expansion
            Path(out_fp).parent.mkdir(exist_ok=True, mode=0o770)
            self.meth_calls_fobj_dict[motif] = gzip.open(out_fp, 'wt')

    def _write_header_line(self):
        for fobj in self.meth_calls_fobj_dict.values():
            fobj.write(self.header_no_newline + '\n')


class BismarkWriter(McallWriter):
    def __init__(self, calls_by_chrom_motif_fp: str, motifs: List[str],
                 chrom: str):
        """Generate output as specified by Bismark methylation extractor

        Args:
            calls_by_chrom_motif_fp: must contain
                [motif] and [chrom] fields
            motifs: list of motifs covered in the run
            chrom: chrom ID
        """
        super().__init__(calls_by_chrom_motif_fp=calls_by_chrom_motif_fp,
                         motifs=motifs,
                         header_no_newline='',
                         chrom=chrom)
        self.symbol_table = {
            'CG': {
                mflags.is_methylated: ('+', 'Z'),
                mflags.is_unmethylated: ('-', 'z')
            },
            'CHG': {
                mflags.is_methylated: ('+', 'X'),
                mflags.is_unmethylated: ('-', 'x')
            },
            'CHH': {
                mflags.is_methylated: ('+', 'H'),
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
