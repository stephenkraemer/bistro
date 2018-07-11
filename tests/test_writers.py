"""Tests for all writer classes

BedWriter, BismarkWriter, McallWriterABC, StratifiedBedWriter

Acceptance tests are part of test_mcall_cmd.py
"""

import gzip
import re

import os.path as op
import numpy as np


from typing import Any, cast
from collections import namedtuple
from io import TextIOWrapper
from pathlib import Path
from unittest.mock import MagicMock, call
from textwrap import dedent
from dataclasses import dataclass

from mqc.pileup.pileup import MotifPileup
from mqc.writers import BedWriter, BismarkWriter, McallWriterABC, StratifiedBedWriter
from mqc.utils import get_resource_abspath
from mqc.flag_and_index_values import (
    methylation_status_flags as mflags,
    meth_status_indices as mstat_ids,
    strat_call_indices as scall_ids,
)


TESTS_DIR = op.dirname(__file__)
DEFAULT_CONFIG_FILE = get_resource_abspath('config.default.toml')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'


def test_bed_writer_computes_bed_lines_and_writes_to_one_file_per_motif(mocker):
    """
    input: MotifPileup
    - compute bed line
    side_effect: write BED line to the mcalls file for the given motif
    """

    index_pos_stub_cg = namedtuple('IndexPositionStub', 'chrom start motif strand end')(
            chrom='1',
            start=10,
            end=11,
            motif='CG',
            strand='+',
    )

    motif_pileup_stub_cg = namedtuple('MotifPileupStub', 'idx_pos beta_value n_meth n_total')(
            idx_pos=index_pos_stub_cg,
            beta_value=0.6,
            n_meth=3,
            n_total=5
    )

    index_pos_stub_chg = namedtuple('IndexPositionStub', 'chrom start motif strand end')(
            chrom='1',
            start=10,
            end=11,
            motif='CHG',
            strand='+',
    )

    motif_pileup_stub_chg = namedtuple('MotifPileupStub', 'idx_pos beta_value n_meth n_total')(
            idx_pos=index_pos_stub_chg,
            beta_value=1,
            n_meth=5,
            n_total=5
    )
    config_stub = {'paths':  {'meth_calls_basepath': 'does/not/matter'},
                   'sample': {'name': 'hsc_1'},
                   'run':    {'motifs': 'CG CHG CHH'.split()}}

    file_mocks = [MagicMock(), MagicMock(), MagicMock()]
    mock_open = MagicMock(side_effect=file_mocks)
    mocker.patch("mqc.writers.gzip.open", mock_open)

    writer = BedWriter(config_stub, chrom='1')
    # noinspection PyTypeChecker
    writer.process(motif_pileup_stub_cg)
    # noinspection PyTypeChecker
    writer.process(motif_pileup_stub_chg)

    expected_line_cg = f"1 10 11 CG . + {3/5:.8f} 3 5\n".replace(' ', '\t')
    expected_line_chg = f"1 10 11 CHG . + {1:.8f} 5 5\n".replace(' ', '\t')

    expected_header_no_endl = '\t'.join(['#chrom', 'start', 'end',
                                         'motif', 'score', 'strand',
                                         'beta_value', 'n_meth', 'n_total'])
    order_of_file_openening = [re.search(r'(CG|CHG|CHH)', x[0][0]).group(1)
                               for x in mock_open.call_args_list]
    file_mock_motif_mapping = {x: file_mocks[i]
                               for i, x in enumerate(order_of_file_openening)}
    assert file_mock_motif_mapping['CG'].write.call_args_list == [call(expected_header_no_endl + '\n'),
                                                                  call(expected_line_cg)]
    assert file_mock_motif_mapping['CHG'].write.call_args_list == [call(expected_header_no_endl + '\n'),
                                                                   call(expected_line_chg)]


def test_mcall_writer_abc_setup(tmpdir):
    """Test setup code provided through ABC

       Open mcall files, write header and provide file objects through
       dict via setup() method
    """
    class WriterStub(McallWriterABC):
        def process(self, motif_pileup):
            pass

    meth_calls_template_fp = Path(tmpdir) / 'subdir/test_[motif]_chr-[chrom].bed.gz'

    header_line = '#chrom\tstart\tend'  # must start with '#'
    writer_stub = WriterStub(
            calls_by_chrom_motif_fp=meth_calls_template_fp.as_posix(),
            motifs='CG CHH'.split(),
            header_no_newline=header_line,
            chrom='1')
    writer_stub.setup()

    assert meth_calls_template_fp.parent.exists()
    assert isinstance(writer_stub.meth_calls_fobj_dict['CG'],
                      TextIOWrapper)
    assert isinstance(writer_stub.meth_calls_fobj_dict['CHH'],
                      TextIOWrapper)

    # delete writer to close file objects
    del writer_stub

    with gzip.open(meth_calls_template_fp.with_name('test_CG_chr-1.bed.gz')
                   .as_posix(), 'rt') as fin:
        assert fin.read() == header_line + '\n'
    with gzip.open(meth_calls_template_fp.with_name('test_CHH_chr-1.bed.gz')
                   .as_posix(), 'rt') as fin:
        assert fin.read() == header_line + '\n'


def test_bismark_writer_writes_bismark_format_to_motif_files(tmpdir):
    """BismarkWriter has to correctly handle QC issues itself

    Other writers use precalculated and precorrected meth. calling stats

    Test that BismarkWriter discards itself:
    - qc_fail_flag
    - phred_fail_flag
    - overlap_flag

    Test with different motifs from one PileupGenerator
    """

    MotifPileupStub = namedtuple('MotifPileupStub', 'idx_pos reads')
    IndexPositionStub = namedtuple('IndexPositionStub', 'chrom start motif')
    AlignedSegmentStub = namedtuple('AlignedSegmentStub', 'query_name')

    @dataclass(frozen=True)
    class BSSeqPileupReadStub:  # pylint: disable=too-few-public-methods
        meth_status_flag: int
        alignment: AlignedSegmentStub
        overlap_flag: int = 0
        qc_fail_flag: int = 0
        phred_fail_flag: int = 0
        trimm_flag: int = 0

    # Test that undefined meth. calls are discarded
    motif_pileup_cg_1 = MotifPileupStub(
            idx_pos=IndexPositionStub(chrom='1', start=10, motif='CG'),
            reads=[
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    alignment=AlignedSegmentStub(query_name='readname1')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_na,
                                    alignment=AlignedSegmentStub(query_name='readname_na')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_ref,
                                    alignment=AlignedSegmentStub(query_name='readname_ref')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_snp,
                                    alignment=AlignedSegmentStub(query_name='readname_snp'))
            ]
    )

    # Test that QC issue and overlap cases are discarded
    motif_pileup_cg_2 = MotifPileupStub(
            idx_pos=IndexPositionStub(chrom='1', start=11, motif='CG'),
            reads=[
                BSSeqPileupReadStub(meth_status_flag=mflags.is_unmethylated,
                                    alignment=AlignedSegmentStub(query_name='readname2')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    alignment=AlignedSegmentStub(query_name='readname3')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    overlap_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_overlap')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    trimm_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_trimm')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    qc_fail_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_qcfail')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    qc_fail_flag=1, trimm_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_multiple1')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    phred_fail_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_multiple2')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    phred_fail_flag=1, trimm_flag=1,
                                    alignment=AlignedSegmentStub(query_name='readname_multiple3')),
            ]
    )

    motif_pileup_chh = MotifPileupStub(
            idx_pos=IndexPositionStub(chrom='1', start=15, motif='CHH'),
            reads=[
                BSSeqPileupReadStub(meth_status_flag=mflags.is_unmethylated,
                                    alignment=AlignedSegmentStub(query_name='readname1')),
                BSSeqPileupReadStub(meth_status_flag=mflags.is_methylated,
                                    alignment=AlignedSegmentStub(query_name='readname2'))
            ]
    )

    meth_calls_template_fp = str(tmpdir / 'subdir/mcalls_hsc1_[motif]_chr-[chrom].bed.gz')
    motifs = 'CG CHH CHG'.split()
    bismark_writer = BismarkWriter(
            calls_by_chrom_motif_fp=meth_calls_template_fp,
            motifs=motifs,
            chrom='1')

    bismark_writer.setup()
    # noinspection PyTypeChecker
    bismark_writer.process(motif_pileup_cg_1)
    # noinspection PyTypeChecker
    bismark_writer.process(motif_pileup_cg_2)
    # noinspection PyTypeChecker
    bismark_writer.process(motif_pileup_chh)

    # delete writer object to close mcall file objects
    del bismark_writer

    expected_cg_mcall_file_content = dedent("""\
    readname1 + 1 10 Z
    readname2 - 1 11 z
    readname3 + 1 11 Z
    """.replace(' ', '\t'))

    expected_chh_mcall_file_content = dedent("""\
    readname1 - 1 15 h
    readname2 + 1 15 H
    """.replace(' ', '\t'))

    cg_call_file = (meth_calls_template_fp
                    .replace('[chrom]', '1')
                    .replace('[motif]', 'CG'))
    with gzip.open(cg_call_file, 'rt') as fin:
        res = fin.read()
    assert res == expected_cg_mcall_file_content

    chh_call_file = cg_call_file.replace('CG', 'CHH')
    with gzip.open(chh_call_file, 'rt') as fin:
        res = fin.read()
    assert res == expected_chh_mcall_file_content


StratifiedBedIndexPositionStub = namedtuple('StratifiedBedIndexPositionStub', 'chrom start motif strand end')


class StratifiedBedMotifPileupStub:
    def __init__(self, meth_counts_arr: np.ndarray, idx_pos: StratifiedBedIndexPositionStub) -> None:
        self.idx_pos = idx_pos
        self.meth_counts_arr = meth_counts_arr
        with np.errstate(divide='raise', invalid='ignore'):
            self.strat_beta_arr = (self.meth_counts_arr[:, mstat_ids.n_meth] /
                                   self.meth_counts_arr[:, mstat_ids.n_total])
        self.n_meth = self.meth_counts_arr[scall_ids.all, mstat_ids.n_meth]
        self.n_total = self.meth_counts_arr[scall_ids.all, mstat_ids.n_total]
        self.beta_value = self.strat_beta_arr[scall_ids.all]


def test_stratified_bed_writer(tmpdir: Any) -> None:
    """Test StratifiedBedWriter: setup and process multiple pileups"""
    index_position_stub_cg1 = StratifiedBedIndexPositionStub(chrom='11', start=10, end=11,
                                                             motif='CG', strand='+')
    index_position_stub_cg2 = StratifiedBedIndexPositionStub(chrom='11', start=11, end=12,
                                                             motif='CG', strand='-')
    index_position_stub_chh = StratifiedBedIndexPositionStub(chrom='11', start=100, end=101,
                                                             motif='CHH', strand='+')

    cg_motif_pileup_stub1 = StratifiedBedMotifPileupStub(
            meth_counts_arr=np.array([[0, 2, 1, 2, 1, 4, 5],
                                      [0, 2, 2, 2, 3, 4, 7]], dtype='i4').T,
            idx_pos=index_position_stub_cg1)
    cg_motif_pileup_stub2 = StratifiedBedMotifPileupStub(
            meth_counts_arr=np.array([[0, 2, 1, 2, 1, 4, 5],
                                      [0, 2, 2, 4, 3, 6, 9]], dtype='i4').T,
            idx_pos=index_position_stub_cg2)
    chh_motif_pileup_stub = StratifiedBedMotifPileupStub(
            meth_counts_arr=np.array([[1, 2, 1, 2, 2, 4, 6],
                                      [1, 2, 2, 2, 4, 4, 8]], dtype='i4').T,
            idx_pos=index_position_stub_chh
    )

    # for interactive development:
    # output_by_chrom_motif = '/home/stephen/temp' + (
    #     '/subdir-should-be-created/recursively/strat-calls_[motif]_[chrom].bed.gz')
    output_by_chrom_motif = str(tmpdir.join('subdir-should-be-created/recursively'
                                            '/strat-calls_[motif]_[chrom].bed.gz'))

    strat_bed_writer = StratifiedBedWriter(calls_by_chrom_motif_fp=output_by_chrom_motif,
                                           motifs=['CG', 'CHH'], chrom='11')
    strat_bed_writer.setup()
    strat_bed_writer.process(cast(MotifPileup, cg_motif_pileup_stub1))
    strat_bed_writer.process(cast(MotifPileup, cg_motif_pileup_stub2))
    strat_bed_writer.process(cast(MotifPileup, chh_motif_pileup_stub))
    del strat_bed_writer  # to flush output files

    with gzip.open(output_by_chrom_motif.replace('[chrom]', '11')
                   .replace('[motif]', 'CG'), 'rt') as fin:
        computed_cg_header, computed_cg_line1, computed_cg_line2 = fin.readlines()

    with gzip.open(output_by_chrom_motif.replace('[chrom]', '11')
                   .replace('[motif]', 'CHH'), 'rt') as fin:
        computed_chh_header, computed_chh_line = fin.readlines()
    # -

    expected_header = '\t'.join(
            ['#chrom', 'start', 'end', 'motif', 'score', 'strand',
             'beta_value', 'n_meth', 'n_total',
             'c_bc_beta_value', 'c_bc_n_meth', 'c_bc_n_total',
             'c_bc_rv_beta_value', 'c_bc_rv_n_meth', 'c_bc_rv_n_total',
             'w_bc_beta_value', 'w_bc_n_meth', 'w_bc_n_total',
             'w_bc_rv_beta_value', 'w_bc_rv_n_meth', 'w_bc_rv_n_total',
             'mate1_beta_value', 'mate1_n_meth', 'mate1_n_total',
             'mate2_beta_value', 'mate2_n_meth', 'mate2_n_total',
             ]) + '\n'
    expected_cg_line1 = compute_expected_calls_str(cg_motif_pileup_stub1)
    expected_cg_line2 = compute_expected_calls_str(cg_motif_pileup_stub2)
    expected_chh_line = compute_expected_calls_str(chh_motif_pileup_stub)

    assert computed_cg_header == computed_chh_header == expected_header
    assert computed_cg_line1 == expected_cg_line1
    assert computed_cg_line2 == expected_cg_line2
    assert computed_chh_line == expected_chh_line


def compute_expected_calls_str(motif_pileup_stub: StratifiedBedMotifPileupStub) -> str:
    idx_pos = motif_pileup_stub.idx_pos
    expected_cg_calls_idx_str = (f'{idx_pos.chrom} {idx_pos.start} {idx_pos.end}'
                                 f' {idx_pos.motif} . {idx_pos.strand}').replace(' ', '\t')
    expected_cg_calls_calls_substr1 = '\t'.join([
        f'{motif_pileup_stub.strat_beta_arr[scall_ids.all]:.8f}',
        str(motif_pileup_stub.meth_counts_arr[scall_ids.all, mstat_ids.n_meth]),
        str(motif_pileup_stub.meth_counts_arr[scall_ids.all, mstat_ids.n_total])
    ])
    expected_cg_calls_calls_substr2 = '\t'.join(np.hstack([
        np.char.mod('%.8f', motif_pileup_stub.strat_beta_arr[0:scall_ids.all, np.newaxis]),
        motif_pileup_stub.meth_counts_arr[0:scall_ids.all, [mstat_ids.n_meth]].astype(str),
        motif_pileup_stub.meth_counts_arr[0:scall_ids.all, [mstat_ids.n_total]].astype(str)
    ]).flatten())
    expected_cg_calls_str = (expected_cg_calls_idx_str + '\t' + expected_cg_calls_calls_substr1
                             + '\t' + expected_cg_calls_calls_substr2 + '\n')
    return expected_cg_calls_str
