"""Tests for QC fail filters

- phred filter
- mapq filter
"""
from typing import cast

from mqc.pileup.pileup import MotifPileup
from mqc.qc_filters import PhredFilter, MapqFilter

from mqc.flag_and_index_values import (
    qc_fail_flags as qflag,
    methylation_status_flags as mflags)


class MotifPileupStub:
    def __init__(self, reads):
        self.reads = reads


def test_phred_filter_tags_reads_with_below_threshold_phred_score_with_phred_fail_flag() -> None:
    """Test PhredFilter

    PhredFilter sets the read.phred_fail_flag. It does not modify the
    read.qc_fail_flag. The phred_fail_flag is set independent of other
    flag values of the read.
    """

    class PileupreadStub:
        def __init__(self, phred_score, overlap_flag=0, trimm_flag=0,
                     qc_fail_flag=0, phred_fail_flag=0):
            self.qc_fail_flag = qc_fail_flag
            self.phred_fail_flag = phred_fail_flag
            self.baseq_at_pos = phred_score
            self.overlap_flag = overlap_flag
            self.meth_status_flag = mflags.is_methylated
            self.trimm_flag = trimm_flag
    pread = PileupreadStub

    motif_pileup = MotifPileupStub([
        pread(phred_score=1), pread(phred_score=10), pread(phred_score=12),
        pread(phred_score=1, overlap_flag=1), pread(phred_score=10, trimm_flag=1),
        pread(20), pread(30)])

    config = { 'basic_quality_filtering': {'min_phred_score': 20}}
    phred_filter = PhredFilter(config)

    phred_filter.process(cast(MotifPileup, motif_pileup))
    for curr_read in motif_pileup.reads:
        if curr_read.baseq_at_pos < 20:
            assert curr_read.phred_fail_flag == 1 and curr_read.qc_fail_flag == 0
        else:
            assert curr_read.phred_fail_flag == 0 and curr_read.qc_fail_flag == 0

    motif_pileup = MotifPileupStub([
        pread(phred_score=1, qc_fail_flag=1),
        pread(phred_score=1, overlap_flag=1, qc_fail_flag=1),
        pread(20, qc_fail_flag=1)])

    config = { 'basic_quality_filtering': {'min_phred_score': 20}}
    phred_filter = PhredFilter(config)

    phred_filter.process(cast(MotifPileup, motif_pileup))
    for curr_read in motif_pileup.reads:
        if curr_read.baseq_at_pos < 20:
            assert curr_read.phred_fail_flag == 1 and curr_read.qc_fail_flag == 1
        else:
            assert curr_read.phred_fail_flag == 0 and curr_read.qc_fail_flag == 1

def test_mapq_filter_tags_reads_with_below_threshold_mapq_with_specific_qc_fail_flag() -> None:

    class AlignedSegmentStub:
        def __init__(self, mapq):
            self.mapping_quality = mapq

    class PileupreadStub:
        def __init__(self, mapq, qc_fail_flag=0):
            self.alignment = AlignedSegmentStub(mapq)
            self.qc_fail_flag = qc_fail_flag
    pread = PileupreadStub

    motif_pileup = MotifPileupStub([pread(1), pread(10), pread(20), pread(30)])

    config = {'basic_quality_filtering': {'min_mapq': 20}}
    mapq_filter = MapqFilter(config)

    mapq_filter.process(cast(MotifPileup, motif_pileup))

    for curr_read in motif_pileup.reads:
        if curr_read.alignment.mapping_quality >= 20:
            assert curr_read.qc_fail_flag & qflag.mapq_fail == 0
        else:
            assert curr_read.qc_fail_flag & qflag.mapq_fail > 0
