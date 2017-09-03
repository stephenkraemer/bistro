"""
input:

"""
import pytest

from mqc.overlap import OverlapHandler
from mqc.flag_and_index_values import methylation_status_flags, qc_fail_flags
mflags = methylation_status_flags
qflags = qc_fail_flags

class AlignedSegmentStub:
    def __init__(self, query_name):
        self.query_name = query_name

class PileupReadStub:
    def __init__(self, name, mflag, phred=30, qcflag=0, oflag=0, trimm_flag=0):
        self.overlap_flag = oflag
        self.qc_fail_flag = qcflag
        self.trimm_flag = trimm_flag
        self.meth_status_flag = mflag
        self.baseq_at_pos = phred
        self.alignment = AlignedSegmentStub(name)
pread = PileupReadStub

class MotifPileupStub:
    def __init__(self, reads):
        self.reads = reads

READS_DICT = {'name1_1': pread('name1', mflags.is_unmethylated),
              'name1_2': pread('name1', mflags.is_methylated),
              'name2_1': pread('name2', mflags.is_methylated),
              'name3_1': pread('name3', mflags.is_unmethylated),
              'name4_1': pread('name4', mflags.is_unmethylated),
              'name4_2': pread('name4', mflags.is_methylated),
              'name5_1': pread('name5', mflags.is_methylated, phred=30),
              'name5_2': pread('name5', mflags.is_methylated, phred=10),
              'name6_1': pread('name6', mflags.is_unmethylated, phred=20),
              'name6_2': pread('name6', mflags.is_unmethylated, phred=40),
              # no overlap in next pair, because one read is trimmed
              'name7_1': pread('name7', mflags.is_unmethylated, trimm_flag=1),
              'name7_2': pread('name7', mflags.is_unmethylated, trimm_flag=0),
              # na meth reads are not considered
              'name8_1': pread('name8', mflags.is_unmethylated),
              'name8_2': pread('name8', mflags.is_na),
              # reads with pre-existing qc fail flag are not considered
              'name9_1': pread('name9', mflags.is_unmethylated, qcflag=1),
              'name9_2': pread('name9', mflags.is_unmethylated, qcflag=0),
              }

MOTIF_PILEUP = MotifPileupStub(READS_DICT.values())
OVERLAP_HANDLER = OverlapHandler()
OVERLAP_HANDLER.process(MOTIF_PILEUP)


def test_reads_without_overlap_are_not_touched():
    assert READS_DICT['name2_1'].qc_fail_flag == 0
    assert READS_DICT['name2_1'].overlap_flag == 0
    assert READS_DICT['name3_1'].qc_fail_flag == 0
    assert READS_DICT['name3_1'].overlap_flag == 0


def test_mismatched_overlaps_are_tagged_with_qc_fail_flag():
    assert READS_DICT['name1_1'].qc_fail_flag == qflags.overlap_fail
    assert READS_DICT['name1_2'].qc_fail_flag == qflags.overlap_fail
    assert READS_DICT['name4_1'].qc_fail_flag == qflags.overlap_fail
    assert READS_DICT['name4_2'].qc_fail_flag == qflags.overlap_fail


def test_for_overlaps_with_matching_meth_status_the_read_with_lower_phred_gets_the_overlap_flag():
    assert READS_DICT['name5_1'].overlap_flag == 0
    assert READS_DICT['name5_2'].overlap_flag == 1
    assert READS_DICT['name6_1'].overlap_flag == 1
    assert READS_DICT['name6_2'].overlap_flag == 0

def test_reads_with_preexisting_qcfail_flag_are_not_considered():
    assert READS_DICT['name9_1'].overlap_flag == 0
    assert READS_DICT['name9_2'].overlap_flag == 0

def test_reads_with_NAmeth_are_not_considered():
    assert READS_DICT['name8_1'].overlap_flag == 0
    assert READS_DICT['name8_2'].overlap_flag == 0

def test_reads_with_trimm_flag_are_not_considered():
    assert READS_DICT['name7_1'].overlap_flag == 0
    assert READS_DICT['name7_2'].overlap_flag == 0

