"""Test Trimmer"""
import numpy as np
import pandas as pd
import pytest

from mqc.trimming import Trimmer
from mqc.flag_and_index_values import (
    bsseq_strand_indices as b_inds,
    methylation_status_flags as mflags,
)

class AlignedSegmentStub:
    def __init__(self, tlen):
        self.template_length = tlen

class PileupReadStub:
    def __init__(self, pos, flen, strand,
                 exp_tr_flag,
                 mflag=mflags.is_methylated, qcflag=0, trimm_flag=0):
        self.meth_status_flag = mflag
        self.exp_tr_flag = exp_tr_flag
        self.pos_in_read = pos
        self.alignment = AlignedSegmentStub(flen)
        self.qc_fail_flag = qcflag
        self.trimm_flag = trimm_flag
        self.bsseq_strand_ind = strand
pread = PileupReadStub

class MotifPileupStub:
    def __init__(self, reads):
        self.reads = reads

class CuttingSitesStub:
    def __init__(self):
        # dataframe index used to find flen dimension
        self.df = pd.DataFrame(index=pd.MultiIndex.from_product(
            ['c_bc c_bc_rv'.split(), range(30)], names='bs_strand flen'.split()))

    # noinspection PyMethodMayBeStatic
    def as_array(self):
        arr = np.zeros((4, 110 + 1, 2))
        arr[b_inds.w_bc, 100, 0] = 10
        arr[b_inds.w_bc, 100, 1] = 91
        arr[b_inds.w_bc, 110, 0] = 10
        arr[b_inds.w_bc, 110, 1] = 101
        arr[b_inds.c_bc, 110, 0] = 0
        arr[b_inds.c_bc, 110, 1] = 81
        return arr

config_stub = {'trimming': {'max_flen_considered_for_trimming': 110}}
read_properties = [
    # trimming boundaries are zero-based slice definitions of the plateau
    # left-closed
    {'pos': 9, 'flen': 100, 'strand': b_inds.w_bc, 'exp_tr_flag': 1},
    {'pos': 10, 'flen': 100, 'strand': b_inds.w_bc, 'exp_tr_flag': 0},
    {'pos': 50, 'flen': 100, 'strand': b_inds.w_bc, 'exp_tr_flag': 0},

    # right-open
    {'pos': 100, 'flen': 110, 'strand': b_inds.w_bc, 'exp_tr_flag': 0},
    {'pos': 101, 'flen': 110, 'strand': b_inds.w_bc, 'exp_tr_flag': 1},

    # Different strands have different cuttings sites
    {'pos': 3, 'flen': 110, 'strand': b_inds.c_bc, 'exp_tr_flag': 0},
    {'pos': 90, 'flen': 110, 'strand': b_inds.c_bc, 'exp_tr_flag': 1},

    # the next two reads only differ in the qc_fail_flag or meth_na flag
    # the failing/NA read should also be processed
    {'pos': 100, 'flen': 100, 'strand': b_inds.w_bc, 'exp_tr_flag': 1},
    {'pos': 100, 'flen': 100, 'strand': b_inds.w_bc, 'qcflag': 1, 'exp_tr_flag': 1},
    {'pos': 100, 'flen': 100, 'strand': b_inds.w_bc, 'mflag': mflags.is_na, 'exp_tr_flag': 1},

    # this read exceeds max_flen
    {'pos': 101, 'flen': 200, 'strand': b_inds.w_bc, 'qcflag': 0, 'exp_tr_flag': 1},
]

base_reads = [pread(**i) for i in read_properties]

motif_pileup = MotifPileupStub(base_reads)

cutting_sites = CuttingSitesStub()
# noinspection PyTypeChecker
trimmer = Trimmer(cutting_sites=cutting_sites)
# noinspection PyTypeChecker
trimmer.process(motif_pileup)


@pytest.mark.parametrize('idx', range(8))
def test_out_of_bounds_read_positions_are_discarded(idx):
    # order of reads in MotifPileup.reads is currently not guaranteed
    # use original read list to find reads by indexing
    curr_read = base_reads[idx]
    assert curr_read.trimm_flag == curr_read.exp_tr_flag

@pytest.mark.parametrize('idx', (8, 9))
def test_qcfail_NAmeth_reads_are_also_processed(idx):
    # order of reads in MotifPileup.reads is currently not guaranteed
    # use original read list to find reads by indexing
    curr_read = base_reads[idx]
    assert curr_read.trimm_flag == curr_read.exp_tr_flag

@pytest.mark.parametrize('idx', (10,))
def test_flen_exceeding_max_flen_are_trimmed_like_max_flen(idx):
    curr_read = base_reads[idx]
    assert curr_read.trimm_flag == curr_read.exp_tr_flag

