from collections import namedtuple, defaultdict
from copy import deepcopy
from itertools import product

from mqc.mbias import MbiasCounter

import mqc.flag_and_index_values as mfl

#TODO: I currently test that any non-zero qc_fail_flag leads to discard from M-bias stats counting. When I update the behavior so that phred score fails are kept in the stats, the tests here also need to be updated accordingly
b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags

MAX_FLEN = 500
CONFIG = defaultdict(dict)
CONFIG['data_properties']['max_read_length_bp'] = 101
CONFIG['trimming']['max_flen_considered_for_trimming'] = MAX_FLEN
CONFIG['run']['motifs'] = ['cg', 'chg']

AlignmentStub = namedtuple('AlignmentStub',
                           'template_length')

IndexPositionStub = namedtuple('IndexPositionStub',
                               'motif')

class BSSeqPileupReadStub:
    def __init__(self, strand_idx, tlen, pos, mflag, fail):
        self.qc_fail_flag = fail
        self.bsseq_strand_ind = strand_idx
        self.meth_status_flag = mflag
        self.alignment = AlignmentStub(template_length = tlen)
        self.pos_in_read = pos
    # TODO: change to method
    def __iter__(self):
        if self.meth_status_flag == m_flags.is_methylated:
            meth_status_index = 0
        elif self.meth_status_flag == m_flags.is_unmethylated:
            meth_status_index = 1
        else:
            raise ValueError("This is a failed read, "
                             "can't give event_class indices")

        return iter((self.bsseq_strand_ind,
                     self.alignment.template_length,
                     self.pos_in_read,
                     meth_status_index))

bstub = BSSeqPileupReadStub

class MotifPileupStub:
    def __init__(self, motif, reads):
        self.idx_pos = IndexPositionStub(motif=motif)
        self.reads = reads

class TestMbiasCounterMotifPileupProcessing:

    def test_updates_counter_array(self):
        reads = [
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50,
                  mflag=m_flags.is_methylated, fail=0),
            bstub(strand_idx=b_inds.c_bc_rv, tlen=120, pos=50,
                  mflag=m_flags.is_unmethylated, fail=0)
        ]

        mbias_counter = MbiasCounter(CONFIG)

        # add reads to CG and CHG motif stratum
        for motif_ind, curr_read in product([0,1], reads):
                mbias_counter.counter_array[(motif_ind,) + tuple(curr_read)] = 2

        motif_pileup_cg = MotifPileupStub(motif = 'cg',
                                          reads=reads)
        motif_pileup_chg = MotifPileupStub(motif = 'chg',
                                           reads=reads)
        mbias_counter.process(motif_pileup_cg)
        mbias_counter.process(motif_pileup_chg)

        # assert counts for CG
        assert mbias_counter.counter_array[(0,) + tuple(reads[0])] == 3
        assert mbias_counter.counter_array[(0,) + tuple(reads[1])] == 3

        # assert counts for CHG
        assert mbias_counter.counter_array[(1,) + tuple(reads[0])] == 3
        assert mbias_counter.counter_array[(1,) + tuple(reads[1])] == 3

    def test_discards_reads_w_bad_qcflag_or_na_bsstrand_or_mcall_fail_ref_snp(self):
        """ Reads are discarded if

            - they have a qc_fail_flag
            - their bsseq strand could not be determined
            - they have methylation calling status: NA, SNP or Ref
        """
        reads = [
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_methylated, fail=1),  # qc fail
            bstub(strand_idx=b_na_ind,    tlen=100, pos=50, mflag=m_flags.is_methylated, fail=0), # bsseq strand not identifiable
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_na, fail=0), # na meth
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_snp, fail=0),  # snp meth
            bstub(strand_idx=b_inds.w_bc, tlen=100, pos=50, mflag=m_flags.is_ref, fail=0),  # ref meth
            bstub(strand_idx=b_inds.c_bc, tlen=200, pos=90, mflag=m_flags.is_unmethylated, fail=0),  # this one should count
        ]
        motif_pileup = MotifPileupStub(motif = 'cg', reads=reads)
        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.process(motif_pileup)
        assert mbias_counter.counter_array[(0,) + tuple(reads[-1])] == 1
        mbias_counter.counter_array[(0,) + tuple(reads[-1])] = 0
        assert (mbias_counter.counter_array == 0).all()

    def test_threshold_exceeding_flens_are_added_to_max_flen_bin(self):
        too_long_flen_read_properties = dict(
            strand_idx=b_inds.w_bc_rv, tlen=MAX_FLEN + 200, pos=20,
            mflag=m_flags.is_methylated, fail=0)
        reads = [bstub(**too_long_flen_read_properties)]
        motif_pileup = MotifPileupStub(motif = 'chg', reads=reads)
        mbias_counter = MbiasCounter(CONFIG)
        mbias_counter.process(motif_pileup)

        capped_flen_read_properties = deepcopy(too_long_flen_read_properties)
        capped_flen_read_properties['tlen'] = MAX_FLEN
        capped_read_idx = tuple(bstub(**capped_flen_read_properties))
        assert mbias_counter.counter_array[(1,) + capped_read_idx] == 1

