import numpy as np
from unittest.mock import MagicMock

from mqc.mcaller import MethCaller, StratifiedMethCaller

import mqc.flag_and_index_values as mfl

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
qflags = mfl.qc_fail_flags
strat_call_ids = mfl.strat_call_indices
mstat_ind = mfl.meth_status_indices

class MotifPileupStub:
    def __init__(self, reads, beta_value=None):
        self.reads = reads
        self.beta_value = beta_value
        # index position not required

class PileupreadStub:
    def __init__(self, meth_status_flag, bsseq_strand_ind,
                 read=1, trimm_flag=0, overlap_flag=0, qc_fail_flag=0,
                 phred_fail_flag=0):
        self.phred_fail_flag = phred_fail_flag
        self.trimm_flag = trimm_flag
        self.overlap_flag = overlap_flag
        self.meth_status_flag = meth_status_flag
        self.bsseq_strand_ind = bsseq_strand_ind
        self.qc_fail_flag = qc_fail_flag
        self.alignment = MagicMock(spec_set=['is_read1'])
        self.alignment.is_read1 = True if read == 1 else False

BASE_READ_PROPERTIES = [
    {'meth_status_flag'   : m_flags.is_unmethylated,
     'bsseq_strand_ind'   : b_inds.c_bc},
    {'meth_status_flag'   : m_flags.is_methylated,
     'bsseq_strand_ind'   : b_inds.c_bc_rv},
    {'meth_status_flag'   : m_flags.is_methylated,
     'bsseq_strand_ind'   : b_inds.w_bc},
    {'meth_status_flag'   : m_flags.is_unmethylated,
     'bsseq_strand_ind'   : b_inds.w_bc_rv},
]

BASE_READS = [PileupreadStub(**property_dict)
              for property_dict in BASE_READ_PROPERTIES]
BASE_MOTIF_PILEUP = MotifPileupStub(reads=BASE_READS)


# noinspection PyTypeChecker,PyUnresolvedReferences
class TestMethCaller:
    def test_computes_stats_from_usable_reads(self):
        meth_caller = MethCaller()
        meth_caller.process(BASE_MOTIF_PILEUP)
        assert BASE_MOTIF_PILEUP.beta_value == 0.5
        assert BASE_MOTIF_PILEUP.n_meth == 2
        assert BASE_MOTIF_PILEUP.n_total == 4

    def test_discards_qc_fail_flag_and_phred_fail_flag(self):

        additional_read_properties = [
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.w_bc,
             'qc_fail_flag'       : qflags.mapq_fail},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'qc_fail_flag'       : qflags.sam_flag_fail},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'qc_fail_flag'       : qflags.sam_flag_fail,
             'phred_fail_flag'    : 1,
             },
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'phred_fail_flag'    : 1,
             },
        ]
        additional_reads = [PileupreadStub(**property_dict)
                            for property_dict in additional_read_properties]

        all_reads = BASE_READS + additional_reads

        motif_pileup_with_qc_fails = MotifPileupStub(reads=all_reads)

        meth_caller = MethCaller()
        meth_caller.process(motif_pileup_with_qc_fails)

        assert motif_pileup_with_qc_fails.beta_value == 0.6
        assert motif_pileup_with_qc_fails.n_meth == 3
        assert motif_pileup_with_qc_fails.n_total == 5

    def test_discards_reads_with_overlap_flag(self):

        additional_read_properties = [
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.w_bc_rv},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.w_bc_rv,
             'overlap_flag'       : 1},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc,
             'overlap_flag'       : 1},
        ]

        additional_reads = [PileupreadStub(**property_dict)
                            for property_dict in additional_read_properties]
        all_reads = BASE_READS + additional_reads

        motif_pileup_with_overlap_flags = MotifPileupStub(reads=all_reads)

        meth_caller = MethCaller()
        meth_caller.process(motif_pileup_with_overlap_flags)

        assert motif_pileup_with_overlap_flags.n_meth == 3
        assert motif_pileup_with_overlap_flags.n_total == 6
        assert motif_pileup_with_overlap_flags.beta_value == 3/6

    def test_discards_unusable_meth_status_flags(self):

        additional_read_properties = [
            {'meth_status_flag'   : m_flags.is_na,
             'bsseq_strand_ind'   : b_inds.w_bc_rv},
            {'meth_status_flag'   : m_flags.is_ref,
             'bsseq_strand_ind'   : b_inds.w_bc},
            {'meth_status_flag'   : m_flags.is_snp,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc}
        ]

        additional_reads = [PileupreadStub(**property_dict)
                            for property_dict in additional_read_properties]
        all_reads = BASE_READS + additional_reads

        motif_pileup = MotifPileupStub(reads=all_reads)

        meth_caller = MethCaller()
        meth_caller.process(motif_pileup)

        assert motif_pileup.n_meth == 3
        assert motif_pileup.n_total == 5
        assert motif_pileup.beta_value == 3/5

    def test_discards_trimmed_events(self):
        additional_read_properties = [
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.w_bc},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc,
             'trimm_flag'         : 1},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'trimm_flag'         : 1},
        ]

        additional_reads = [PileupreadStub(**property_dict)
                            for property_dict in additional_read_properties]
        all_reads = BASE_READS + additional_reads

        motif_pileup = MotifPileupStub(reads=all_reads)

        meth_caller = MethCaller()
        meth_caller.process(motif_pileup)

        assert motif_pileup.n_meth == 2
        assert motif_pileup.n_total == 7
        assert motif_pileup.beta_value == 2/7


# noinspection PyUnresolvedReferences

class TestStratifiedMethCaller:

    qc_issue_read_list = [
            # bad: trimm + qc
            {'meth_status_flag' : m_flags.is_methylated,
             'bsseq_strand_ind' : b_inds.w_bc_rv,
             'read'             : 2,
             'trimm_flag'       : 1,
             'qc_fail_flag'     : 1},

            {'meth_status_flag' : m_flags.is_unmethylated,
             'bsseq_strand_ind' : b_inds.c_bc,
             'read'          : 1},
            {'meth_status_flag' : m_flags.is_methylated,
             'bsseq_strand_ind' : b_inds.c_bc_rv,
             'read'          : 1},

            # bad: qc
            {'meth_status_flag' : m_flags.is_methylated,
             'bsseq_strand_ind' : b_inds.c_bc,
             'read'             : 2,
             'qc_fail_flag'     : 1},

            {'meth_status_flag' : m_flags.is_methylated,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'          : 1},
            {'meth_status_flag' : m_flags.is_unmethylated,
             'bsseq_strand_ind' : b_inds.w_bc_rv,
             'read'          : 2},

            # bad: trimm
            {'meth_status_flag' : m_flags.is_unmethylated,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1,
             'trimm_flag'       : 1},

            # bad: phred fail
            {'meth_status_flag' : m_flags.is_unmethylated,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1,
             'phred_fail_flag'       : 1},


            # bad: phred fail, qc fail
            {'meth_status_flag' : m_flags.is_unmethylated,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1,
             'phred_fail_flag'  : 1,
             'qc_fail_flag'     : 1,
             },

            # bad: meth status
            {'meth_status_flag' : m_flags.is_snp,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1},
            {'meth_status_flag' : m_flags.is_ref,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1},
            {'meth_status_flag' : m_flags.is_na,
             'bsseq_strand_ind' : b_inds.w_bc,
             'read'             : 1},

        ]

    def test_creates_strat_mcalls_array_discarding_bad_reads(self):
        """Only counts reliable meth. calls

        Events are discard if
        - they are NA, Ref, SNP
        - they have a qc_fail_flag, phred_fail_flag, trimm_flag

        The overlap flag is treated specially due to the stratification
        of the calls, and has a dedicated test below
        """
        meth_caller = StratifiedMethCaller()

        read_properties = self.qc_issue_read_list

        reads = [PileupreadStub(**property_dict) for property_dict in read_properties]
        motif_pileup = MotifPileupStub(reads)

        # noinspection PyTypeChecker
        meth_caller.process(motif_pileup)

        assert motif_pileup.meth_counts_arr[strat_call_ids.mate2, mstat_ind.n_meth] == 0
        assert motif_pileup.meth_counts_arr[strat_call_ids.w_bc_rv, mstat_ind.n_total] == 1
        assert motif_pileup.meth_counts_arr[strat_call_ids.c_bc, mstat_ind.n_meth] == 0
        assert motif_pileup.meth_counts_arr[strat_call_ids.mate1, mstat_ind.n_meth] == 2
        assert motif_pileup.meth_counts_arr[strat_call_ids.mate1, mstat_ind.n_total] == 3

    def test_overlapping_reads_handling_depends_on_stratum(self):
        """Overlapping reads are both counted for their respective
        mate and bs_strand, but counted only once for the 'all' stratum"""
        meth_caller = StratifiedMethCaller()

        read_properties = [
                            {'meth_status_flag' : m_flags.is_methylated,
                             'bsseq_strand_ind' : b_inds.c_bc,
                             'read'          : 1,
                             'overlap_flag': 1
                             },
                           ]

        reads = [PileupreadStub(**property_dict) for property_dict in read_properties]
        motif_pileup = MotifPileupStub(reads)

        # noinspection PyTypeChecker
        meth_caller.process(motif_pileup)

        assert motif_pileup.meth_counts_arr[strat_call_ids.all, mstat_ind.n_meth] == 0
        assert motif_pileup.meth_counts_arr[strat_call_ids.all, mstat_ind.n_total] == 0
        assert motif_pileup.meth_counts_arr[strat_call_ids.c_bc, mstat_ind.n_meth] == 1
        assert motif_pileup.meth_counts_arr[strat_call_ids.c_bc, mstat_ind.n_total] == 1
        assert motif_pileup.meth_counts_arr[strat_call_ids.mate1, mstat_ind.n_meth] == 1
        assert motif_pileup.meth_counts_arr[strat_call_ids.mate1, mstat_ind.n_total] == 1


    def test_creates_beta_value_array_with_zero_division_handling(self):
        """Zero division is saved as nan, not as inf"""

        meth_caller = StratifiedMethCaller()

        read_properties = [
                            {'meth_status_flag' : m_flags.is_unmethylated,
                             'bsseq_strand_ind' : b_inds.c_bc,
                             'read'          : 1},
                            {'meth_status_flag' : m_flags.is_methylated,
                             'bsseq_strand_ind' : b_inds.c_bc,
                             'read'          : 1},
                            {'meth_status_flag' : m_flags.is_methylated,
                             'bsseq_strand_ind' : b_inds.w_bc,
                             'read'          : 1},
                            {'meth_status_flag' : m_flags.is_unmethylated,
                             'bsseq_strand_ind' : b_inds.w_bc_rv,
                             'read'          : 2},
                           ]

        reads = [PileupreadStub(**property_dict) for property_dict in read_properties]
        motif_pileup = MotifPileupStub(reads)
        # noinspection PyTypeChecker
        meth_caller.process(motif_pileup)

        assert motif_pileup.strat_beta_arr[strat_call_ids.c_bc]    == 0.5
        assert motif_pileup.strat_beta_arr[strat_call_ids.mate1]   == 2/3
        assert motif_pileup.strat_beta_arr[strat_call_ids.w_bc_rv] == 0    # 0  nmeth/ 1 ntotal
        assert motif_pileup.strat_beta_arr[strat_call_ids.all]     == 0.5
        assert np.isnan(motif_pileup.strat_beta_arr[strat_call_ids.c_bc_rv])  # 0 / 0

    def test_global_meth_call_attributes_are_set(self):

        meth_caller = StratifiedMethCaller()

        read_properties = self.qc_issue_read_list

        reads = [PileupreadStub(**property_dict) for property_dict in read_properties]
        motif_pileup = MotifPileupStub(reads)
        # noinspection PyTypeChecker
        meth_caller.process(motif_pileup)

        assert motif_pileup.n_meth == 2
        assert motif_pileup.n_total == 4
        assert motif_pileup.beta_value == 0.5


