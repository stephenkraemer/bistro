from mqc.mcaller import MethCaller

import mqc.flag_and_index_values as mfl

b_inds = mfl.bsseq_strand_indices
b_na_ind = mfl.bsseq_strand_na_index
m_flags = mfl.methylation_status_flags
qflags = mfl.qc_fail_flags

class MotifPileupStub:
    def __init__(self, reads, beta_value=None):
        self.reads = reads
        self.beta_value = beta_value
        # index position not required

class PileupreadStub:
    def __init__(self, meth_status_flag, bsseq_strand_ind,
                 trimm_flag=0, overlap_flag=0, qc_fail_flag=0):
        self.trimm_flag = trimm_flag
        self.overlap_flag = overlap_flag
        self.meth_status_flag = meth_status_flag
        self.bsseq_strand_ind = bsseq_strand_ind
        self.qc_fail_flag = qc_fail_flag

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

class TestMethCaller:
    def test_computes_stats_from_usable_reads(self):
        meth_caller = MethCaller()
        meth_caller.process(BASE_MOTIF_PILEUP)
        assert BASE_MOTIF_PILEUP.beta_value == 0.5
        assert BASE_MOTIF_PILEUP.n_meth == 2
        assert BASE_MOTIF_PILEUP.n_unmeth == 2

    def test_discards_qc_fail_flag(self):

        additional_read_properties = [
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.c_bc},
            {'meth_status_flag'   : m_flags.is_methylated,
             'bsseq_strand_ind'   : b_inds.w_bc,
             'qc_fail_flag'       : qflags.mapq_fail},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'qc_fail_flag'       : qflags.phred_score_fail},
            {'meth_status_flag'   : m_flags.is_unmethylated,
             'bsseq_strand_ind'   : b_inds.c_bc_rv,
             'qc_fail_flag'       : qflags.sam_flag_fail},
        ]
        additional_reads = [PileupreadStub(**property_dict)
                            for property_dict in additional_read_properties]

        all_reads = BASE_READS + additional_reads

        motif_pileup_with_qc_fails = MotifPileupStub(reads=all_reads)

        meth_caller = MethCaller()
        meth_caller.process(motif_pileup_with_qc_fails)

        assert motif_pileup_with_qc_fails.beta_value == 0.6
        assert motif_pileup_with_qc_fails.n_meth == 3
        assert motif_pileup_with_qc_fails.n_unmeth == 2

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
        assert motif_pileup_with_overlap_flags.n_unmeth == 3
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
        assert motif_pileup.n_unmeth == 2
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
        assert motif_pileup.n_unmeth == 5
        assert motif_pileup.beta_value == 2/7
