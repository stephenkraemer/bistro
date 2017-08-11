from collections import namedtuple

bsseq_strand_indices = namedtuple('BsseqIndex', 'w_bc c_bc w_bc_rv c_bc_rv')(
    c_bc=0,
    c_bc_rv=1,
    w_bc=2,
    w_bc_rv=3,
)

bsseq_strand_na_index = -1

methylation_status_flags = namedtuple('MethylationStatusFlags',
                                      ['is_na', 'is_methylated', 'is_unmethylated', 'is_snp', 'is_ref'])(
    is_na=16,
    is_methylated=8,
    is_unmethylated=4,
    is_snp=2,
    is_ref=1
)

qc_fail_flags = namedtuple('QcFailFlags', 'sam_flag_fail phred_score_fail mapq_fail')(
    sam_flag_fail = 1,
    phred_score_fail = 2,
    mapq_fail = 4
)
