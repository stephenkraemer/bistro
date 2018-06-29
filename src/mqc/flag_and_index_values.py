"""Enums for array indices used throughout the package

Note
----
Using namedtuples instead of python enums because of superior performance

For example:

    from enum import IntEnum
    from collections import namedtuple

    l = [1, 2, 3]

    bs_strand = IntEnum('bs_strand_idx', 'c w cr wr'.split())
    %timeit l[bs_strand.c]

    class BsStrand(IntEnum):
        c = 0
        w = 1
        cr = 2
        wr = 3
    %timeit l[BsStrand.c]

    bs_strand_namedtuple = namedtuple('bs_strand_idx2', 'c w cr wr'.split())(
        c = 0,
        w = 1,
        cr = 2,
        wr = 3
    )
    %timeit l[bs_strand_namedtuple.w]

    151 ns ± 6.07 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
    147 ns ± 3.43 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
    88 ns ± 4.17 ns per loop (mean ± std. dev. of 7 runs, 10000000 loops each)
"""

from collections import namedtuple

bsseq_strand_indices = namedtuple('BsseqIndex', 'c_bc c_bc_rv w_bc w_bc_rv')(
    c_bc=0,
    c_bc_rv=1,
    w_bc=2,
    w_bc_rv=3,
)
bsseq_strand_na_index = -1


# The indices for the BS-strands must be kept in sync between strat_call_indices and
# bsseq_strand_indices, as I have code relying on this
strat_call_indices = namedtuple('StratCallIndices', 'c_bc c_bc_rv w_bc w_bc_rv mate1 mate2 all')(
    c_bc=0,
    c_bc_rv=1,
    w_bc=2,
    w_bc_rv=3,
    mate1=4,
    mate2=5,
    all=6,
)

meth_status_indices = namedtuple('MethStatIndices', 'n_meth n_total')(
    n_meth=0,
    n_total=1
)

methylation_status_flags = namedtuple('MethylationStatusFlags',
                                      ['is_na', 'is_methylated', 'is_unmethylated', 'is_snp', 'is_ref'])(
    is_na=16,
    is_methylated=8,
    is_unmethylated=4,
    is_snp=2,
    is_ref=1
)

qc_fail_flags = namedtuple('QcFailFlags', 'sam_flag_fail'
                                          ' phred_score_fail'
                                          ' mapq_fail'
                                          ' missing_info_fail'
                                          ' softclipped'
                                          ' pos_in_read_exceeds_tlen'
                                          ' overlap_fail')(
    sam_flag_fail = 1,
    phred_score_fail = 2,
    mapq_fail = 4,
    missing_info_fail = 8,
    softclipped = 16,
    pos_in_read_exceeds_tlen = 32,
    overlap_fail = 64
)
