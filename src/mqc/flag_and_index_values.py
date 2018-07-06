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

# noinspection PyPep8

from typing import NamedTuple

bsseq_strand_na_index = -1

# The indices for the BS-strands must be kept in sync between strat_call_indices and
# bsseq_strand_indices, as I have code relying on this
class BsseqIndex(NamedTuple):
    c_bc: int = 0
    c_bc_rv: int = 1
    w_bc: int = 2
    w_bc_rv: int = 3
bsseq_strand_indices = BsseqIndex()


# The indices for the BS-strands must be kept in sync between strat_call_indices and
# bsseq_strand_indices, as I have code relying on this
class StratCallIndices(NamedTuple):
    c_bc: int = 0
    c_bc_rv: int = 1
    w_bc: int = 2
    w_bc_rv: int = 3
    mate1: int = 4
    mate2: int = 5
    all: int = 6
strat_call_indices = StratCallIndices()


class MethStatIndices(NamedTuple):
    n_meth: int = 0
    n_total: int = 1
meth_status_indices = MethStatIndices()


class MethylationStatusFlags(NamedTuple):
    is_na: int = 16
    is_methylated: int = 8
    is_unmethylated: int = 4
    is_snp: int = 2
    is_ref: int = 1
methylation_status_flags = MethylationStatusFlags()


class QcFailFlags(NamedTuple):
    sam_flag_fail: int = 1
    phred_score_fail: int = 2
    mapq_fail: int = 4
    missing_info_fail: int = 8
    softclipped: int = 16
    pos_in_read_exceeds_tlen: int = 32
    overlap_fail: int = 64
qc_fail_flags = QcFailFlags()
