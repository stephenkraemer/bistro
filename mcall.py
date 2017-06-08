"""Add ndarray with methylation calls to MotifPileup

This module is designed to be cythonized later

The choice of using an ndarray to represent methylation levels instead
of hiding the implementation in a MethLevels object is deliberate. The
gained simplicity is worth the coupling and problems in refactoring.

I may use a typing alias to help with refactoring.

Methylation calls may be done
- for the whole motif
- per position
- both

Methylation calls may be stratified
- by "all reads"
- by sequencing strand
- by mate

The ndarray will always be a multidimensional array containing counts for
all these scenarios, but only the ones specified to be computed will actually
be filled

It is up to the clients to make sure to only read values which have actually
been computed. Again this design was chosen deliberately because it makes
the code much simpler, especially with regard to porting this code and
its clients to cython later on

Currently, only motifs with up to three bases are allowed, therefore the shape
of the array is (4, 7, 3), where the index names and levels are in order:

position
  - motif
  - position1 
  - position2
  - position3

stratification level
  - all reads
  - mate1
  - mate2
  - c
  - crv
  - w
  - wrv

meth parameter
  - n_meth
  - n_total
  - beta
"""

import mqc
import numpy as np
from typing import Iterable
import mqc.flag_and_index_values.methylation_status_flags as mflags

all_reads_idx = 0
mate_idx_offset = 1
bsseq_strand_idx_offset = 3

meth_idx = 0
# currently, we first use an unmeth count column and then reuse it as
# total counts column
unmeth_idx = 1
total_idx = 1

def call_meth_level(motif_pileup: mqc.MotifPileup,
                    do_call_per_mate,
                    do_call_per_bsseq_strand,
                    do_call_per_pos):
    """Compute methylation levels for a pileup position"""
    # TODO-performance: currently, the total motif meth levels are always computed

    meth_calls_arr = np.array([4, 7, 3], dtype=np.uint16)

    if do_call_per_pos:
        # if we call per position, we will compute the motif meth levels
        # based on the per position calls after the loop
        for pos_in_motif, watson_base, pileup_reads in \
                motif_pileup.column_wise_generator():
            call_meth_level_foo(
                    pileup_reads,
                    meth_calls_arr=meth_calls_arr,
                    pos_idx=pos_in_motif + 1,
                    do_call_per_mate=do_call_per_mate,
                    do_call_per_bsseq_strand=do_call_per_bsseq_strand)
        # TODO-performance: computation is done across all strata
        #                   perhaps better to explicitely count totals too
        meth_calls_arr[0] = meth_calls_arr[1:].sum()

    else:  # only call on the whole motif
        call_meth_level_foo(motif_pileup.all_pileupreads_generator(),
                            meth_calls_arr = meth_calls_arr,
                            pos_idx = 0,
                            do_call_per_mate=do_call_per_mate,
                            do_call_per_bsseq_strand=do_call_per_bsseq_strand)



def call_meth_level_foo(pileup_reads: Iterable[mqc.BSSeqPileupRead],
                        meth_calls_arr,
                        pos_idx,
                        do_call_per_mate,
                        do_call_per_bsseq_strand) -> None:

    """On the computation of the n_total columns
    
    There are two alternative solutions:
    1. Use a function (inlined when ported to cython)
       def update_meth_calls_arr(meth_calls_arr, pos_idx, mate_idx, meth_status)
       - the function would check if we have an meth or unmeth event
       - it would update only the meth column or the meth and total column 
    
    2. Use the meth status indices (0 = meth, 1 = unmeth) to first fill the 
       n_total column with unmeth events, and sum up n_meth and n_unmeth columns
       at the end. <-- this is the current solution, but upon porting the code
       to cython, solution 1 shall be implemented
    """

    for read in pileup_reads:
        # TODO-refactor: with SNP calling in mind, change to two flags:
        #                meth_status and meth_level, the latter only defined
        #                if appropriate (or similar)

        has_meth_status = False
        if read.meth_status_flag == mflags.is_methylated:
            has_meth_status = True
            meth_status_idx = meth_idx
        elif read.meth_status_flag == mflags.is_unmethylated:
            has_meth_status = True
            meth_status_idx = unmeth_idx

        if not read.qc_fail_flag and not read.trimm_flag and has_meth_status:
            # for per mate and per strand calls, we don't discard overlapping reads
            if do_call_per_mate:
                meth_calls_arr[pos_idx, read.mate + mate_idx_offset, meth_status_idx] += 1
            if do_call_per_bsseq_strand:
                strat_idx = read.bsseq_strand_ind + bsseq_strand_idx_offset
                meth_calls_arr[pos_idx, strat_idx, meth_status_idx] += 1
            # for the "all reads" beta value, we discard overlapping reads
            if not read.overlap_flag:
                meth_calls_arr[pos_idx, all_reads_idx, meth_status_idx] += 1

    # will be replaced by inline function as detailed above when porting to cython
    meth_calls_arr[pos_idx, :, total_idx] = (
        meth_calls_arr[pos_idx, :, meth_idx] + meth_calls_arr[pos_idx, :, total_idx])

    return None
