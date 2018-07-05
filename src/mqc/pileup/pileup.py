"""Classes to iteratively process index positions

Uses Iterator-Visitor design

An iterator class is chosen which collects information at every pileup position
in form of a MotifPileup. This MotifPileup is then handed to the process()
method of the different visitors of the iterator.

The visitors can be divided into two groups:
1. Taggers: taggers modify attributes of BSSeqPileupReads
2. Counters: counters iterate over the BSSeqPileupReads at a MotifPileup and
             gather, save and write out information

The Tagger.process() methods will be called in order of appearance in the
Iterator .taggers attribute. Often it will be important to pay attention to the
 correct order! E.g. you must tag read positions which have to be trimmed before
computing overlaps between mate1 and mate2.

A typical selection and order of taggers will be:
  1. Basic exclusion filters
     1. flag
     2. mapq
     3. refskip, indel or N at pileup position
     4. softclip at pileup position
     5. position outside of TLEN
  2. Mate exclusion filtering (if option is set)
  3. Fragment length filtering (discard small fragments?)
     1. currently automatically implemented through adjusted M-bias trimming,
        which discards small fragments
  4. Conversion error calling
  5. Trimming
     1. Constant base trimming (e.g. gap repairs)
     2. Adjusted trimming (e.g. adjusted M-bias trimming)
  6. Overlap handling (may adjust phred scores, depending on algorithm)
  7. Phred score filtering
  8. Methylation calling on remaining reads

Counters include DistributionCounters, such as the beta value, coverage, Mbias
and Phred score counter, but also the MethylationCaller class

Different iterators may provide different information through the BSSeqPileupreads
which are collected in the MotifPileup.

By default, a StepwisePileupIterator is provided, which only provides information about the
read aligning to the pileup site, and not about its mate. This corresponds to the standard htslib
pileup functionality.

Implementation of a CachedPileupIterator which would provide information
about both mates at the same time and about the global beta
values at all CH and CG positions in the fragment is planned but not done yet.
This would for example be helpful for more sophisticated conversion error calling
and qc filtering in general, as well as better fragment length determination
or some kinds of epipolymorphism analysis
"""
import numpy as np

import pysam
from itertools import chain, repeat
from typing import Iterator, List, NamedTuple

from mqc.index import IndexPosition
from mqc.pileup.bsseq_pileup_read import pileups, BSSeqPileupRead


class MotifPileup:

    meth_counts_arr: np.ndarray
    strat_beta_arr: np.ndarray
    beta_value: float
    n_meth: int
    n_total: int

    def __init__(self, reads: List[BSSeqPileupRead], idx_pos: IndexPosition) -> None:
        self.idx_pos = idx_pos
        self.reads = reads


def stepwise_pileup_generator(index_positions: Iterator[IndexPosition],
                              alignment_file: pysam.AlignmentFile,
                              ) -> Iterator[MotifPileup]:
    """MotifPileup Generator"""
    # pylint: disable=stop-iteration-return
    # (pylint false positive)

    idx_pos_iterable = iter(index_positions)
    curr_idx = next(idx_pos_iterable)

    pileup_columns = alignment_file.pileup(reference=curr_idx.chrom,
                                           start=curr_idx.start,
                                           end=None,
                                           truncate=True)

    # Functional form required
    # https://github.com/python/mypy/issues/4349
    # noinspection PyPep8Naming
    EmptyPileupColumn = NamedTuple('EmptyPileupColumn', [('reference_pos', int)])
    empty_pileup_column = EmptyPileupColumn(-1)
    empty_pileup_column_iterable = repeat(empty_pileup_column)

    pileup_columns_iterator = chain(pileup_columns,
                                    empty_pileup_column_iterable)
    curr_pileup_column = next(pileup_columns_iterator)
    curr_pileup_pos = curr_pileup_column.reference_pos

    while True:
        try:
            if curr_pileup_pos == -1:
                yield MotifPileup(reads=[], idx_pos=curr_idx)
                for curr_idx in idx_pos_iterable:
                    yield MotifPileup(reads=[], idx_pos=curr_idx)
                break
            elif curr_idx.start > curr_pileup_pos:
                curr_pileup_column = next(pileup_columns_iterator)
                curr_pileup_pos = curr_pileup_column.reference_pos
                continue
            elif curr_idx.start < curr_pileup_pos:
                yield MotifPileup(reads=[], idx_pos=curr_idx)
                curr_idx = next(idx_pos_iterable)
                continue
            elif curr_idx.start == curr_pileup_pos:
                pileup_reads = pileups(curr_pileup_column, curr_idx.watson_base)
                yield MotifPileup(reads=pileup_reads, idx_pos=curr_idx)
                curr_idx = next(idx_pos_iterable)
                curr_pileup_column = next(pileup_columns_iterator)
                curr_pileup_pos = curr_pileup_column.reference_pos
        except StopIteration:
            return  # generator will now raise StopIteration
