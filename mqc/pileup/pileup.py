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
     1. currently automatically implemented through adjusted M-bias trimming, which discards small fragments
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

import itertools
import pysam
import mqc.pileup.bsseq_pileup_read as bsread
from mqc.pileup.bsseq_pileup_read import BSSeqPileupRead

from mqc.index import IndexPosition
from typing import Iterable, Tuple, List


class MotifPileup:
    """Access to all relevant information about an index position
    
    Provides BSSeqPileupReads at all C/G motif bases, with iterators to 
    either iterate over all reads in chained fashion, or to iteratore over
    the individual motif bases separately
    """

    def __init__(self,
                 motif_pileups,
                 idx_pos: IndexPosition):
        if motif_pileups is None:
            self.motif_pileups = [[]] * len(idx_pos.watson_motif)
        else:
            self.motif_pileups = motif_pileups
        self.idx_pos = idx_pos

    def column_wise_generator(self) -> Tuple[
        int, str, List[BSSeqPileupRead]]:
        for pos_idx, motif_base, pileup_reads in enumerate(zip(
                self.idx_pos.watson_motif, self.motif_pileups)):
            if motif_base in ['C', 'G']:
                yield (pos_idx, motif_base, pileup_reads)

    def all_pileupreads_generator(self):
        return itertools.chain.from_iterable(self.motif_pileups)


def stepwise_pileup_generator(index_file, bam_path) -> Iterable[MotifPileup]:
    """ Extract MotifPileups at index positions from the stream of pileups at all positions
    Iterates over all index positions on a chromosome. Where coverage is available,
    MotifPileup objects with the BSSeqPileupReads at the index positions are created.
    Positions without coverage are allowed and will result in emtpy MotifPileup objects.

    For every index position, the MotifPileup object is passed to the process()
    method of 1. all taggers (in order) and then 2. all counters (in order)
    """

    # Get pileups at all chromosome positions
    # Either get the PileupColumn from pysam/htslib
    # Or an empty pileup column (None) if the position has no coverage
    curr_idx_pos = next(index_file)  # first position in index file
    all_pileup_columns = all_position_pileup_generator(
        bam_path,
        ref=curr_idx_pos.chrom,
        start=curr_idx_pos.start)
    for pileup_column, curr_pileup_pos in all_pileup_columns:

        # Process hits
        # Either return a MotifPileup with PileupReads at all C/G positions
        # Or return an empty MotifPileup if the BAM file ends within
        # the current index position
        if curr_pileup_pos == curr_idx_pos.start:
            curr_base = curr_idx_pos.watson_motif[0]
            if curr_base in ['C', 'G'] and not pileup_column is None:
                motif_pileups = [
                    bsread.pileups(pileup_column)]
            else:
                motif_pileups = [[]]
            for curr_base in curr_idx_pos.watson_motif[1:]:
                try:
                    pileup_column, curr_pileup_pos = next(all_pileup_columns)
                except StopIteration:
                    # PileupColumn within the current index is not available
                    # Return empty MotifPileup for this index position.
                    yield MotifPileup(motif_pileups=None,
                                      idx_pos=curr_idx_pos)
                    # stop iteration over pileup columns
                    break
                if curr_base in ['C', 'G'] and not pileup_column is None:
                    motif_pileups.append(
                        bsread.pileups(pileup_column))
                else:
                    # pileup_column None or non C/G base
                    motif_pileups.append([])
            yield MotifPileup(motif_pileups=motif_pileups,
                              idx_pos=curr_idx_pos)

            try:
                curr_idx_pos = next(index_file)
            except StopIteration:  # end of index reached, end generator
                return

    # The pileup column generator is exhausted (end of BAM file)
    # But we still have index positions -> return empty MotifPileups
    # for all index positions
    for curr_idx_pos in index_file:
        yield MotifPileup(motif_pileups=None,
                          idx_pos=curr_idx_pos)


def all_position_pileup_generator(bam_path, ref, start):
    """ Return pileup results (pysam.PileupColumn or None) for all chromosome positions following the first index position

    The htslib pileup engine skips positions with 0 coverage
    This generator will yield a MotifPileup also at these positions
    """
    sam_file = pysam.AlignmentFile(bam_path)
    pileup_columns = sam_file.pileup(reference=ref,
                                     start=start,
                                     end=None,
                                     truncate=True)

    last_pos = start - 1
    for pileup_column in pileup_columns:
        curr_pos = pileup_column.reference_pos
        while last_pos + 1 < curr_pos:
            yield (None, last_pos + 1)
            last_pos += 1
        last_pos = curr_pos
        yield (pileup_column, curr_pos)

    sam_file.close()
