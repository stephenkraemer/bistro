"""
input:
- sequence of IndexPositions
- BAM file

output
- one MotifPileup per IndexPosition
- with reads=List[BSSeqPileupread]   if coverage
       reads = []                    else

scenarios

if coverage starts after first index position[s]:
    return empty MotifPileups until first index position
if coverage ends before last index position[s]:
   return empty MotifPileups for last index positions

if coverage is lost between index positions:
   correctly return MotifPileups

if coverage stretches beyond index positions in both directions:
   correctly return MotifPileups

integration problem

pysam.AlignmentFile.pileup
- skips 0 coverage positions
- reference position is 0 based
- takes start argument. what happens if there is no alignment at start (most probably silent skip)

How to test

- 2 index positions
- PileupColumnStubs, following each other, with gaps in between

index according to scenario

expected behavior

MotifPileup is called with - the appropriate index positions
the PileupColumnStubs.ref_position matches the index position
if PileColumnStubs: pileups(PileupColumnStubs)
MotifPIleup(pileups(PileupColumnStubs, index poistion)
else:
MotifPile(None, index position)
"""
from unittest.mock import MagicMock

from mqc.pileup.pileup import stepwise_pileup_generator

class PileupColumnStub:
    def __init__(self, ref_pos):
        self.reference_pos = ref_pos

class IndexPositionStub:
    def __init__(self, start):
        self.chrom = 'chr1',
        self.start = start
        self.end = start + 1
        self.watson_base = 'C'

class MotifPileupStub:
    def __init__(self, idx_pos, reads):
        self.reads = reads
        self.idx_pos = idx_pos
    def __eq__(self, other):
        reads_equal = (self.reads == other.reads)
        idx_pos_equal =  (self.idx_pos == other.idx_pos)
        if reads_equal and idx_pos_equal:
            return True
        else:
            return False

class AlignmentFileStub:
    def __init__(self, coverage):
        self.coverage = coverage

    def pileup(self, *args, **kwargs):
        pileup_columns = [PileupColumnStub(i)
                          for i, cov in enumerate(self.coverage)
                          if cov > 0]
        return pileup_columns

class TestMotifPileupGenerator:
    def test_motif_pileups_are_found_when_coverage_stretches_beyond_index_positions_at_both_sides(self, mocker, monkeypatch):


        bsseq_pileups_mock = mocker.patch('mqc.pileup.pileup.pileups')
        bsseq_pileups_mock.side_effect = lambda pileup_column_stub, watson_base: 'pileups_' + str(pileup_column_stub.reference_pos)
        monkeypatch.setattr('mqc.pileup.pileup.MotifPileup', MotifPileupStub)

        alignment_file = AlignmentFileStub(coverage = [1, 1, 1, 1, 1, 1, 1])
        index_positions = [IndexPositionStub(i) for i in (2, 4)]

        computed_res = list(stepwise_pileup_generator(index_positions=index_positions,
                                                      alignment_file=alignment_file))
        expected_res = [MotifPileupStub(idx_pos=idx_pos,
                                        reads = 'pileups_' + str(idx_pos.start))
                        for idx_pos in index_positions]

        assert computed_res == expected_res

    def test_empty_motif_pileups_are_returned_when_index_pos_start_or_end_in_regions_without_coverage(self, mocker, monkeypatch):

        bsseq_pileups_mock = mocker.patch('mqc.pileup.pileup.pileups')
        bsseq_pileups_mock.side_effect = lambda pileup_column_stub, watson_base: 'pileups_' + str(pileup_column_stub.reference_pos)
        monkeypatch.setattr('mqc.pileup.pileup.MotifPileup', MotifPileupStub)
        coverage = [0, 0, 0, 1, 1, 0, 0]
        alignment_file = AlignmentFileStub(coverage)
        index_positions = [IndexPositionStub(i) for i in (1, 3, 5, 6)]

        computed_res = list(stepwise_pileup_generator(index_positions=index_positions,
                                                      alignment_file=alignment_file))

        expected_res = []
        for idx_pos in index_positions:
            if coverage[idx_pos.start]:
                expected_res.append(MotifPileupStub(idx_pos=idx_pos,
                                                    reads='pileups_' + str(
                                                        idx_pos.start)))
            else:
                expected_res.append(MotifPileupStub(idx_pos=idx_pos,
                                                    reads=[]))

        assert computed_res == expected_res

    def test_empty_motif_pileups_are_returned_for_index_positions_within_uncovered_regions_within_the_bam(self, mocker, monkeypatch):
        bsseq_pileups_mock = mocker.patch('mqc.pileup.pileup.pileups')
        bsseq_pileups_mock.side_effect = lambda pileup_column_stub, watson_base: 'pileups_' + str(pileup_column_stub.reference_pos)
        monkeypatch.setattr('mqc.pileup.pileup.MotifPileup', MotifPileupStub)
        coverage = [1, 1, 1, 0, 0, 1, 1]
        alignment_file = AlignmentFileStub(coverage)
        index_positions = [IndexPositionStub(i) for i in (1, 4, 5)]

        computed_res = list(stepwise_pileup_generator(index_positions=index_positions,
                                                      alignment_file=alignment_file))

        expected_res = []
        for idx_pos in index_positions:
            if coverage[idx_pos.start]:
                expected_res.append(MotifPileupStub(idx_pos=idx_pos,
                                                    reads='pileups_' + str(
                                                        idx_pos.start)))
            else:
                expected_res.append(MotifPileupStub(idx_pos=idx_pos,
                                                    reads=[]))

        assert computed_res == expected_res
