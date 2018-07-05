from typing import List

from pysam import PileupRead, PileupColumn

class BSSeqPileupRead(PileupRead):

    @property
    def baseq_at_pos(self) -> int: ...

    @property
    def overlap_flag(self) -> int: ...

    @overlap_flag.setter
    def overlap_flag(self, value: int): ...

    @property
    def trimm_flag(self) -> int: ...


    @trimm_flag.setter
    def trimm_flag(self, value: int): ...

    @property
    def bsseq_strand_ind(self) -> int: ...

    @property
    def meth_status_flag(self) -> int: ...

    @property
    def observed_watson_base(self) -> int: ...

    @property
    def qc_fail_flag(self) -> int: ...

    @qc_fail_flag.setter
    def qc_fail_flag(self, value: int): ...

    @property
    def pos_in_read(self) -> int: ...

def pileups(pileup_column: PileupColumn, expected_watson_base: str) -> List[BSSeqPileupRead]: ...
