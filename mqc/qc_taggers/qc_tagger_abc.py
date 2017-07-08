from abc import ABCMeta, abstractmethod
import mqc

class QcTagger(meta=ABCMeta):

    @abstractmethod
    def process(self, motif_pileup: mqc.MotifPileup) -> None:
        pass