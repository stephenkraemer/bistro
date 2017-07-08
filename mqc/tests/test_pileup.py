import mqc
from mqc.tests.conftest import bam_path, index_file_path, sample_name

class SimpleCounter:
    def __init__(self):
        self.coverage = []
        pass

    def process(self, motif_pileup: mqc.MotifPileup):
        len(motif_pileup.all_pileupreads_generator())

simple_counter = SimpleCounter()

pileup_iterator = mqc.StepwisePileupIterator(
        bam_path=bam_path,
        index_file_path=index_file_path,
        taggers=[],
        counters=[simple_counter]
)

