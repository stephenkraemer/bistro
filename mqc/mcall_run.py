from abc import abstractmethod
import mqc
import multiprocessing as mp
from functools import partial
from typing import List, Iterable, Callable

def mcall_run(config, index_file_paths, bam_path):
    first_run = MbiasDeterminationRun(config, cutting_sites=None)
    first_run.run_parallel()

    adjusted_cutting_sites =
    second_run = OneRun().run_parallel(writer=)
    join_run_results(first_run, second_run)



class PileupRun:
    def __init__(self, config, cutting_sites=None):
        self.config = config
        self.visitors = self.setup_visitors()

    def run_parallel(self):
        with mp.Pool(processes=config['parallelization']['n_cores']) as pool:
            pool_it = pool.imap(self.atomic_mcall_run_per_chromosome, self.index_file_paths)
            for chrom_counters in pool_it:
                self.counters = map(lambda x,y: x + y, self.counters, chrom_counters)

    def atomic_mcall_run_per_chromosome(self, index_file_path):

        index_file = mqc.index.IndexFile(index_file_path)
        motif_pileup_generator = mqc.pileup.stepwise_pileup_generator(
                index_file=index_file,
                bam_path=self.bam_path
        )

        for motif_pileup in motif_pileup_generator:
            for v in self.visitors:
                v.process(motif_pileup)

    @abstractmethod
    def setup_visitors(self):
        pass


class MbiasDeterminationRun(PileupRun):
    def setup_visitors(self):
        self.visitors = [
            mqc.MbiasCounter(self.config)
        ]

class QcAndMethCallingRun(OneRun):
    def set_up_visitors(self, adjusted_cutting_sites):
        self.visitors = [
            mqc.Trimmer(config, adjusted_cutting_sites),
            mqc.OverlapHandler(config),
            mqc.PhredScoreFilter(config),
            mqc.MethylationCaller(config),
            mqc.BetaValueCounter(config),
            mqc.CoverageCounter(config)
        ]

"""
def get_qc_taggers(cutting_sites, config):
    qc_taggers = [
        partial(mqc.trimm_reads, cutting_sites=mqc.MinimalCuttingSites[config['cutting_sites']]),
        partial(mqc.trimm_overlaps, config=config),
        partial(mqc.filter_by_phred_score, config=config)
    ]
    return qc_taggers

def get_counters(config, trimming_status, sample_metadata):
    # add sample metadata to output dataframes
    counters = [
        mqc.MbiasCounter(config),
        mqc.BetaValueCounter(config, trimming_status, sample),
        mqc.CoverageCounter(config, trimming_status, sample)
    ]
    return counters
"""
