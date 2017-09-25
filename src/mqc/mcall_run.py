"""Manage runs across a set of index positions"""

import multiprocessing as mp
import pickle
from abc import abstractmethod, ABCMeta
from collections import OrderedDict
from copy import deepcopy
from typing import Dict, List

import pysam
from mqc.coverage import CoverageCounter
from mqc.index import IndexFile
from mqc.mbias import MbiasCounter, FixedRelativeCuttingSites
from mqc.mcaller import MethCaller, StratifiedMethCaller
from mqc.beta_value import StratifiedBetaCounter
from mqc.overlap import OverlapHandler
from mqc.pileup.pileup import stepwise_pileup_generator
from mqc.qc_filters import PhredFilter, MapqFilter
from mqc.trimming import Trimmer
from mqc.visitors import Counter, Visitor
from mqc.writers import BedWriter


# from mqc.mbias import MbiasData, AdjustedMbiasCuttingSites, MbiasCounter
# from mqc.beta_values import StratifiedBetaValueCounter
# from mqc.coverage import CoverageCounter


def collect_stats(config):

    mbias_run = MbiasDeterminationRun(config, cutting_sites=None)
    mbias_run.run_parallel()
    for counter in mbias_run.summed_up_counters.values():
        counter.save_dataframe()


def run_mcalling(config):
    if config['run']['use_mbias_fit']:
        with open(config['paths']['adjusted_cutting_sites_obj_p'], 'rb') as fin:
            cutting_sites = pickle.load(fin)
    else:
        cutting_sites =  FixedRelativeCuttingSites(config)
    second_run = QcAndMethCallingRun(config,
                                     cutting_sites=cutting_sites)
    second_run.run_parallel()

    for counter in second_run.summed_up_counters.values():
        counter.save_dataframe()

class PileupRun(metaclass=ABCMeta):
    """ABC for classes managing individual data collection runs

    Here, run means one pass over some or all index positions, creating one
    MotifPileup for every position and using this information to do futher
    processing.

    Every run is associated with a set of configuration variables, most notably
    there can only be one definition of the CuttingSites for every run.

    Attributes
    ----------
    summed_up_counters: dict
        Stores the results for all Counters used in the PileupRun, defaults
        to {} if no Counters were applied
    n_cores: int
        number of cores available for parallelization (from config)

    Notes
    -----
    *Parallelization:* The Runs may be parallelized over separate index
    files, commonly one index file per chromosome. One index file may
    contain different sequence contexts (currently [CHH, CHG, CG] - may
    become unrestricted in the future).

    *How it works:* Runs are realized using an Iterator-Visitor-like
    pattern. For simple cases, the following algorithm is used to perform
    the run: ::

        For all index positions:
            create a MotifPileup
            for curr_visitor in PileupRun.visitors:
                curr_visitor.process(MotifPileup)

    The algorithm described above is implemented in
     :func:`~mqc.PileupRun._run_over_index_file`.
    For simple tasks such as methylation calling, this is sufficient. To
    implement complex tasks, one may have to overwrite the basic
    implementation of :func:`~mqc.PileupRun._run_over_index_file`.

    *Notes on implementing PileupRun subclasses:* There is only one abstract
    method, which is used to set up the visitors (PileupRun._get_visitors).
    Often, this function will simply

    1. initialize visitors by calling their respective constructors with
       the config dict and potentially a CuttingSites instance as arguments
    2. return a list of visitors

    """

    def __init__(self, config, cutting_sites=None):

        self.summed_up_counters = {}

        # All config vars used directly by methods of this class
        self.n_cores = config['run']['cores']
        self.index_files = config['run']['index_files']
        self.bam_path = config['run']['bam']
        self.cutting_sites = cutting_sites

        # Save whole config dict, because it is used
        # by the visitor class constructors
        self.config = config

    def run_parallel(self):
        # TODO: other path for cores=1
        """Parallelize over index files

        Note
        ----
        - currently when performing tasks such as methylation calling,
          it is assumed that one index file will correspond to one chromosome.
          If this is not the case, parallelization of e.g. writing methylation
          calls to file becomes more complex, because one has to maintain the
          correct order of the output.

        """

        with mp.Pool(processes=self.n_cores) as pool:
            counter_dicts_per_idx_file = pool.map(
                self._single_run_over_index_file,
                self.index_files)
            self.sum_up_counters(counter_dicts_per_idx_file)


        # counter = self._single_run_over_index_file(self.index_files[0])
        # return counter

    def sum_up_counters(
            self, counter_dicts_per_idx_file: List[Dict[str, Counter]]):
        """Sum up Counters retrieved from run_parallel

        Parameters
        ----------
        counter_dicts_per_idx_file:
            Retrieved from the parallel processing of the individual index
            files. Because all Counters were instantiated by the same
            :func:~mqc.PileupRun._get_visitors() method, the names,
            types and attributes of the Counters returned from every
            individual computation can be expected to be identical
        """

        # Note that you can't simply use the first encountered Counter
        # instance (first_counter) to initialize summed_up_counters[curr_name]
        # This would mean that (id(summed_up_counters[curr_name].counter_array)
        #     == id(first_counter.counter_array)
        self.summed_up_counters = {name: deepcopy(counter)
                                   for name, counter in counter_dicts_per_idx_file[0].items()}

        for curr_counters_dict in counter_dicts_per_idx_file[1:]:
            for curr_name, curr_counter in curr_counters_dict.items():
                self.summed_up_counters[curr_name].counter_array += (
                    curr_counter.counter_array)


    def _single_run_over_index_file(self, index_file_path: str) -> \
            Dict[str, 'Counter']:
        """Iterate over index file, generate MotifPileups and pass to Visitors

        Parameters
        ----------
        index_file_path
            IndexFile instance representing (usually) chromosome indices.
            See documentation if you need to use other index 'levels',
            because this may be problematic (for the moment) when runnning
            computations in parallel

        Returns
        -------
        dict
            Dict of counters (note: note ordered!). While arbitrary
            combinations of visitors may be applied during the PileupRun,
            only counters will be returned for further processing.
        """

        # TODO: put the index file path template into one central place to be able to change it
        # TODO: solve this with regex?
        # index_file_path_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
        chrom = index_file_path.split('_')[-1].replace('.bed.gz', '')

        visitors: Dict[str, Visitor] = self._get_visitors(chrom)

        alignment_file = pysam.AlignmentFile(self.bam_path)
        index_file = IndexFile(index_file_path)

        motif_pileup_generator = stepwise_pileup_generator(
            index_positions=index_file,
            alignment_file=alignment_file)

        for motif_pileup in motif_pileup_generator:
            for curr_visitor in visitors.values():
                curr_visitor.process(motif_pileup)

        counters: Dict = _filter_for_counters(visitors)

        return counters

    @abstractmethod
    def _get_visitors(self, chrom) -> Dict[str, Visitor]:
        """Every run specifies its algorithm by providing a list of visitors

        Returns
        -------
        OrderedDict[Visitors]
            Returns ordered dict of new instances of different Visitors.
            The same Visitor may in principle be used multiple times, although
            that is not a common case.
        """
        pass


class MbiasDeterminationRun(PileupRun):
    """Go over all index positions and collect M-bias stats

    *Implementation notes*

    - overlap handling is not required because M-bias stats are stratified
    by sequencing strand. Overlapping read pairs represent two events on two
    different strands which should both be counted.
    - trimming is not required because we want to collect stats even in regions
    where we know that bias exists (for verification, if nothing else)

    Possible improvements
    ---------------------
    In the future, this class may be able sample the index positions to reduce
    the number of computations. (or this may be built into PileupRun...)
    """

    def _get_visitors(self, chrom) -> Dict[str, Visitor]:
        return OrderedDict(
            mbias_counter=MbiasCounter(self.config),
        )


class QcAndMethCallingRun(PileupRun):
    """Methylation calling with QC filtering and stats collection"""

    def _get_visitors(self, chrom) -> Dict[str, Visitor]:
        visitors = OrderedDict(
            mapq_filter=MapqFilter(self.config),
            trimmer=Trimmer(self.config, self.cutting_sites),
            overlap_handler=OverlapHandler(),
            phred_filter=PhredFilter(self.config),
        )

        if self.config['run']['strat_beta_dist']:
            visitors['meth_caller'] = StratifiedMethCaller()
            visitors['beta_counter'] = StratifiedBetaCounter(self.config)
        else:
            visitors['meth_caller'] = MethCaller()

        visitors['mcall_writer'] = BedWriter(self.config, chrom=chrom)
        visitors['coverage_counter'] = CoverageCounter(self.config)

        return visitors


def _filter_for_counters(visitors) -> Dict:
    """Given a dict of Visitors, return a dict of Counters"""
    counters = {visitor_name: visitor_obj
                for visitor_name, visitor_obj in visitors.items()
                if isinstance(visitor_obj, Counter)}
    # Note: must return empty dict when no counters present
    return counters
