"""Manage runs across a set of index positions

Implementation notes
--------------------
- no OrderedDict type available for type annotation, using Dict
"""

import multiprocessing as mp
from abc import abstractmethod, ABCMeta
from collections import OrderedDict
from typing import Dict, List

from mqc.index import IndexFile
from mqc.pileup.pileup import stepwise_pileup_generator
from mqc.visitors import Counter, Visitor


# from mqc.mbias import MbiasData, AdjustedMbiasCuttingSites, MbiasCounter
# from mqc.beta_values import StratifiedBetaValueCounter
# from mqc.coverage import CoverageCounter


def mcall_run(config):
    """Conceptual draft of a function to perform methylation calling and QC"""
    first_run = MbiasDeterminationRun(config, cutting_sites=None)
    first_run.run_parallel()
    mbias_counter = first_run.summed_up_counters['mbias_counter']
    mbias_data = MbiasData(mbias_counter, config)
    adjusted_cutting_sites = AdjustedMbiasCuttingSites(
        mbias_data, calling_mode='adjusted', config=config)
    second_run = QcAndMethCallingRun(config,
                                     adjusted_cutting_sites.get_array())
    second_run.run_parallel()
    beta_value_counter = second_run.summed_up_counters['beta_counter']
    # further processing


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
        self.n_cores = config['parallelization']['n_cores']
        self.index_files = config['run']['index_files']
        if not isinstance(self.index_files[0], IndexFile):
            raise TypeError('Index files must be given as IndexFile objects')
        self.bam_path = config['run']['bam_path']
        self.cutting_sites = cutting_sites

        # Save whole config dict, because it is used
        # by the visitor class constructors
        self.config = config

    def run_parallel(self):
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

        # Set summed_up_counters to empty counters
        self.summed_up_counters = _filter_for_counters(self._get_visitors())
        for curr_counters_dict in counter_dicts_per_idx_file:
            for curr_name, curr_counter in curr_counters_dict.items():
                self.summed_up_counters[curr_name].counter_array += (
                    curr_counter.counter_array)

        # Note that you can't simply use the first encountered Counter
        # instance (first_counter) to initialize summed_up_counters[curr_name]
        # This would mean that (id(summed_up_counters[curr_name].counter_array)
        #     == id(first_counter.counter_array)

    def _single_run_over_index_file(self, index_file: IndexFile) -> \
            Dict[str, 'Counter']:
        """Iterate over index file, generate MotifPileups and pass to Visitors

        Parameters
        ----------
        index_file
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

        visitors: Dict[str, Visitor] = self._get_visitors()

        motif_pileup_generator = stepwise_pileup_generator(
            index_file=index_file,
            bam_path=self.bam_path
        )

        for motif_pileup in motif_pileup_generator:
            for curr_visitor in visitors.values():
                curr_visitor.process(motif_pileup)

        counters: Dict = _filter_for_counters(visitors)

        return counters

    @abstractmethod
    def _get_visitors(self) -> Dict[str, Visitor]:
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

    Implementation notes
    --------------------
    - overlap handling is not required because M-bias stats are stratified
    by sequencing strand, and thus automatically also by mate.
    - trimming is not required because we want to collect stats even in regions
    where we know that bias exists (for verification, if nothing else)

    Possible improvements
    ---------------------
    In the future, this class may be able sample the index positions to reduce
    the number of computations.
    """

    def _get_visitors(self) -> Dict[str, Visitor]:
        return OrderedDict(
            mbias_counts=MbiasCounter(self.config),
        )


class QcAndMethCallingRun(PileupRun):
    """Methylation calling with QC filtering and stats collection"""

    def _get_visitors(self) -> Dict[str, Visitor]:
        return OrderedDict(
            # trimmer=Trimmer(self.config, self.cutting_sites),
            # ol_handler=OverlapHandler(self.config),
            # mcall_writer=BedMcallWriter(self.config),
            # beta_counter=StratifiedBetaValueCounter(self.config),
            # coverage_counter=CoverageCounter(self.config,
            #                                  trimming_status='adjusted'),
        )


def _filter_for_counters(visitors) -> Dict:
    """Given a dict of Visitors, return a dict of Counters"""
    counters = {visitor_name: visitor_obj
                for visitor_name, visitor_obj in visitors.items()
                if isinstance(visitor_obj, Counter)}
    # Note: must return empty dict when no counters present
    return counters
