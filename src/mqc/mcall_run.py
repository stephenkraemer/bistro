"""Commands processing a BAM file by iterating over MotifPileups

Provides PileupRun ABC as template for algorithms that parse a BAM
file by processing MotifPileups iteratively. To create a valid Runner
by subclassing PileupRun, the abstract method _get_visitors must be
implemented. This method determines the order of application of different
Visitors to the MotifPileup, which defines the processing required by
the user. The PileupRun ABC provides the logic for running the processing
in parallel across multiple chromosomes.

You can create own Runners by subclassing PileupRun in your own modules.

This module provides the base runners shipped with mqc, currently:
- MbiasDeterminationRun :: collects M-bias stats
- QcAndMethCallingRun :: performs methylation calling with various
      output formats and QC filters

SNP calling will be added, and functionality to gather information in the
biscuit epiread format.
"""

import json
import multiprocessing as mp
import sys
from abc import abstractmethod, ABCMeta
from collections import OrderedDict
from copy import deepcopy
from pathlib import Path
from typing import Dict, List, Any, Optional

import pandas as pd
import pysam

from mqc.index import IndexFile
from mqc.mbias import MbiasCounter, CuttingSitesReplacementClass
from mqc.mcaller import MethCaller, StratifiedMethCaller
from mqc.overlap import OverlapHandler
from mqc.pileup.pileup import stepwise_pileup_generator
from mqc.qc_filters import PhredFilter, MapqFilter
from mqc.trimming import Trimmer
from mqc.visitors import Counter, Visitor
from mqc.writers import BedWriter, BismarkWriter, StratifiedBedWriter


# from mqc.mbias import MbiasData, AdjustedMbiasCuttingSites, MbiasCounter
# from mqc.beta_values import StratifiedBetaValueCounter
# from mqc.coverage import CoverageCounter

ConfigDict = Dict[str, Any]


def collect_mbias_stats(config: ConfigDict) -> None:
    """Runner function for stats command

    Collects M-bias stats as M-bias Counter pickle

    Will be changed to M-bias dataframe output in the future, the
    current output just makes development on M-bias stats processing
    easier...
    """

    mbias_run = MbiasDeterminationRun(config, cutting_sites=None)
    mbias_run.run_parallel()
    for counter in mbias_run.summed_up_counters.values():
        counter.save()


def run_mcalling(config: ConfigDict) -> None:
    """Methylation calling runner function"""

    trimm_command, trimm_param = config['run']['trimming'].split('::')
    max_read_length = config['run']['max_read_length']

    if trimm_command == 'read_end':
        # TODO: fill
        raise NotImplementedError

    elif trimm_command == 'frag_end':
        param_err_message = (f'Cannot interprete the parameter {trimm_param} '
                             f'for command {trimm_command}')
        try:
            trimm_param_dict = json.loads(trimm_param)
            assert set(trimm_param_dict.keys()) == {*'c_bc c_bc_rv w_bc w_bc_rv'.split()}, \
                'Fragment end-based trimming specification does not have all BS-Seq strands'
        except AttributeError:
            print(param_err_message, file=sys.stderr)
            raise
        if not isinstance(trimm_param_dict['c_bc'], list):
            raise ValueError(param_err_message)

        cutting_sites = CuttingSitesReplacementClass.from_rel_to_frag_end_cutting_sites(
            cut_site_spec=trimm_param_dict, max_read_length=max_read_length)

    elif trimm_command == 'cutting_sites':

        assert Path(trimm_param).exists(), f'Can not find cutting sites df {trimm_param}'
        suffix = Path(trimm_param).suffix
        assert suffix in ['.p', '.feather', '.tsv'], \
            f'Cutting sites dataframe file does not have a supported format ({suffix})'

        if suffix == '.p':
            cutting_sites_df = pd.read_pickle(trimm_param)
        elif suffix == '.feather':
            # TODO-important: fill
            raise NotImplementedError
        else:  # tsv
            # TODO-important: fill
            raise NotImplementedError

        cutting_sites = CuttingSitesReplacementClass(cutting_sites_df=cutting_sites_df,
                                                     max_read_length=max_read_length)

    # noinspection PyUnboundLocalVariable
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

    def __init__(self, config: ConfigDict,
                 cutting_sites: Optional[CuttingSitesReplacementClass] = None) -> None:

        self.summed_up_counters: Dict[str, Counter] = {}

        # All config vars used directly by methods of this class
        self.n_cores = config['run']['cores']
        self.index_files = config['run']['index_files']
        self.bam_path = config['run']['bam']
        self.cutting_sites = cutting_sites

        # Save whole config dict, because it is used
        # by the visitor class constructors
        self.config = config

    def run_parallel(self) -> None:
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

    def sum_up_counters(
            self, counter_dicts_per_idx_file: List[Dict[str, Counter]]) -> None:
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

    def _single_run_over_index_file(self, index_file_path: str) -> Dict[str, Counter]:
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
    def _get_visitors(self, chrom: str) -> Dict[str, Visitor]:
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

    def _get_visitors(self, chrom: str) -> Dict[str, Visitor]:
        # noinspection PyTypeChecker
        return OrderedDict(
            mbias_counter=MbiasCounter(self.config),
        )


class QcAndMethCallingRun(PileupRun):
    """Methylation calling with QC filtering and stats collection

    Notes:
    - phred filtering is usually done after overlap handling,
      because overlap handling can be used to adjust phred scores.
      However, this is currently not implemented in the default
      overlap handling.
    """

    def _get_visitors(self, chrom: str) -> Dict[str, Visitor]:
        """Assemble visitors for methylation calling run

        First, apply different QC-related visitors which tag the reads
        which should be discarded.

        Then add the McallWriters as specified in the output_formats
        config option. For most current and future writers,
        it will be necessary to add a Visitor preparing the
        data to be written before the McallWriter. E.g. bed ->
        MethCaller, stratified_bed -> StratifiedMethCaller, or in
        the future VCF -> SnpCaller, epiread -> EpireadAssembler. Only
        the BismarkWriter does not depend on any prior Visitor.

        Writing both bed and stratified bed output is currently
        not implemented, because the stratified bed also contains
        all information from the standard bed, just with some
        additional informations and a different header / column order.

        This would be easy to implement if requested.
        """

        if self.cutting_sites is not None:
            visitors = OrderedDict(
                mapq_filter=MapqFilter(self.config),
                trimmer=Trimmer(cutting_sites=self.cutting_sites),
                overlap_handler=OverlapHandler(),
                phred_filter=PhredFilter(self.config),
            )
        else:
            raise TypeError('Methylation calling expects a valid CuttingSites instance')

        output_formats_list: List[str] = self.config['run']['output_formats']

        if 'stratified_bed' in output_formats_list:
            visitors['meth_caller'] = StratifiedMethCaller()
            visitors['stratified_bed_writer'] = StratifiedBedWriter(
                    calls_by_chrom_motif_fp=self.config['paths']['strat_bed_calls_by_chrom_motif'],
                    motifs=self.config['run']['motifs'],
                    chrom=chrom).setup()
            if 'bed' in output_formats_list:
                # the StratifiedBedWriter provides the necessary MotifPileup attributes
                print('WARNING: creating BED and stratified BED output at the same time.'
                      ' These files are redundant - did you specify both options together by mistake?')
                visitors['bed_writer'] = BedWriter(self.config, chrom=chrom)
        elif 'bed' in output_formats_list:
            visitors['meth_caller'] = MethCaller()
            visitors['bed_writer'] = BedWriter(self.config, chrom=chrom)

        if 'bismark' in output_formats_list:
            visitors['bismark_writer'] = BismarkWriter(
                calls_by_chrom_motif_fp=self.config['paths']['bismark_calls_by_chrom_motif'],
                motifs=self.config['run']['motifs'],
                chrom=chrom
            ).setup()

        return visitors


def _filter_for_counters(visitors: Dict[str, Visitor]) -> Dict[str, Counter]:
    """Filter dict of visitors, return only the counters"""
    counters = {visitor_name: visitor_obj
                for visitor_name, visitor_obj in visitors.items()
                if isinstance(visitor_obj, Counter)}
    # Note: must return empty dict when no counters present
    return counters
