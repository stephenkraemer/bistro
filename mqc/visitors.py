"""Visitor base classes"""

import itertools
import numpy as np
import pandas as pd

from abc import ABCMeta, abstractmethod
from typing import List, Union

from mqc.pileup.pileup import MotifPileup
from collections import Iterable


class Visitor(metaclass=ABCMeta):
    """Just used for type annotations"""
    @abstractmethod
    def process(self, motif_pileup: MotifPileup):
        pass


class Counter(metaclass=ABCMeta):
    """Store, update and retrieve count-based statistics

    Subclasses may for example store
    - beta value distribution data
    - coverage distribution data
    - M-bias stats
    - etc.

    Attributes
    ----------
        dim_names:
            names of counter dimensions, will be used as column names
            in counter dataframe
        dim_levels:
            List[Union[list, tuple, range]] containing
            levels for all dimensions.
        counter_array:
            Numpy array to hold counts
    """

    def __init__(self,
                 dim_names: List[str],
                 dim_levels: List[List[Union[str, int]]],
                 counter_array: np.ndarray):
        """ Initialization of attributes common to all Counters

        Every subclass should use its init function to:
        - initialize its counter array
        - call super.__init__()
        - save all config variables required for further use as attributes

        """

        if not isinstance(dim_names, list):
            raise TypeError('Counter expects list of dimension names')

        if not all([isinstance(elem, (list, tuple, range))
                    for elem in dim_levels]):
            raise TypeError('Counter expects dimension levels argument'
                            'to be list of list/tuple/range, with one list per'
                            'dimension')

        self.dim_names = dim_names
        self.dim_levels = dim_levels
        self.counter_array = counter_array
        self._counter_dataframe: pd.DataFrame = pd.DataFrame()

    @abstractmethod
    def process(self, motif_pileup: MotifPileup):
        """Update counter inplace based on information in MotifPilepup

        Notes
        -----
        The process method will in general follow this algorithm:
        1. retrieve attributes of the MotifPileup or the contained PileupReads
        2. transform attributes into integer indices, e.g.

            - introduce upper bound on values
            - bin values
            - map values (e.g. string values or flag values) to integer indices

        3. Increment counter at element specified by
           the obtained integer indices
        """
        pass

    def get_dataframe(self):
        """Get counts in dataframe format

        The dataframe is cached after the first computation, so repeated
        calls are fine.

        Returns
        -------
        pd.DataFrame

            Dataframe with one column per array dimension, named after
            :attr:`~mqc.Counter.dim_names`, and an additional column
            for the counts, named 'Counts'
        """
        if self._counter_dataframe.empty:
            self._counter_dataframe = self._compute_dataframe()
        return self._counter_dataframe

    def _compute_dataframe(self):
        indices_per_dim = [np.arange(x) for x in self.counter_array.shape]
        index_level_tuples = [
            list(zip(curr_indices, curr_levels))
            for curr_indices, curr_levels
            in zip(indices_per_dim, self.dim_levels)]

        rows = self._create_dataframe_rows(index_level_tuples)
        column_names = self.dim_names + ['Counts']
        df = (pd.DataFrame(rows, columns=column_names)
              .sort_values(by=self.dim_names))
        df.index = range(0, df.shape[0])
        return df

    def _create_dataframe_rows(self, index_level_tuples):
        rows = []
        for idx_combi_with_levels in itertools.product(*index_level_tuples):
            # must be tuple for np array indexing
            idx_tuple = tuple(idx for idx, level in idx_combi_with_levels)
            levels = [level for idx, level in idx_combi_with_levels]
            curr_row = levels + [self.counter_array[idx_tuple]]
            rows.append(curr_row)
        return rows
