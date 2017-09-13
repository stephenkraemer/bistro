"""Visitor base classes"""

import itertools
import numpy as np
import pandas as pd
import os.path as op
from os import makedirs

from abc import ABCMeta, abstractmethod
from typing import List, Union

from mqc.pileup.pileup import MotifPileup
from mqc.utils import convert_array_to_df
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
        - call super().__init__()
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
            for the counts, named 'counts'
        """
        if self._counter_dataframe.empty:
            self._counter_dataframe = self._compute_dataframe()
        return self._counter_dataframe

    def _compute_dataframe(self):
        return convert_array_to_df(arr=self.counter_array,
                                   dim_levels=self.dim_levels,
                                   dim_names=self.dim_names,
                                   value_column_name='counts')

    def save_dataframe(self, save_location: str):
        """
        Saves the counter's dataframe to a specified location.
        File type is determined by given path.
        :param save_location: Dataframe save location for pickle, tsv or csv file.
        """
        ft = save_location.split('.')[-1]

        makedirs(op.dirname(save_location), exist_ok=True, mode=0o777)

        if ft == 'p':
            self.get_dataframe().to_pickle(save_location)
        elif ft == 'csv':
            self.get_dataframe().reset_index().to_csv(save_location, header=True, index=False)
        elif ft == 'tsv':
            self.get_dataframe().reset_index().to_csv(save_location, sep='\t', header=True, index=False)
        else:
            print(f"Unsupported file type: \"{save_location}\"")

