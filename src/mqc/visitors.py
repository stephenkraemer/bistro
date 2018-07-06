"""Visitor base classes"""

import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import os.path as op
from os import makedirs

from abc import ABCMeta, abstractmethod
from typing import List, Union

from mqc.pileup.pileup import MotifPileup
from mqc.utils import convert_array_to_df



class Visitor(metaclass=ABCMeta):
    """Just used for type annotations"""
    @abstractmethod
    def process(self, motif_pileup: MotifPileup) -> None:
        pass


class Counter(Visitor, metaclass=ABCMeta):
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
                 counter_array: np.ndarray,
                 save_stem: str) -> None:
        """ Initialization of attributes common to all Counters

        Every subclass should use its init function to:
        - initialize its counter array
        - call super().__init__()
        - save all config variables required for further use as attributes

        """

        if not isinstance(dim_names, list):
            raise TypeError('Counter expects list of dimension names')

        self.dim_names = dim_names
        self.dim_levels = dim_levels
        self.counter_array = counter_array
        self._counter_dataframe: pd.DataFrame = pd.DataFrame()

        self.save_stem = save_stem


    @abstractmethod
    def process(self, motif_pileup: MotifPileup) -> None:
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

    def get_dataframe(self) -> pd.DataFrame:
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

    def _compute_dataframe(self) -> pd.DataFrame:
        return convert_array_to_df(arr=self.counter_array,
                                   dim_levels=self.dim_levels,
                                   dim_names=self.dim_names,
                                   value_column_name='counts')

    def save_dataframe(self) -> None:
        """
        Saves the counter's dataframe as tsv and pickle to the stem specified in __init__.

        """
        makedirs(op.dirname(self.save_stem), exist_ok=True, mode=0o770)

        self.get_dataframe().to_pickle(self.save_stem+'.p')
        self.get_dataframe().reset_index().to_csv(self.save_stem+'.tsv', sep='\t', header=True, index=False)

    def save(self) -> None:
        """Save Counter as pickle"""
        out_path = Path(self.save_stem).with_suffix('.p')
        out_path.parent.mkdir(parents=True, exist_ok=True, mode=0o770)
        with out_path.open('wb') as fout:
            pickle.dump(self, fout)

