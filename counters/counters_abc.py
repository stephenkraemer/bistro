from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import pandas as pd
import mqc
from typing import List, Tuple


class Counter(meta=ABCMeta):
    @property
    @abstractmethod
    def index_names(self):
        pass

    @property
    @abstractmethod
    def index_levels(self):
        pass

    def __init__(self,
                 array_shape: List, dtype,
                 dim_levels,
                 dim_names = None):
        """
        
        :param array_shape: 
        :param dtype: 
        :param labels: List[Union[List[str], int]],
                       one element per dimension, allowed values:
                       List[str]: one value per level, i.e. string labels for all levels
                       int: stands for: use integer labels, starting at int for this dimension
        """
        self.counter_array = self.initialize_counter_array()

    @abstractmethod
    def initialize_counter_array(self):
        """Make sure that appropriate dtype is selected"""
        pass

    def process(self, motif_pileup: mqc.MotifPileup) -> None:
        indices: List[Tuple[int]]
        indices = self.get_indices(motif_pileup)
        for idx_tuple in indices:
            self.counter_array[idx_tuple] += 1

    def counter_array_to_dataframe(self):
        "will use self.labels!"
        self.counter_dataframe = pd.DataFrame()

    @abstractmethod
    def get_indices(self, motif_pileup: mqc.MotifPileup) -> List[Tuple[int]]:
        pass
