import numpy as np
import pandas as pd

from mqc.visitors import Counter
from mqc.pileup.pileup import MotifPileup
import pytest


class CounterStub(Counter):
    def __init__(self, dim_levels, dim_names, arr):
        super().__init__(dim_levels=dim_levels,
                         dim_names=dim_names,
                         counter_array=arr)

    def process(self, motif_pileup: MotifPileup):
        pass


def df_not_correct_message(computed_df, expected_df):
    return (f"Computed dataframe is wrong\n"
            f"Expected dataframe: \n{repr(computed_df)}\n"
            f"Computed dataframe: \n{repr(expected_df)}\n")


class TestArrayToDataFrameConversion:
    def test_converts_multidim_arr_to_df_with_sorted_index_variables(self):
        arr = np.array([[1, 2],
                        [3, 4]])
        dim_names = ['str_index_var', 'int_index']
        dim_levels = ['b a'.split(), range(2, 4)]

        # Note that these rows are sorted
        expected_rows = [['a', 2, 3], ['a', 3, 4], ['b', 2, 1], ['b', 3, 2]]
        expected_df = pd.DataFrame(expected_rows,
                                   columns=dim_names + ['Counts'])
        expected_df.index = range(0, 4)

        counter_stub = CounterStub(dim_names=dim_names, dim_levels=dim_levels,
                                   arr=arr)

        computed_df = counter_stub.get_dataframe()

        assert expected_df.equals(computed_df), df_not_correct_message(
            computed_df, expected_df)

        assert expected_df.index.equals(computed_df.index)

    def test_converts_onedim_arr_to_df(self):
        arr = np.array([1, 2])
        dim_names = ['str_index_var']
        dim_levels = ['b a'.split()]

        expected_rows = [['a', 2], ['b', 1]]
        expected_df = pd.DataFrame(expected_rows,
                                   columns=dim_names + ['Counts'])
        expected_df.index = range(0, 2)

        counter_stub = CounterStub(dim_names=dim_names, dim_levels=dim_levels,
                                   arr=arr)
        computed_df = counter_stub.get_dataframe()

        assert expected_df.equals(computed_df), df_not_correct_message(
            computed_df, expected_df)
        assert expected_df.index.equals(computed_df.index)

    def test_raises_if_dim_levels_are_not_list_of_list(self):
        with pytest.raises(TypeError):
            counter_stub = CounterStub(dim_levels=[1, 2, 3],
                                       dim_names=['level_a'],
                                       arr=np.array([1, 2]))

    def test_raises_if_dim_names_is_not_list_of_str(self):
        """Mainly relevant for 1D counters
        ... where name may be passed as str by accident"""
        with pytest.raises(TypeError):
            counter_stub = CounterStub(dim_levels=[[1, 2, 3]],
                                       dim_names='level_a',
                                       arr=np.array([1, 2]))