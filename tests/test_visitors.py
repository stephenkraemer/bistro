import numpy as np
import pandas as pd
import tempfile
import os.path as op

from mqc.visitors import Counter
from mqc.pileup.pileup import MotifPileup
import pytest


class CounterStub(Counter):
    def __init__(self, dim_levels, dim_names, arr, save_stem=None):
        super().__init__(dim_levels=dim_levels,
                         dim_names=dim_names,
                         counter_array=arr,
                         save_stem=save_stem)

    def process(self, motif_pileup: MotifPileup):
        pass


def df_not_correct_message(computed_df, expected_df):
    return (f"Computed dataframe is wrong\n"
            f"Expected dataframe: \n{repr(expected_df)}\n"
            f"Computed dataframe: \n{repr(computed_df)}\n")


def test_dataframe_saving():
    with tempfile.TemporaryDirectory() as tmpdir:
        df_file_stem = 'df_file'

        arr = np.array([[1, 2],
                        [3, 4]])
        dim_names = ['str_index_var', 'int_index']
        dim_levels = ['b a'.split(), range(2, 4)]

        counter_stub = CounterStub(dim_names=dim_names,
                                   dim_levels=dim_levels,
                                   arr=arr,
                                   save_stem=op.join(tmpdir, df_file_stem))
        counter_stub.save_dataframe()

        assert op.exists(op.join(tmpdir, df_file_stem + '.p')) and op.exists(
            op.join(tmpdir, df_file_stem + '.tsv'))


class TestArrayToDataFrameConversion:

    def test_converts_onedim_arr_to_df(self):
        arr = np.array([1, 2])
        dim_names = ['str_index_var']
        dim_levels = ['b a'.split()]

        expected_rows = [['b', 1], ['a', 2]]
        expected_df = (pd.DataFrame(expected_rows,
                                    columns=dim_names + ['counts'])
                       .set_index(dim_names))

        counter_stub = CounterStub(dim_names=dim_names, dim_levels=dim_levels,
                                   arr=arr)
        computed_df = counter_stub.get_dataframe()

        assert expected_df.equals(computed_df), df_not_correct_message(
            computed_df, expected_df)
        assert expected_df.index.equals(computed_df.index)


    def test_converts_multidim_arr_to_df(self):
        arr = np.array([
            [[1, 2],
             [3, 4]],
            [[5, 6],
             [7, 8]]])

        dim_names = ['str_index_var', 'int_index', 'cat_index']
        dim_levels = ['b a'.split(),
                      range(2, 4),
                      pd.Categorical(['e', 'd'],
                                     categories=['e', 'd'], ordered=True)]

        # Note that these rows are sorted, assumign that string dim levels are
        # given as categoricals in sorting order
        expected_rows = [['b', 2, 'e', 1],
                         ['b', 2, 'd', 2],
                         ['b', 3, 'e', 3],
                         ['b', 3, 'd', 4],
                         ['a', 2, 'e', 5],
                         ['a', 2, 'd', 6],
                         ['a', 3, 'e', 7],
                         ['a', 3, 'd', 8]]
        expected_df = (pd.DataFrame(expected_rows,
                                    columns=dim_names + ['counts'])
                       .set_index(dim_names))

        counter_stub = CounterStub(dim_names=dim_names, dim_levels=dim_levels,
                                   arr=arr)

        computed_df = counter_stub.get_dataframe()

        assert expected_df.equals(computed_df), df_not_correct_message(
            computed_df, expected_df)

        assert expected_df.index.equals(computed_df.index)

    def test_raises_if_dim_names_is_not_list_of_str(self):
        """Mainly relevant for 1D counters
        ... where name may be passed as str by accident"""
        with pytest.raises(TypeError):
            counter_stub = CounterStub(dim_levels=[[1, 2, 3]],
                                       dim_names='level_a',
                                       arr=np.array([1, 2]))
