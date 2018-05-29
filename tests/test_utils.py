import json
import os
from collections import OrderedDict

from mqc.utils import (
    open_gzip_or_plain_file,
    sel_expand,
    update_nested_dict,
    NamedIndexSlice,
    assert_match_between_variables_and_index_levels,
    subset_dict,
    hash_dict, TmpChdir)
import pytest
import gzip
import pandas as pd
import numpy as np

class TestOpenGzipOrPlainFile:
    def test_fails_when_file_not_found(self):
        with pytest.raises(FileNotFoundError):
            open_gzip_or_plain_file('/path/to/chrom.fa.gz')


    def test_fails_when_file_empty(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        tmp_path.write('')
        error_message = f'File {tmp_path} is empty'
        with pytest.raises(IOError) as e:
            open_gzip_or_plain_file(str(tmp_path))
        assert str(e.value) == error_message


    def test_fails_with_no_access_permissions(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        tmp_path.write(content)
        tmp_path.chmod(0o000)
        error_message = f"Can't open {tmp_path}"
        with pytest.raises(OSError) as e:
            open_gzip_or_plain_file(str(tmp_path))
        assert str(e.value) == error_message


    def test_recognizes_plain_text_files(self, tmpdir):
        tmp_fa_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        tmp_fa_path.write(content)
        read_content = (open_gzip_or_plain_file(str(tmp_fa_path))
                        .read())
        assert content == read_content


    def test_recognizes_gzip_reads_as_unicode_str(self, tmpdir):
        tmp_path = tmpdir.join('chrom1.fa')
        content = 'ACGC\nACGC\n'
        with gzip.open(str(tmp_path), 'wt') as fobj:
            fobj.write(content)
        read_content = (open_gzip_or_plain_file(str(tmp_path))
                        .read())
        assert content == read_content


@pytest.mark.xfail(skip=True)
def test_sel_expand():
    raise NotImplemented

# Test update_nested_dict
# ==============================================================================

@pytest.fixture()
def base_dict():
    return {
        'key1': {
            'key1.1': 1.1,
            'key1.2': {
                'key1.3': 1.3
            }
        },
        'key2': {
            'key2.1': 2.1,
            'key2.2': 2.2,
        },
        'key3': 3.0,
        'key4': 4.0,
    }


@pytest.fixture()
def custom_dict():
    return {
        'key1': {
            'key1.2': {
                'key1.3': 1.32,
            }
        },
        'key2': {
            'key2.1': 2.12,
        },
        'key3': 3.2
    }


@pytest.fixture()
def expected_dict():
    return {
        'key1': {
            'key1.1': 1.1,
            'key1.2': {
                'key1.3': 1.32
            }
        },
        'key2': {
            'key2.1': 2.12,
            'key2.2': 2.2,
        },
        'key3': 3.2,
        'key4': 4.0,
    }


class TestUpdateNestedDict:

    def test_updates_nested_keys(
            self, base_dict, custom_dict, expected_dict):
        computed_dict = update_nested_dict(base_dict, custom_dict)
        assert computed_dict == expected_dict

    def test_allows_new_keys_if_new_keys_allowed(
            self, base_dict, custom_dict, expected_dict):
        custom_dict['key5'] = 5.0
        expected_dict['key5'] = 5.0
        computed_dict = update_nested_dict(
            base_dict, custom_dict, allow_new_keys=True)
        assert computed_dict == expected_dict

    def test_raises_keyerror_if_new_key_by_default(
            self, base_dict, custom_dict, expected_dict):
        custom_dict['key5'] = 5.0
        expected_dict['key5'] = 5.0
        with pytest.raises(KeyError):
            computed_dict = update_nested_dict(
                base_dict, custom_dict)

    def test_raises_typeerror_if_replacement_type_does_not_match(
            self, base_dict, custom_dict, expected_dict):
        custom_dict['key2']['key2.1'] = '2'
        with pytest.raises(TypeError):
            computed_dict = update_nested_dict(base_dict, custom_dict)
            assert computed_dict == expected_dict

    def test_returned_dict_is_deepcopy(
        self, base_dict):
        computed_dict = update_nested_dict(base_dict, {})
        assert id(computed_dict) != id(base_dict)
        assert id(computed_dict['key1']) != id(base_dict['key1'])


# Test named index slicing
# ==============================================================================

nidxs = NamedIndexSlice


@pytest.fixture()
def df_for_named_index_slice_tests():
    multi_idx = pd.MultiIndex.from_product(
        [[1, 2, 3],
         ['a', 'b', 'c'],
         [1, 2, 3],
         [1, 2, 3]],
        names='level0 level1 level2 level3'.split()
    )

    df = pd.DataFrame(np.random.randn(3**4, 2),
                      index=multi_idx,
                      columns=['col1', 'col2'],
                      )
    return df

named_index_slice_index_value_params = [
    OrderedDict(
        level1=['a', 'b'],
        level0=3,
        level2=np.concatenate(
            [np.tile([False], 9),
             np.tile([True], 3**4 - 9)]),
        level3=slice(1, 2)
    ),
    OrderedDict(
        level1='a',
        level0=3,
    ),
    OrderedDict(
        level0=3,
    ),
    OrderedDict(
        level1='a',
        level0=slice(None),
    ),
]


@pytest.mark.parametrize(
    'index_value_orddict', named_index_slice_index_value_params)
class TestNamedIndexSlice:

    def test_retrieves_correct_values(
            self, df_for_named_index_slice_tests,
            index_value_orddict):

        index_values = pd.Series(index_value_orddict)
        df = df_for_named_index_slice_tests

        nidxs_slice = df.loc[nidxs(**index_values), :]

        pd.testing.assert_frame_equal(
            nidxs_slice,
            df.loc[tuple(index_values.sort_index().values), :])


    def test_allows_assignment_to_correct_slice(
            self, df_for_named_index_slice_tests,
            index_value_orddict):

        index_values = pd.Series(index_value_orddict)
        df = df_for_named_index_slice_tests

        df.loc[nidxs(**index_values), :] = 10

        assert all(df.loc[tuple(index_values.sort_index().values), :] == 10)


def test_assert_match_between_variables_and_index_levels():
    midx = pd.MultiIndex.from_product([list('abc'), [1, 2, 3], list('cde')],
                                      names=['col0', 'col1', 'col2'])
    df = pd.DataFrame(index=midx)
    df['values'] = 0

    assert_match_between_variables_and_index_levels(
        variables=['col2', 'col1'], index_level_names=df.index.names
    )

    with pytest.raises(ValueError):
        assert_match_between_variables_and_index_levels(
            variables=['col2', 'COL0'], index_level_names=df.index.names)


# noinspection PyMethodMayBeStatic
class TestSubsetDict:
    @pytest.mark.parametrize('keys,correct_subset', [
        (['a'], {'a': 'a'}),
        (['b', 1], {'b': 'b', 1: 1}),
        (1, {1: 1}),
        ('b', {'b': 'b'}),
        (('z', 9), {('z', 9): 'z9'}),
    ])
    def test_subsets_dict_to_keys(self, keys, correct_subset):
        d = {'a': 'a', 'b': 'b', 1: 1, 2: 2, ('z', 9): 'z9'}
        assert subset_dict(keys, d) == correct_subset

    def raises_if_key_is_missing(self):
        d = {'a': 'a', 'b': 'b', 1: 1, 2: 2}
        with pytest.raises(KeyError):
             subset_dict(['e', 1], d)

    def raises_if_not_given_list_of_str(self):
        d = {'a': 'a', 'b': 'b', 1: 1, 2: 2}
        with pytest.raises(TypeError):
            subset_dict(TestNamedIndexSlice, d)


@pytest.mark.parametrize(
    'unsorted_d,sorted_d',
    [(
            OrderedDict([('b', 1),
                     ('a', OrderedDict([('d', 3), ('c', 2)]))]),
            OrderedDict([('a', OrderedDict([('c', 2), ('d', 3)])),
                         ('b', 1)]),
    )])
def test_hash_dict(unsorted_d, sorted_d):
    # Test that nested, unsorted dict has same hash as pre-sorted dict
    # to make sure that hash_dict correctly sorts the dict keys
    # when stringifying
    assert hash(json.dumps(sorted_d)) == hash_dict(unsorted_d)



def test_tempchdir(tmpdir):

    os.chdir(tmpdir)
    home_fp = os.path.expanduser('~')

    with TmpChdir(home_fp):
        res = 1 + 1
        assert os.getcwd() == home_fp
    assert os.getcwd() == str(tmpdir)

    with pytest.raises(KeyError):
        with TmpChdir(home_fp):
            res = 1 + 1
            assert os.getcwd() == home_fp
            raise KeyError

    assert os.getcwd() == str(tmpdir)

