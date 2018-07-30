import itertools
import json
from collections import Hashable, OrderedDict
from copy import deepcopy
from pathlib import Path

import magic
import gzip

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from typing import (
    List, Set, Union, Sequence, Any, Dict, Optional, IO)


def convert_array_to_df(arr: np.ndarray,
                        dim_levels: Sequence[Sequence[Any]], dim_names: List[str],
                        value_column_name: str):
    """ Convert ndarray to DataFrame

    Parameters
    ----------
    arr: np.ndarray
    dim_levels:
        One sequence per index level, detailing the unique index levels
        *The index levels must be given in sorted order*
    dim_names
    value_column_name

    Returns
    -------
    pd.DataFrame
        One index level per dimension, one column with counts, named
        value_column_name

    Typical choices for the dim_level iterables are pd.Categorical, List[str],
    range, and List[pd.Interval]

    *Note* that this does not sort the index. If you gave all dim_levels in sorted
    order the index will however be in sorted form automatically (although this
    is not marked in the DataFrame)
    """
    if arr.size == 0:
        raise ValueError('NDarray is empty')

    # Create MultiIndex for multidimensional counter, and Index for 1D counter
    if arr.ndim > 1:
        midx = pd.MultiIndex.from_product(dim_levels, names=dim_names)
        arr.shape = (arr.size,)
        df = pd.DataFrame(arr, index=midx, columns=[value_column_name])
    else:
        idx = pd.Index(dim_levels[0], name=dim_names[0])
        df = pd.DataFrame(arr, index=idx, columns=[value_column_name])

    return df


def _create_dataframe_rows(arr, index_level_tuples):
    rows = []
    for idx_combi_with_levels in itertools.product(*index_level_tuples):
        # must be tuple for np array indexing
        idx_tuple = tuple(idx for idx, level in idx_combi_with_levels)
        levels = [level for idx, level in idx_combi_with_levels]
        curr_row = levels + [arr[idx_tuple]]
        rows.append(curr_row)
    return rows


def open_gzip_or_plain_file(filepath, mode='rt'):
    """Recognizes and opens plain and gzip text files

    Throws informative exceptions if
    - file does not exist
    - file is empty

    Uses python-magic to identify the file type; python-magic is a wrapper
    for libmagic, so this should be quite robust.

    Args:
        filepath (str): relative or absolute path to file
        mode (str): Defaults to 'rt'

    Returns:
        file-like object (text file)
    """

    # Note: There is a race condition here, as the file may
    # be deleted between those checks and the opening of the file
    # Improve or don't use where such cases are a problem

    # These checks work even if there are no access permissions to the file


    # TODO-low-prio: What happens if there are no access permissions to the dir
    #      containing the file? Will these checks fail with cryptic errors?
    if not os.path.exists(filepath):
        raise FileNotFoundError(f'File {filepath} not found')
    if os.path.getsize(filepath) == 0:
        raise IOError(f'File {filepath} is empty')

    try:
        filetype_mime = magic.from_file(filepath, mime=True)
    except OSError:
        raise OSError(f"Can't open {filepath}")

    # Need to account for different gzip media types
    # see: https://tools.ietf.org/id/draft-levine-application-gzip-03.html
    gzip_media_types = ['application/gzip', 'application/gzip-compressed',
                        'application/gzipped', 'application/x-gzip',
                        'application/x-gzip-compressed', 'gzip/document']
    if filetype_mime in gzip_media_types:
        try:
            fobj = gzip.open(filepath, mode)
        except OSError:
            raise OSError(f"Can't open {filepath}")
    else:
        try:
            fobj = open(filepath, mode)
        except OSError:
            raise OSError(f"Can't open {filepath}")

    return fobj


def get_resource_abspath(basename):
    return resource_filename('mqc.resources', basename)


def sel_expand(template, **kwargs):
    input_is_path = False
    if isinstance(template, Path):
        template = str(template)
        input_is_path = True
    fields = kwargs.keys()
    values = [kwargs[f] for f in fields]
    values = [[val] if isinstance(val, (int, str)) else val
              for val in values]
    value_combinations = itertools.product(*values)

    def get_expanded_template(template, fields, comb):
        for field, value in zip(fields, comb):
            template = template.replace('{' + field + '}', value)
        return template

    res = [get_expanded_template(template, fields, comb) for comb in value_combinations]

    if input_is_path:
        res = [Path(elem) for elem in res]

    if len(res) == 1:
        return res[0]
    else:
        return res


def update_nested_dict(base_dict, custom_dict, allow_new_keys=False) -> dict:
    """Update nested dict

    Parameters
    ----------
    base_dict
        Base dict acts as template. It specifies the allowed type for
        existing keys. If allow_new_keys=False, it also specifies the set
        of allowed keys
    custom_dict
        Overwrites base dict values

    Returns
    -------
    dict
        Merging result as new dict (deepcopy)

    Raises
    ------
    TypeError
        if type of the update value does not match the original value
    KeyError
        if custom_dict contains a key not present in base_dict (also
        within nested dicts) (if allow_new_keys=False)
    """

    base_dict = deepcopy(base_dict)

    for key, custom_value in custom_dict.items():

        default_value = base_dict.get(key)
        if default_value is None:
            if allow_new_keys:
                base_dict[key] = custom_value
                continue
            else:
                raise KeyError(f"Custom dict contains unknown key '{key}'")

        if type(custom_value) != type(default_value):
            raise TypeError(f"Can't overwrite type {type(default_value)}"
                            f" with type {type(custom_value)} for"
                            f"variable '{key}'")

        if isinstance(custom_value, dict):
            recursively_updated_dict = update_nested_dict(
                default_value, custom_value)
            base_dict[key] = recursively_updated_dict
        else:  # key points to 'scalar' value, not mapping
            base_dict[key] = custom_dict[key]

    return base_dict


# noinspection PyPep8Naming
def NamedIndexSlice(**kwargs) -> tuple:
    """Kwargs based index slice: name the index levels and specify the slices

    Args:
        **kwargs: level to slice object mapping

    Returns:
        function accepting a dataframe as input and returning an indexing
        tuple, for use in .loc indexing

    Examples:
        nidxs = NamedIndexSlice
        df.loc[nidxs(level1='a', level2=[1, 2], level5=[True, True, False]), :]

    Note:
        This function uses the ability of the .loc indexing to accept a
        callable. This callable is called with the dataframe as only
        argument and must return an object that can be evaluated by loc.

        This function returns a function which will produce a tuple with
        correct assignment of slice object to level, inserting slice(None)
        where necessary.

        Consecutive trailing slice(None) objects are removed from the
        indexing tuple. This is done to deal with issue
        https://github.com/pandas-dev/pandas/issues/18631 which essentially
        means: indexing with one or more scalars, e.g. ('a', ) or ('a', 'b')
        -> index levels are dropped. Indexing with one or more scalars
        followed by slice(None) or a list of values -> the index levels are
        not dropped. By leaving out trailing slice(None) objects index
        levels are dropped for scalar only indexing, which is the same
        behavior one would have achieved with standard loc indexing (where
        usually trailing slice(None) or : specs (in IndexSlice) are omitted
    """

    def fn(df):
        slicing_list = [kwargs.get(index_name, slice(None))
                        for index_name in df.index.names]
        for i in reversed(range(len(slicing_list))):
            if slicing_list[i] == slice(None):
                slicing_list.pop(i)
            else:
                break

        return tuple(slicing_list)

    return fn


def assert_match_between_variables_and_index_levels(
        variables: Union[List, Set], index_level_names: np.ndarray):
    """Check that a list or set of strings matches index levels

    Raises: ValueError

    The variables do not have to be ordered
    """

    variables = set(variables)
    index_level_names = set(index_level_names)

    if variables & set(index_level_names) != variables:
        unmatched_variables = variables - set(index_level_names)
        raise ValueError(
            'The following plotting variables are not in the index: ',
            unmatched_variables)

def subset_dict(key_or_keys: Union[Hashable, Sequence[Hashable]], d) -> dict:
    """Subset dict by one or more keys

    Args:
        key_or_keys: All hashable objects are interpreted as
            'one key'. This means that a tuple is *a single key*.
            Several keys should be given as sequence of hashables.
        d: the input dictionary

    Returns:
        A new dictionary (currently *not* a deepcopy), containing
        only the keys specified by key_or_keys
    """
    try:
        if isinstance(key_or_keys, Hashable):
            return {key_or_keys: d[key_or_keys]}
        elif isinstance(key_or_keys, Sequence):
            return {k: d[k] for k in key_or_keys}
        else:
            raise TypeError(f'Subset_dict cant handle input type {type(key_or_keys)}')
    except KeyError:
        raise KeyError("Error while subsetting dict."
                       f" Can't subset dict with keys {key_or_keys},"
                       " because at least one key is missing")

def hash_dict(d: dict) -> int:
    """Hash a dict: all keys must be str

    *ONLY USE WHEN ALL KEYS ARE STRING*

    Args:
        d: may be nested

    Returns:
        hash value

    Notes:
        This works by converting the dict to a str representation, and the
        hashing the str. To be independent of the key order in the dict,
        the  dicts must be sorted reproducibly, also in a nested dict. I
        use the json.dumps(sort_keys=True) functionality. I am not fully
        sure whether this is robust enough. Test for every use case.
        See also:
        https://stackoverflow.com/questions/5884066/hashing-a-dictionary
    """
    # TODO: verify that nested sorting is robust
    # TODO: verify that all keys are str
    return hash(json.dumps(d, sort_keys=True))

import os

class TmpChdir:
    """Context manager to step into directory temporarily"""
    def __init__(self, path):
        self.old_dir = os.getcwd()
        self.new_dir = path

    def __enter__(self):
        os.chdir(self.new_dir)

    def __exit__(self, *args):
        os.chdir(self.old_dir)


def csv_file_gen(file_obj: IO[Any], fieldnames: List[str],
                 field_type_dict=Optional[Dict[str, type]], sep: str ='\t') \
        -> Dict[str, Union[float, int, str]]:
    """Return field dict per line in file

    Args:
        field_type_dict: All field types are string by default.
            Optionally, types for one or more fields may be
            defined in the field_type_dict.

    Returns:
        Dict mapping field names to their values, with appropriate
        types. The items are ordered according to the field order.
        This relies on python 3.7+ dicts, and does not use OrderedDict.

    The input lines are stripped before processing.
    """

    for line in file_obj:
        fields_dict_ordered = dict(zip(fieldnames, line.strip().split(sep)))

        if field_type_dict is not None:
            for field_name, field_type_obj in field_type_dict.items():
                fields_dict_ordered[field_name] = field_type_obj(
                        fields_dict_ordered[field_name])

        yield fields_dict_ordered
