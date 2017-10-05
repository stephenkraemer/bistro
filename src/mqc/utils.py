import itertools
import magic
import os
import gzip

import numpy as np
import pandas as pd
from pkg_resources import resource_filename


def convert_array_to_df(arr, dim_levels, dim_names, value_column_name):
    """ Convert ndarray to DataFrame

    Parameters
    ----------
    arr: np.ndarray
    dim_levels: List[Iterable]
        One iterable per index level, detailing the unique index levels
        *The index levels must be given in sorted order*
    dim_names: List[str]
    value_column_name: str

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
    midx = pd.MultiIndex.from_product(dim_levels, names=dim_names)
    arr.shape = (arr.size,)
    df = pd.DataFrame(arr, index=midx, columns=[value_column_name])
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
