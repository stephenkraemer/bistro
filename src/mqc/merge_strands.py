"""Merge strands in symmetric methylation calling motifs

Deals with merging of methylation calls (BED and stratified-BED format)
and of index files.

Contains the functions called by meth_calls merge and index merge
"""

import gzip
import logging
from itertools import islice, chain, count, starmap
from toolz import partition
from operator import itemgetter
from typing import Dict, NamedTuple, List, Iterable, Any, Union, IO, Iterator, Tuple

from mqc.utils import csv_file_gen
from mqc.index import possible_index_field_names

class BedMethCallsMerger:
    f"""Merge strand-resolved CpG methylation stats

    Args:
        strand_resolved_meth_calls: strand-resolved methylation calls created
            with bistro. bed or stratified_bed format.
        merged_meth_calls: file path for merged calls

    The following index columns are recognized:
    {possible_index_field_names}

    All other columns are considered data columns and will be aggregated
    per motif. All field names containing the substrings 'n_meth' or
    'n_total' are considered ints, all field_names containing the substring
    'beta_value' are considered floats.

    The interval boundaries are merged, all other fields are taken
    from the plus strand.
    """

    def __init__(self, strand_resolved_meth_calls: str, merged_meth_calls: str,
                 verbose=False) -> None:
        self.strand_resolved_meth_calls = strand_resolved_meth_calls
        self.merged_meth_calls = merged_meth_calls
        self.logger = logging.Logger(__name__ + '.BedMethCallsMerger')
        ch = logging.StreamHandler()
        if verbose:
            formatter = logging.Formatter(
                    '%(asctime)s - %(name)s - %(message)s')
            ch.setFormatter(formatter)
            ch.setLevel(logging.DEBUG)
        else:
            formatter = logging.Formatter('%(message)s')
            ch.setFormatter(formatter)
            ch.setLevel(logging.INFO)
        self.logger.addHandler(ch)

    def run(self) -> None:
        self.logger.info('Merging strand-resolved methylation calls: starting now...')
        with gzip.open(self.strand_resolved_meth_calls, 'rt') as fin, \
                gzip.open(self.merged_meth_calls, 'wt') as fout:

            header = fin.readline()

            field_names = header.rstrip().split('\t')
            field_type_dict = self._data_column_to_type_mapping(field_names)
            self.logger.debug(f'The following column-to-type mappings are used: {field_type_dict}')
            index_cols = [field_name for field_name in field_names
                          if field_name in possible_index_field_names]
            self.logger.debug(f'The following columns were identified as index columns: {index_cols}')
            n_index_cols = len(index_cols)
            # assert that index columns are at the beginning of the file
            assert field_names[0:len(index_cols)] == index_cols
            assert set(field_type_dict.keys()) == set(field_names[n_index_cols:])

            fout.write(header)
            field_dict_gen = csv_file_gen(file_obj=fin, fieldnames=field_names,
                                          field_type_dict=field_type_dict, sep='\t')
            for line in self._create_merged_calls_gen(field_dict_gen, n_index_cols):
                fout.write(line)
        self.logger.info('Successfully finished merging.')

    @staticmethod
    def _data_column_to_type_mapping(field_names: List[str]) -> Dict[str, type]:
        field_type_dict: Dict[str, type] = {}
        int_field_substrs = ['n_meth', 'n_unmeth', 'n_total']
        float_field_substrs = ['beta_value']
        for field_name in field_names:
            if any([int_field_substr in field_name
                    for int_field_substr in int_field_substrs]):
                field_type_dict[field_name] = int
            elif any([float_field_substr in field_name
                      for float_field_substr in float_field_substrs]):
                field_type_dict[field_name] = float
        return field_type_dict

    def _create_merged_calls_gen(self, field_dict_gen: Iterator[Dict[str, Any]],
                                 n_index_cols: int) -> Iterator[str]:
        """Consume lines in pairs and merge them

        Args:
            field_dict_gen: Generator returning lines of the BED file
                as dictionary with correct types
            n_index_cols: Number of index columns in BED file (as
                opposed to columns containing data)

        Returns:
            merged line (with newline already added)

        Implementation notes:
        - Expects the following column names: #chrom, start, end
        - ZeroDivision is replaced by nan
        """
        for i in count():
            if i % 300_000 == 0:
                self.logger.debug(f'At position {i}')
            try:
                field_dict1_ordered = next(field_dict_gen)
                field_dict2_ordered = next(field_dict_gen)
            except StopIteration:
                break

            def merge_stats(stats1: Tuple, stats2: Tuple) -> Tuple[str, str, str]:
                n_meth = stats1[1] + stats2[1]
                n_total = stats1[2] + stats2[2]
                if n_total == 0:
                    beta_value_str = 'nan'
                else:
                    beta_value_str = f'{n_meth / n_total:.8f}'

                return beta_value_str, str(n_meth), str(n_total)

            # meth. stats for bed and stratified bed are always beta_value
            # n_meth n_total for different strata. Create chunks of each 3 subsequent
            # stats for line1 and line2, then merge them
            merged_stats = (starmap(merge_stats,
                                    zip(partition(3, islice(field_dict1_ordered.values(), n_index_cols, None)),
                                        partition(3, islice(field_dict2_ordered.values(), n_index_cols, None)))))

            line = '\t'.join(chain(
                    (field_dict1_ordered['#chrom'], field_dict1_ordered['start'], field_dict2_ordered['end']),
                    islice(field_dict1_ordered.values(), 3, n_index_cols),
                    *merged_stats
            ))

            yield line + '\n'
