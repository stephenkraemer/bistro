"""
input:

MotifPileupStub
---------------
idx_pos: IndexPositionStub
beta_value
n_meth
n_unmeth

IndexPositionStub
-----------------
chrom
start
motif
strand

output:
-------

FileObjectStub
--------------

line written to file

chrom start start+1 motif '.' strand beta n_meth n_unmeth
"""
import re

import pytest
import os.path as op

from collections import namedtuple
from unittest.mock import mock_open, MagicMock, call

from mqc.writers import BedWriter

TESTS_DIR = op.dirname(__file__)
DEFAULT_CONFIG_FILE = op.join(TESTS_DIR, '../config.default.toml')
SAMPLE_NAME = 'hsc_rep1'
SAMPLE_META = 'population=hsc,rep=1,model=blk6'

def test_bed_writer_computes_bed_lines_and_writes_to_one_file_per_motif(mocker):
    """
    input: MotifPileup
    - compute bed line
    side_effect: write BED line to the mcalls file for the given motif
    """

    index_pos_stub_cg = namedtuple('IndexPositionStub', 'chrom start motif strand end')(
        chrom = '1',
        start = 10,
        end = 11,
        motif = 'CG',
        strand = '+',
    )

    # TODO: switch to n_total
    motif_pileup_stub_cg = namedtuple('MotifPileupStub', 'idx_pos beta_value n_meth n_unmeth')(
        idx_pos = index_pos_stub_cg,
        beta_value = 0.6,
        n_meth = 3,
        n_unmeth = 2
    )

    index_pos_stub_chg = namedtuple('IndexPositionStub', 'chrom start motif strand end')(
        chrom = '1',
        start = 10,
        end = 11,
        motif = 'CHG',
        strand = '+',
    )

    # TODO: switch to n_total
    motif_pileup_stub_chg = namedtuple('MotifPileupStub', 'idx_pos beta_value n_meth n_unmeth')(
        idx_pos = index_pos_stub_chg,
        beta_value = 1,
        n_meth = 5,
        n_unmeth = 0
    )
    config_stub = {'paths': {'meth_calls_basepath': 'does/not/matter'},
                   'meta': {'name': 'hsc_1'},
                   'run': {'motifs': 'CG CHG CHH'.split()}}

    file_mocks = [MagicMock(), MagicMock(), MagicMock()]
    mock_open = MagicMock(side_effect = file_mocks)
    mocker.patch("mqc.writers.gzip.open", mock_open)


    writer = BedWriter(config_stub, chrom='1')
    writer.process(motif_pileup_stub_cg)
    writer.process(motif_pileup_stub_chg)

    expected_line_cg = "1 10 11 CG . + 0.6 3 2\n".replace(' ', '\t')
    expected_line_chg = "1 10 11 CHG . + 1 5 0\n".replace(' ', '\t')

    expected_header = '\t'.join(['#chrom', 'start', 'end',
                                 'motif', 'score', 'strand',
                                 'beta_value', 'n_meth', 'n_unmeth'])
    order_of_file_openening = [re.search(r'(CG|CHG|CHH)', x[0][0]).group(1)
                               for x in mock_open.call_args_list]
    file_mock_motif_mapping = {x: file_mocks[i]
                               for i, x in enumerate(order_of_file_openening)}
    assert file_mock_motif_mapping['CG'].write.call_args_list == [call(expected_header), call(expected_line_cg)]
    assert file_mock_motif_mapping['CHG'].write.call_args_list == [call(expected_header), call(expected_line_chg)]
