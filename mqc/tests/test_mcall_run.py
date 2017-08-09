import numpy as np
import pytest
from collections import defaultdict, OrderedDict
from unittest.mock import Mock, MagicMock, call

from mqc.index import IndexFile
from mqc.mcall_run import PileupRun
from mqc.visitors import Counter, Visitor


def get_counter_stub(arr, attr_str):
    return MagicMock(spec_set=['counter_array',
                               'config_attribute',
                               'process',
                               '__class__'],
                     counter_array=arr.copy(),
                     config_attribute=attr_str,
                     __class__=Counter)


class PileupRunStubWithVisitorsAndCounters(PileupRun):
    def _get_visitors(self, chrom):
        # empty array must match dimensions of arrays to be added
        empty_arr = np.array([[0, 0], [0, 0]])
        return OrderedDict(
            counter1=get_counter_stub(empty_arr, 'counter1'),
            counter2=get_counter_stub(empty_arr, 'counter2'),
            visitor=MagicMock(__class__=Visitor),
        )


class PileupRunStubNoCounters(PileupRun):
    def _get_visitors(self, chrom):
        # empty array must match dimensions of arrays to be added
        empty_arr = np.array([[0, 0], [0, 0]])
        return OrderedDict(
            visitor=MagicMock(__class__=Visitor)
        )


@pytest.fixture
def config():
    config = defaultdict(dict)
    config['run']['cores'] = 3
    index_file_mock = Mock(spec=IndexFile)
    config['run']['index_files'] = [index_file_mock] * 3
    config['run']['bam'] = 'test_bam_path'
    return config


@pytest.fixture()
def base_arr():
    return np.array([[1, 2], [1, 2]])


@pytest.fixture()
def counter_dicts_per_idx_file(base_arr):
    counter1_mock1 = get_counter_stub(base_arr,     'counter1')
    counter1_mock2 = get_counter_stub(2 * base_arr, 'counter1')
    counter2_mock1 = get_counter_stub(3 * base_arr, 'counter2')
    counter2_mock2 = get_counter_stub(4 * base_arr, 'counter2')

    counter_dicts_per_idx_file = [
        OrderedDict(counter1=counter1_mock1, counter2=counter2_mock1),
        OrderedDict(counter1=counter1_mock2, counter2=counter2_mock2),
    ]

    return counter_dicts_per_idx_file


class TestSumCounters:
    def test_adds_up_counter_arrays(self, base_arr, config,
                                    counter_dicts_per_idx_file):
        pileup_run_stub = PileupRunStubWithVisitorsAndCounters(config)
        pileup_run_stub.sum_up_counters(counter_dicts_per_idx_file)

        true_counter1_summed_arr = 3 * base_arr
        true_counter2_summed_arr = 7 * base_arr

        res = pileup_run_stub.summed_up_counters
        test_counter_arr1 = res['counter1'].counter_array
        test_counter_arr2 = res['counter2'].counter_array

        assert (np.array_equal(test_counter_arr1, true_counter1_summed_arr))
        assert (np.array_equal(test_counter_arr2, true_counter2_summed_arr))

    def test_leaves_counter_attributes_intact(self, base_arr, config,
                                              counter_dicts_per_idx_file):
        pileup_run_stub = PileupRunStubWithVisitorsAndCounters(config)
        pileup_run_stub.sum_up_counters(counter_dicts_per_idx_file)
        res = pileup_run_stub.summed_up_counters
        assert res['counter1'].config_attribute == 'counter1'
        assert res['counter2'].config_attribute == 'counter2'

    def test_returns_empty_dict_when_no_counters_present(self, config):
        pileup_run_stub = PileupRunStubNoCounters(config)
        pileup_run_stub.sum_up_counters([{}, {}])
        assert pileup_run_stub.summed_up_counters == {}

    def test_ignores_visitors_which_are_not_counters(
            self, counter_dicts_per_idx_file, config):
        pileup_run_stub = PileupRunStubWithVisitorsAndCounters(config)
        pileup_run_stub.sum_up_counters(counter_dicts_per_idx_file)
        assert (list(pileup_run_stub.summed_up_counters.keys())
                == ['counter1', 'counter2'])

def motif_pileup_mock(start):
    m = MagicMock()
    m.idx_pos.start = start
    return m

class TestSingleRun:

    def test_for_all_motif_pileups_calls_process_method_on_all_visitors(
            self, mocker, config):

        counter_mock = MagicMock(spec_set=['process', '__class__'],
                                 __class__=Counter)
        visitor_mock = MagicMock(spec_set=['process', '__class__'],
                                 __class__=Visitor)

        # Index file is irrelevant because stepwise_pileup_generator is patched
        mocker.patch('mqc.mcall_run.IndexFile')
        mocker.patch('mqc.mcall_run.pysam.AlignmentFile')

        motif_pileup_stub1 = motif_pileup_mock(1)
        motif_pileup_stub2 = motif_pileup_mock(2)
        mocker.patch('mqc.mcall_run.stepwise_pileup_generator',
                     return_value=[motif_pileup_stub1,
                                   motif_pileup_stub2])

        class PileupRunStub(PileupRun):
            def _get_visitors(self, chrom):
                # empty array must match dimensions of arrays to be added
                return OrderedDict(
                    counter1=counter_mock,
                    visitor=visitor_mock,
                )
        pileup_run_stub = PileupRunStub(config)

        counters = pileup_run_stub._single_run_over_index_file(
            index_file_path='not used in this test')

        expected_calls = [call.process(motif_pileup_stub1),
                          call.process(motif_pileup_stub2)]
        assert counter_mock.mock_calls == expected_calls
        assert visitor_mock.mock_calls == expected_calls

    def test_returns_dict_of_counters_without_visitors(
            self, mocker, config):
        pileup_run_stub = PileupRunStubWithVisitorsAndCounters(config)

        mocker.patch('mqc.mcall_run.stepwise_pileup_generator',
                     return_value=[motif_pileup_mock(1),
                                   motif_pileup_mock(2)])

        # Index file and alignemtn file is irrelevant because stepwise_pileup_generator is patched
        mocker.patch('mqc.mcall_run.pysam.AlignmentFile')
        mocker.patch('mqc.mcall_run.IndexFile')

        counters = pileup_run_stub._single_run_over_index_file(
            index_file_path='not used due to IndexFile mock')

        assert list(counters.keys()) == ['counter1', 'counter2']
        assert isinstance(counters['counter1'], Counter)



