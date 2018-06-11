import pytest
from mqc.cli import dict_from_kwarg_cli_option


def test_dict_from_kwarg_cli_option():
    value = 'key1=value1,key2=value2'

    computed_dict = dict_from_kwarg_cli_option(value)
    expected_dict = {'key1': 'value1',
                     'key2': 'value2'}

    assert computed_dict == expected_dict


class TestCallTool:
    @pytest.mark.xfail(strict=True, skip=True)
    def test_raises_if_unknown_output_format(self):
        raise NotImplemented

    @pytest.mark.xfail(strict=True, skip=True)
    def test_raises_if_bed_and_stratified_bed_requested_together(self):
        raise NotImplemented
