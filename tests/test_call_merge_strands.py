import gzip
import subprocess
from pathlib import Path
import pytest

from tests.test_mcall_cmd import EXPECTED_RESULTS_MERGED, EXPECTED_RESULTS_DICT2
from tests.utils import to_tsv
from mqc.merge_strands import BedMethCallsMerger


class TestCreateMergedCallsGenerator:

    def test_merges_bed(self):
        merger = BedMethCallsMerger(strand_resolved_meth_calls='_',
                                    merged_meth_calls='_',
                                    verbose=True)
        field_dict_gen = iter([
            {'#chrom': '1', 'start': '10', 'end': '11',
             'motif': 'CG', 'score': '.', 'strand': '+',
             'beta_value': 0.5, 'n_meth': 2, 'n_total': 4},
            {'#chrom': '1', 'start': '11', 'end': '12',
             'motif': 'CG', 'score': '.', 'strand': '+',
             'beta_value': 1, 'n_meth': 1, 'n_total': 1}
        ])
        computed_line = next(merger._create_merged_calls_gen(
                field_dict_gen, n_index_cols=6))
        res_field_dict = {'#chrom': '1', 'start': '10', 'end': '12',
         'motif': 'CG', 'score': '.', 'strand': '+',
         'beta_value': f'{3/5:.8f}', 'n_meth': '3', 'n_total': '5'}
        expected_line = '\t'.join(res_field_dict.values()) + '\n'
        assert computed_line == expected_line

    def test_merged_stratified_bed(self):
        field_dict_gen = iter([
            {'#chrom': '1', 'start': '10', 'end': '11',
             'motif': 'CG', 'score': '.', 'strand': '+',
             'beta_value': 0.5, 'n_meth': 2, 'n_total': 4,
             'c_bc_beta_value': 0.5,
             'c_bc_n_meth': 1,
             'c_bc_n_total': 2,
             'c_bc_rv_beta_value': float('nan'),
             'c_bc_rv_n_meth': 0,
             'c_bc_rv_n_total': 0,
             'w_bc_beta_value': float('nan'),
             'w_bc_n_meth': 0,
             'w_bc_n_total': 0,
             'w_bc_rv_beta_value': 0.5,
             'w_bc_rv_n_meth': 1,
             'w_bc_rv_n_total': 2,
             'mate1_beta_value': 0.5,
             'mate1_n_meth': 1,
             'mate1_n_total': 2,
             'mate2_beta_value': 0.5,
             'mate2_n_meth': 1,
             'mate2_n_total': 2,
             },
            {'#chrom': '1', 'start': '11', 'end': '12',
             'motif': 'CG', 'score': '.', 'strand': '+',
             'beta_value': 1, 'n_meth': 1, 'n_total': 1,
             'c_bc_beta_value': float('nan'),
             'c_bc_n_meth': 0,
             'c_bc_n_total': 0,
             'c_bc_rv_beta_value': 1,
             'c_bc_rv_n_meth': 1,
             'c_bc_rv_n_total': 1,
             'w_bc_beta_value': float('nan'),
             'w_bc_n_meth': 0,
             'w_bc_n_total': 0,
             'w_bc_rv_beta_value': float('nan'),
             'w_bc_rv_n_meth': 0,
             'w_bc_rv_n_total': 0,
             'mate1_beta_value': float('nan'),
             'mate1_n_meth': 0,
             'mate1_n_total': 0,
             'mate2_beta_value': 1,
             'mate2_n_meth': 1,
             'mate2_n_total': 1,
             }
        ])

        merger = BedMethCallsMerger(strand_resolved_meth_calls='_',
                                    merged_meth_calls='_',
                                    verbose=True)

        computed_line = next(merger._create_merged_calls_gen(
                field_dict_gen, n_index_cols=6))
        res_field_dict = {'#chrom': '1', 'start': '10', 'end': '12',
                          'motif': 'CG', 'score': '.', 'strand': '+',
                          'beta_value': f'{3/5:.8f}', 'n_meth': '3', 'n_total': '5',
                          'c_bc_beta_value': f'{1/2:.8f}',
                          'c_bc_n_meth': 1,
                          'c_bc_n_total': 2,
                          'c_bc_rv_beta_value': f'{1:.8f}',
                          'c_bc_rv_n_meth': 1,
                          'c_bc_rv_n_total': 1,
                          'w_bc_beta_value': 'nan',
                          'w_bc_n_meth': 0,
                          'w_bc_n_total': 0,
                          'w_bc_rv_beta_value': f'{1/2:.8f}',
                          'w_bc_rv_n_meth': 1,
                          'w_bc_rv_n_total': 2,
                          'mate1_beta_value':f'{1/2:.8f}',
                          'mate1_n_meth': 1,
                          'mate1_n_total': 2,
                          'mate2_beta_value': f'{2/3:.8f}',
                          'mate2_n_meth': 2,
                          'mate2_n_total': 3,
                          }
        res_field_dict = {k: str(v) for k, v in res_field_dict.items()}
        expected_line = '\t'.join(res_field_dict.values()) + '\n'
        assert computed_line == expected_line


class TestBedMethCallsMerger:
    def test_merges_bed(self, tmpdir):
        # tmpdir_dir = TemporaryDirectory()
        # tmpdir = tmpdir_dir.name
        tmpdir_path = Path(tmpdir)
        input_calls = tmpdir_path.joinpath('calls.bed.gz')
        merged_calls = tmpdir_path.joinpath('merged-calls.bed.gz')
        with gzip.open(input_calls, 'wt') as fobj:
            fobj.write(EXPECTED_RESULTS_DICT2[('CG', 'bed', '1')])
        merger = BedMethCallsMerger(strand_resolved_meth_calls=input_calls,
                                    merged_meth_calls=merged_calls,
                                    verbose=True)
        merger.run()

        with gzip.open(merged_calls, 'rt') as fin:
            computed_merged_calls = fin.read()

        assert computed_merged_calls == EXPECTED_RESULTS_MERGED[('CG', 'bed', '1')]

    def test_merges_stratified_bed(self, tmpdir):
        tmpdir_path = Path(tmpdir)
        input_calls = tmpdir_path.joinpath('calls.bed.gz')
        merged_calls = tmpdir_path.joinpath('merged-calls.bed.gz')
        with gzip.open(input_calls, 'wt') as fobj:
            fobj.write(EXPECTED_RESULTS_DICT2[('CG', 'stratified_bed', '1')])

        merger = BedMethCallsMerger(strand_resolved_meth_calls=input_calls,
                                    merged_meth_calls=merged_calls,
                                    verbose=True)
        merger.run()
        with gzip.open(merged_calls, 'rt') as fin:
            computed_merged_calls = fin.read()
        assert computed_merged_calls == EXPECTED_RESULTS_MERGED[('CG', 'stratified_bed', '1')]


@pytest.mark.acceptance_test
@pytest.mark.parametrize('verbose', [True, False])
def test_merge_tool(tmpdir, verbose):
    tmpdir = Path(tmpdir)
    strand_resolved_meth_calls_path = tmpdir.joinpath('strand_resolved_calls.bed.gz')
    merged_meth_calls_path = tmpdir.joinpath('merged_calls.bed.gz')
    with gzip.open(strand_resolved_meth_calls_path, 'wt') as fobj:
        fobj.write(merge_tool_test_files[('bed', 'strand_resolved')])
    command_list = ['bistro', 'meth_calls', 'merge',
                    '--strand_resolved', strand_resolved_meth_calls_path,
                    '--merged', merged_meth_calls_path]
    if verbose:
        command_list += ['-v']
    subprocess.run(command_list, check=True)
    with gzip.open(merged_meth_calls_path, 'rt') as fobj:
        computed_file_content = fobj.read()
    assert computed_file_content == merge_tool_test_files[('bed', 'merged')]

merge_tool_test_files = {
    ('bed', 'strand_resolved'): to_tsv('''\
    #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
    7	3000224	3000225	CG	.	+	0.50000000	2	4
    7	3000225	3000226	CG	.	-	1.00000000	2	2
    7	3000239	3000240	CG	.	+	1.00000000	1	1
    7	3000240	3000241	CG	.	-	0.50000000	1	2
    7	3000262	3000263	CG	.	+	1.00000000	1	1
    7	3000263	3000264	CG	.	-	nan	0	0
    7	3000302	3000303	CG	.	+	nan	0	0
    7	3000303	3000304	CG	.	-	0.33333333	1	3
    '''),
    ('bed', 'merged'): to_tsv('''\
    #chrom	start	end	motif	score	strand	beta_value	n_meth	n_total
    7	3000224	3000226	CG	.	+	0.66666667	4	6
    7	3000239	3000241	CG	.	+	0.66666667	2	3
    7	3000262	3000264	CG	.	+	1.00000000	1	1
    7	3000302	3000304	CG	.	+	0.33333333	1	3
    '''),
}

