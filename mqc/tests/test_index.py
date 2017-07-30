import textwrap
from collections import defaultdict, OrderedDict
import pytest
import tempfile
import gzip
import os
import subprocess
import pandas as pd

from mqc.index import parallel_index_generation, bed_lines_generator

@pytest.fixture(scope="module")
def fasta_output_dir():
    return tempfile.mkdtemp()


@pytest.fixture()
def fasta_seq():
    return ('GCNNN'
            'ACGTA'
            'AGCCA'
            'GCGNT'
            'CCGAA'
            'TACNCG')

@pytest.fixture(scope="module")
def fasta_file_content():
    return ('>chr1\n'
            'GCNNn\n'  
            'ACGTA\n'  
            'AGCcA\n'  
            'gcgNt\n'  
            'CCGAA\n'  
            'TAcNcg\n')

@pytest.fixture(scope="module")
def fasta_files(fasta_output_dir, fasta_file_content):

    fasta_file_path_chr1 = os.path.join(fasta_output_dir, 'chr1.fa.gz')
    with gzip.open(fasta_file_path_chr1, 'wt') as fobj:
        fobj.write(fasta_file_content)

    fasta_file_path_chrMt = os.path.join(fasta_output_dir, 'chrMt.fa.gz')
    with gzip.open(fasta_file_path_chrMt, 'wt') as fobj:
        fobj.write(fasta_file_content)

    return fasta_file_path_chr1, fasta_file_path_chrMt


@pytest.fixture(scope="module")
def index_df():
    columns = 'chr start end motif score strand triplet_seq seq_context'.split()
    rows = [
        ['chr1', '1',  '2',  'CN',  '.', '+', 'CNN', '.'],
        ['chr1', '6',  '7',  'CG',  '.', '+', 'CGT', 'NACGT'],
        ['chr1', '7',  '8',  'CG',  '.', '-', 'CGT', 'TACGT'],
        ['chr1', '11', '12', 'CHH', '.', '-', 'CTT', 'GGCTT'],
        ['chr1', '12', '13', 'CHH', '.', '+', 'CCA', ''],
        ['chr1', '13', '14', 'CHG',  '.', '+', 'CAG', ''],
        ['chr1', '15', '16', 'CHG',  '.', '-', 'CTG', ''],
        ['chr1', '16', '17', 'CG',  '.', '+', 'CGN', ''],
        ['chr1', '17', '18', 'CG',  '.', '-', 'CGC', ''],
        ['chr1', '20', '21', 'CHG',  '.', '+', 'CCG', ''],
        ['chr1', '21', '22', 'CG',  '.', '+', 'CGA', ''],
        ['chr1', '22', '23', 'CG',  '.', '-', 'CGG', ''],
        ['chr1', '27', '28', 'CN',  '.', '+', 'CNC', ''],
        ['chr1', '29', '30', 'CG',  '.', '+', 'CG$', ''],
        ['chr1', '30', '31', 'CG',  '.', '-', 'CGN', ''],
    ]
    return pd.DataFrame(rows, columns=columns)

def df_to_bed_file_content(df, colnames, n_lines=-1):
    if n_lines == -1:
        n_lines = len(df)
    ser_of_row_strings = (df.iloc[:n_lines, :].loc[:, colnames]
                          .apply(lambda ser: ser.str.cat(sep='\t'), axis=1))
    header_line = '\t'.join(colnames) + '\n'
    return header_line + ser_of_row_strings.str.cat(sep='\n') + '\n'


def run_make_index(output_path_template, fasta_output_dir, addit_run_args):
    args = (['--fasta_path_template', f'{fasta_output_dir}/chr{{chr}}.fa.gz',
            '--output_path_template', output_path_template]
            + addit_run_args)
    subprocess.check_call(['mqc', 'make_index'] + args)


def create_make_index_err_msg(expected_index_file_text,
                              generated_index_file_text):

    err_message = textwrap.dedent(f"""\
        EXPECTED:
        {expected_index_file_text}
        
        COMPUTED:
        {generated_index_file_text}
        
        """)

    return err_message


@pytest.mark.acceptance_test
class TestMqcMakeIndexTest:

    def test_create_index_with_motifs_cg_chh_chg(
            self, fasta_output_dir, fasta_files, index_df, tmpdir):

        output_path_template = f'{tmpdir}/chr{{chr}}_all-motifs.bed.gz'
        run_make_index(output_path_template, fasta_output_dir,
                       addit_run_args=['--cg', '--chh', '--chg'])


        generated_index_file_path = output_path_template.format(chr='1')
        with gzip.open(generated_index_file_path, 'rt') as index_fobj:
            generated_index_file_text = index_fobj.read()

        columns = 'chr start end motif score strand'.split()
        expected_index_file_text = df_to_bed_file_content(index_df, columns)

        err_message = create_make_index_err_msg(expected_index_file_text,
                                                generated_index_file_text)

        assert generated_index_file_text == expected_index_file_text, \
            err_message


    def test_create_index_with_motifs_cg_chh_chg_and_triplet_seq_and_seq_context_annotations(
            self, fasta_output_dir, fasta_files, index_df, tmpdir):

        output_path_template = f'{str(tmpdir)}/chr{{chr}}_all-motifs.bed.gz'
        add_run_args = ['--cg', '--chh', '--chg',
                        '--triplet_seq',
                        '--seq_context', '2']

        run_make_index(output_path_template, fasta_output_dir,
                       add_run_args)

        generated_index_file_path = output_path_template.format(chr='1')
        with gzip.open(generated_index_file_path, 'rt') as index_fobj:
            # I only added seq_context fields for the first lines in the solution
            generated_index_file_text = ''.join(index_fobj.readlines()[0:5])


        columns = 'chr start end motif score strand triplet_seq seq_context'.split()
        expected_index_file_text = df_to_bed_file_content(index_df, columns,
                                                          n_lines=4)

        err_message = create_make_index_err_msg(expected_index_file_text,
                                                generated_index_file_text)

        assert generated_index_file_text == expected_index_file_text, \
            err_message

def mock_fasta_to_index(*args, **kwargs):
    return kwargs

def test_make_index_determines_filepaths_and_calls_parallel_index_generation(tmpdir):

    fa1_tmp = tmpdir.join("chr1.fa.gz").ensure()
    fa2_tmp = tmpdir.join("chrMt.fa.gz").ensure()

    fa2_template = tmpdir.join("chr{chr}.fa.gz")
    index_output_dir = tmpdir.mkdir('indices')
    output_path_template = index_output_dir.join(
        "index_chr{chr}_cg-chh.bed.gz")

    other_args = dict(motifs=['cg', 'chh'],
                      annotations=['triplett_motifs'])

    mock_fasta_to_index_returns = parallel_index_generation(
        fasta_path_template=str(fa2_template),
        output_path_template=str(output_path_template),
        fasta_to_index_fn=mock_fasta_to_index,
        chr_prefix='chr',
        cores=2,
        **other_args)

    call_args1 = dict(fasta_fp=str(fa1_tmp),
                      output_fp=str(index_output_dir.join(
                          "index_chr1_cg-chh.bed.gz")),
                      chr_name='chr1',
                      **other_args)

    call_args2 = dict(fasta_fp=str(fa2_tmp),
                      output_fp=str(index_output_dir.join(
                          "index_chrMt_cg-chh.bed.gz")),
                      chr_name='chrMt',
                      **other_args)

    assert call_args1 in mock_fasta_to_index_returns
    assert call_args2 in mock_fasta_to_index_returns


@pytest.fixture()
def annotations_orddict():
    annotations = OrderedDict()
    annotations['seq_context'] = 2
    annotations['triplet_seq'] = True
    return annotations

class TestBedLinesGenerator:

    @staticmethod
    def get_expected_lines(index_df, expected_cols):
        return [list(tup[1:]) for tup in index_df[expected_cols].itertuples()]

    def test_finds_cytosines_and_classifies_motifs_correctly_even_at_boundaries_and_next_to_Ns(
            self, fasta_seq, index_df, annotations_orddict):

        for k in annotations_orddict.keys():
            annotations_orddict[k] = False

        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['cg', 'chh', 'chg'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = 'chr start end motif score strand'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines == expected_lines

    def test_optionally_annotates_motif_triplet_sequence(
            self, fasta_seq, index_df, annotations_orddict):
        annotations_orddict['seq_context'] = 0
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['cg', 'chh', 'chg'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = 'chr start end motif score strand triplet_seq'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines == expected_lines


    def test_optionally_annotates_motif_seq_context(
            self, fasta_seq, index_df, annotations_orddict):
        annotations_orddict['triplet_seq'] = False
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['cg', 'chh', 'chg'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = 'chr start end motif score strand seq_context'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)

        assert computed_lines[0:3] == expected_lines[0:3]

    def test_order_of_annotation_fields_is_taken_from_annotations_orddict_order(
            self, fasta_seq, index_df):

        # first seq_context, then triplet_seq
        annotations = OrderedDict()
        annotations['seq_context'] = 2
        annotations['triplet_seq'] = True
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['cg', 'chh', 'chg'],
                                                  annotations=annotations,
                                                  chr_name='chr1'))
        expected_cols = 'chr start end motif score strand seq_context triplet_seq'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines[0:3] == expected_lines[0:3]


        # first triplet_seq, then seq_context
        annotations = OrderedDict()
        annotations['triplet_seq'] = True
        annotations['seq_context'] = 2
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['cg', 'chh', 'chg'],
                                                  annotations=annotations,
                                                  chr_name='chr1'))
        expected_cols = 'chr start end motif score strand triplet_seq seq_context'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines[0:3] == expected_lines[0:3]
