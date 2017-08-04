import shutil
import textwrap
from collections import defaultdict, OrderedDict
import pytest
import tempfile
import gzip
import os
import os.path as op
import subprocess
import pandas as pd
from unittest.mock import call

from mqc.index import start_parallel_index_generation, bed_lines_generator, find_chroms, read_fasta

@pytest.fixture()
def fasta_seq():
    return ('GCNNN'
            'ACGTA'
            'AGCCA'
            'GCGNT'
            'CCGAA'
            'TACNCG')

@pytest.fixture(scope="module")
def genome_fasta_path():

    chr_text_template = ('>chr{chr}\n'
                         'GCNNn\n'
                         'ACGTA\n'
                         'AGCcA\n'
                         'gcgNt\n'
                         'CCGAA\n'
                         'TAcNcg\n\n')

    tmpdir = tempfile.mkdtemp()
    fa_path = op.join(tmpdir, 'mm10.fa.gz')

    with gzip.open(fa_path, 'wt') as fobj:
        fobj.write(chr_text_template.format(chr='1'))
        fobj.write(chr_text_template.format(chr='Mt'))

    yield fa_path
    # TODO: is this safe against exceptions?
    shutil.rmtree(tmpdir)


@pytest.fixture()
def index_df():
    columns = '#chrom start end motif score strand triplet_seq seq_context'.split()
    rows = [
        # ['chr1', '1',  '2',  'CN',  '.', '+', 'CNN', '.'],
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
        # ['chr1', '27', '28', 'CN',  '.', '+', 'CNC', ''],
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

#TODO: is this necessary? -> just use tmpdir?
@pytest.fixture(scope='function')
def index_output_dir():
    tmpdir = tempfile.mkdtemp()
    yield tmpdir
    # TODO: is this safe against exceptions?
    shutil.rmtree(tmpdir)

def run_make_index(index_output_dir, genome_fasta, addit_run_args):
    args = (['--genome_fasta', genome_fasta,
            '--output_dir', index_output_dir]
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
            self, index_output_dir, genome_fasta_path, index_df):

        # TODO: change flags --cg etc. to uppercase
        run_make_index(index_output_dir=index_output_dir,
                       genome_fasta=genome_fasta_path,
                       addit_run_args=['--cg', '--chh', '--chg'])

        genome_name = 'mm10'
        motifs_str = 'CG-CHG-CHH'
        output_file_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
        chromMt_index_fp = output_file_template.format(chrom='chrMt')
        with gzip.open(chromMt_index_fp, 'rt') as index_fobj:
            chromMt_index_text = index_fobj.read()

        columns = '#chrom start end motif score strand'.split()
        index_df['#chrom'] = 'chrMt'
        expected_index_file_text = df_to_bed_file_content(index_df, columns)

        err_message = create_make_index_err_msg(expected_index_file_text,
                                                chromMt_index_text)

        assert chromMt_index_text == expected_index_file_text, \
            err_message


    def test_create_index_with_motifs_cg_chh_chg_and_triplet_seq_and_seq_context_annotations(
            self, index_output_dir, genome_fasta_path, index_df):

        add_run_args = ['--cg', '--chh', '--chg',
                        '--triplet_seq',
                        '--seq_context', '2']

        run_make_index(index_output_dir=index_output_dir,
                       genome_fasta=genome_fasta_path,
                       addit_run_args=add_run_args)

        genome_name = 'mm10'
        motifs_str = 'CG-CHG-CHH'
        output_file_template = f"{index_output_dir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
        generated_index_file_path = output_file_template.format(chrom='chr1')
        with gzip.open(generated_index_file_path, 'rt') as index_fobj:
            # I only added seq_context fields for the first lines in the solution
            generated_index_file_text = ''.join(index_fobj.readlines()[0:4])


        columns = '#chrom start end motif score strand triplet_seq seq_context'.split()
        expected_index_file_text = df_to_bed_file_content(index_df, columns,
                                                          n_lines=3)

        err_message = create_make_index_err_msg(expected_index_file_text,
                                                generated_index_file_text)

        assert generated_index_file_text == expected_index_file_text, \
            err_message

    def test_create_cg_only_index_without_annotations(
            self, index_df, tmpdir, genome_fasta_path):

        add_run_args = ['--cg']
        run_make_index(index_output_dir=tmpdir,
                       genome_fasta=genome_fasta_path,
                       addit_run_args=add_run_args)

        genome_name = 'mm10'
        motifs_str = 'CG'
        output_file_template = f"{tmpdir}/{genome_name}_{motifs_str}_{{chrom}}.bed.gz"
        generated_index_file_path = output_file_template.format(chrom='chr1')
        with gzip.open(generated_index_file_path, 'rt') as index_fobj:
            chr1_idx_text = index_fobj.read()

        columns = '#chrom start end motif score strand'.split()
        is_cg = index_df['motif'] == 'CG'
        cg_index_df = index_df.loc[is_cg, :]
        expected_index_file_text = df_to_bed_file_content(cg_index_df, columns)

        err_message = create_make_index_err_msg(expected_index_file_text,
                                                chr1_idx_text)

        assert chr1_idx_text == expected_index_file_text, \
            err_message

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

    def test_discards_cytosines_which_are_not_in_the_specified_motifs(
            self, index_df, fasta_seq):

        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG'],
                                                  annotations={},
                                                  chr_name='chr1'))

        expected_cols = '#chrom start end motif score strand'.split()
        is_cg = index_df['motif'] == 'CG'
        cg_index_df = index_df.loc[is_cg, :]
        expected_lines = self.get_expected_lines(cg_index_df, expected_cols)

        assert computed_lines == expected_lines

    def test_finds_cytosines_and_classifies_motifs_correctly_even_at_boundaries_and_next_to_Ns(
            self, fasta_seq, index_df, annotations_orddict):

        for k in annotations_orddict.keys():
            annotations_orddict[k] = False

        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG', 'CHH', 'CHG'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = '#chrom start end motif score strand'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines == expected_lines

    def test_optionally_annotates_motif_triplet_sequence(
            self, fasta_seq, index_df, annotations_orddict):
        annotations_orddict['seq_context'] = 0
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG', 'CHH', 'CHG'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = '#chrom start end motif score strand triplet_seq'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines == expected_lines


    def test_optionally_annotates_motif_seq_context(
            self, fasta_seq, index_df, annotations_orddict):
        annotations_orddict['triplet_seq'] = False
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG', 'CHH', 'CHG'],
                                                  annotations=annotations_orddict,
                                                  chr_name='chr1'))
        expected_cols = '#chrom start end motif score strand seq_context'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)

        assert computed_lines[0:3] == expected_lines[0:3]

    def test_order_of_annotation_fields_is_taken_from_annotations_orddict_order(
            self, fasta_seq, index_df):

        # first seq_context, then triplet_seq
        annotations = OrderedDict()
        annotations['seq_context'] = 2
        annotations['triplet_seq'] = True
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG', 'CHH', 'CHG'],
                                                  annotations=annotations,
                                                  chr_name='chr1'))
        expected_cols = '#chrom start end motif score strand seq_context triplet_seq'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines[0:3] == expected_lines[0:3]


        # first triplet_seq, then seq_context
        annotations = OrderedDict()
        annotations['triplet_seq'] = True
        annotations['seq_context'] = 2
        computed_lines = list(bed_lines_generator(fasta_seq=fasta_seq,
                                                  motifs=['CG', 'CHH', 'CHG'],
                                                  annotations=annotations,
                                                  chr_name='chr1'))
        expected_cols = '#chrom start end motif score strand triplet_seq seq_context'.split()
        expected_lines = self.get_expected_lines(index_df, expected_cols)
        assert computed_lines[0:3] == expected_lines[0:3]

def test_find_chroms(genome_fasta_path):
    chroms = find_chroms(genome_fasta_path)
    assert chroms == ['chr1', 'chrMt']

def test_read_fasta(genome_fasta_path, fasta_seq):
    chr1_seq = read_fasta(genome_fasta_path, 'chr1')
    chrMt_seq = read_fasta(genome_fasta_path, 'chrMt')
    assert chr1_seq  == fasta_seq
    assert chrMt_seq == fasta_seq

