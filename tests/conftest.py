import os
import pytest
import pytoml
import time
import mqc

TEST_DIR = os.path.dirname(__file__)
bam_path = os.path.join(TEST_DIR, 'test.bam')
index_file_path = os.path.join(TEST_DIR, 'test_index.bed.gz')
sample_name = 'test_sample'

timestamp = time.strftime('%d_%H')
OUTDIR = os.path.join('/home/kraemers/temp', 'mqc_test_' + timestamp)
os.makedirs(OUTDIR, exist_ok=True)
print(f'Saving objects created during testing in {OUTDIR}')

args_dict = {'sample_name': 'hsc_rep1',
             'sample_meta': 'population=hsc,replicate=1',
             'output_dir': OUTDIR}

# config_file = os.path.join(TEST_DIR, 'custom.config.toml')

CONFIG = mqc.config.get_config_dict(args_dict,
                                    config_file_path=None)

for d in ['coverage_dir', 'beta_value_dir', 'mbias_dir']:
    os.makedirs(CONFIG['paths'][d], exist_ok=True)

@pytest.fixture(scope="session", autouse=True)
def test_output_dir():
    return OUTDIR

@pytest.fixture(scope="session", autouse=True)
def config():
    return CONFIG


@pytest.fixture(scope="session", autouse=True)
def pileup_motifs_list(config):
    minimal_cutting_sites = mqc.mbias.MinimalCuttingSites(config)
    motif_pileup_iter = mqc.motif_pileup_generator(bam_path, index_file_path)
    pileup_motifs_and_idxs = []
    for motif_pileups, curr_idx_pos in motif_pileup_iter:
        mqc.pileup.annotate_pileupreads(motif_pileups=motif_pileups, index_position=curr_idx_pos,
                                        cutting_sites=minimal_cutting_sites, config=config)
        pileup_motifs_and_idxs.append([motif_pileups, curr_idx_pos])
    return pileup_motifs_and_idxs




# Marker for incremental tests (http://doc.pytest.org/en/latest/example/simple.html#incremental-testing-test-steps)

def pytest_runtest_makereport(item, call):
    if "incremental" in item.keywords:
        if call.excinfo is not None:
            parent = item.parent
            parent._previousfailed = item


def pytest_runtest_setup(item):
    if "incremental" in item.keywords:
        previousfailed = getattr(item.parent, "_previousfailed", None)
        if previousfailed is not None:
            pytest.xfail("previous test failed (%s)" % previousfailed.name)
