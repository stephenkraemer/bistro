import os
import pytest
import pytoml
import time

import mqc

test_dir = os.path.dirname(__file__)
bam_path = os.path.join(test_dir, 'test.bam')
index_file_path = os.path.join(test_dir, 'test_index.bed.gz')
sample_name = 'test_sample'

@pytest.fixture(scope="session", autouse=True)
def config():
    test_dir = os.path.dirname(__file__)
    config_file = os.path.join(test_dir, '../config.default.toml')
    with open(config_file) as f_toml:
        config_dict = pytoml.load(f_toml)
    return config_dict


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


@pytest.fixture(scope="session", autouse=True)
def test_output_dir():
    timestamp = time.strftime('%d_%H')
    outdir = os.path.join('/home/kraemers/temp', 'mqc_test_' + timestamp + '/')
    os.makedirs(outdir, exist_ok=True)
    print(f'Saving objects created during testing in {outdir}')
    return outdir


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
