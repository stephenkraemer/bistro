"""Simple integration test
Due to the very small size of the test file, no adjusted cutting sites can be determined
Therefore, the results for adjusted trimming represent a corner case test for the situation
where the complete sample has to be discarded
"""

import mqc
from mqc.tests.conftest import bam_path, index_file_path, sample_name

test_output_dir = mqc.tests.conftest.test_output_dir()
config = mqc.tests.conftest.config()

mqc.qc_run.qc_run(
    bam_path=bam_path,
    index_file_path=index_file_path,
    config=config,
    meth_metrics_dir_abs=test_output_dir,
    sample_name=sample_name
)
