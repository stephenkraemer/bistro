"""
# TODO: improve type hinting and imports problem
from . import pileup
# TODO: PileupSegment only required for type hinting
from .snp_calling import get_snp_score
from .meth_pileup import MethPileup
from .pileup_engines import PysamPileupEngine

# import pyximport

# pyximport.install()
# from .segments import PileupSegment, PysamPileupSegment

import mqc.methylation_calling
import mqc.overlap_handling
import mqc.utilities
import mqc.conv_err
"""
from .index import IndexFile, IndexPosition
from .parallel_meth_calling import call_methylation_parallelized
