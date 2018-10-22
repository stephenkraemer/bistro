# from mqc.beta_values import BetaValueData, BetaValuePlotter, \
#     StratifiedBetaValueCounter
from mqc.visitors import Counter, Visitor
# from mqc.coverage import CoverageCounter
from mqc.index import IndexFile
# from mqc.mbias import (MbiasData, MbiasCounter, AdjustedMbiasCuttingSites,
#                        FixedRelativeCuttingSites)
from mqc.mcall_run import PileupRun
# from mqc.pileup.pileup import MotifPileup

from pkg_resources import resource_filename

def get_bistro_addons_snakefile():
    return resource_filename('mqc.resources', 'bistro_workflow_add-ons.smk')
