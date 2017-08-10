import gzip
import os
import os.path as op

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup

class BedWriter(Visitor):

    def __init__(self, config, chrom: str):

        #TODO: avoid hard coding
        meth_calls_path_template = (
            config['paths']['call']['meth_calls_basepath']
            + f"_{config['sample']['name']}_{{motif}}_{chrom}.bed.gz")

        os.makedirs(op.dirname(meth_calls_path_template), exist_ok=True, mode=0o770)

        self.meth_calls_fobj_dict = {
            motif: gzip.open(meth_calls_path_template.format(motif=motif), 'wt')
            for motif in config['run']['motifs']}

        header_no_endl = '\t'.join(['#chrom', 'start', 'end',
                                    'motif', 'score', 'strand',
                                    'beta_value', 'n_meth', 'n_unmeth'])
        for fobj in self.meth_calls_fobj_dict.values():
            fobj.write(header_no_endl + '\n')

    def process(self, motif_pileup: MotifPileup):
        motif = motif_pileup.idx_pos.motif
        line_no_endl = '\t'.join([
            motif_pileup.idx_pos.chrom,
            str(motif_pileup.idx_pos.start),
            str(motif_pileup.idx_pos.end),
            motif,
            '.',
            motif_pileup.idx_pos.strand,
            f"{motif_pileup.beta_value:.6f}",
            str(motif_pileup.n_meth),
            str(motif_pileup.n_unmeth),
        ])
        self.meth_calls_fobj_dict[motif].write(line_no_endl + '\n')

