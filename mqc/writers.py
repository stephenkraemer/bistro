import gzip

from mqc.visitors import Visitor
from mqc.pileup.pileup import MotifPileup

class BedWriter(Visitor):

    def __init__(self, config, chrom: str):
        meth_calls_path_template = (
            config['paths']['meth_calls_basepath']
            + f"{config['meta']['name']}_{{motif}}_{chrom}.bed.gz")

        self.meth_calls_fobj_dict = {
            motif: gzip.open(meth_calls_path_template.format(motif=motif), 'wt')
            for motif in config['run']['motifs']}

        header = '\t'.join(['#chrom', 'start', 'end',
                            'motif', 'score', 'strand',
                            'beta_value', 'n_meth', 'n_unmeth'])
        for fobj in self.meth_calls_fobj_dict.values():
            fobj.write(header)

    def process(self, motif_pileup: MotifPileup):
        motif = motif_pileup.idx_pos.motif
        line_no_endl = '\t'.join([
            motif_pileup.idx_pos.chrom,
            str(motif_pileup.idx_pos.start),
            str(motif_pileup.idx_pos.end),
            motif,
            '.',
            motif_pileup.idx_pos.strand,
            str(motif_pileup.beta_value),
            str(motif_pileup.n_meth),
            str(motif_pileup.n_unmeth),
        ])
        self.meth_calls_fobj_dict[motif].write(line_no_endl + '\n')

