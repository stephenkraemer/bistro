import mqc

import gzip


class IndexFile:
    def __init__(self, bed_abspath):
        self.index_fobj = gzip.open(bed_abspath, 'rt')
        next(self.index_fobj)  # skip header

    def __next__(self):
        next_index_line = next(self.index_fobj)
        return IndexPosition(next_index_line)

    def __iter__(self):
        return self


class IndexPosition:
    def __init__(self, index_line: str):
        fields = index_line.rstrip().split('\t')
        self.chrom = fields[0].replace('chr', '')
        self.start = int(fields[1])
        self.end = int(fields[2])
        # TODO: use name for motif? -> name as motif on the respective strand?
        self.name = fields[3]
        self.score = fields[4]
        self.strand = fields[5]  # '+' or '-'
        # TODO: read in adjacent positions with json string loading
        self.adjacent_conversion_control_positions = fields[6]

        self.original_motif = self.name
        self.watson_motif = (self.name
                             if self.strand == '+'
                             else mqc.utilities.reverse_complement_seq(self.name))

    def __str__(self):
        return '{}:{}-{}, motif: {} (on strand {})'.format(
            self.chrom, self.start, self.end,
            self.original_motif, self.strand)
