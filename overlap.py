from collections import defaultdict
import mqc

HAS_OVERLAP = 1
BETTER_MATE_IN_OVERLAP = 2
WORSE_MATE_IN_OVERLAP = 4


def tag_overlaps(pileups: 'List[mqc.bsseq_pileup_read.BSSeqPileupRead]'):
    read_hash = defaultdict(list)
    for pileupread in pileups:
        if not pileupread.trimm_flag:
            read_hash[pileupread.alignment.query_name].append(pileupread)
    for query_name, reads in read_hash.items():
        if len(reads) == 2:
            process_overlap(reads[0], reads[1])


def process_overlap(read1, read2):
    if read1.observed_watson_base == reverse_complement_seq(read2.observed_watson_base):
            read1.overlap_flag = HAS_OVERLAP & BETTER_MATE_IN_OVERLAP
            read2.overlap_flag = HAS_OVERLAP & WORSE_MATE_IN_OVERLAP
    else:
        if read1.baseq_at_pos > read2.baseq_at_pos:
            read1.overlap_flag = HAS_OVERLAP & BETTER_MATE_IN_OVERLAP
            read2.overlap_flag = HAS_OVERLAP & WORSE_MATE_IN_OVERLAP
        elif read2.baseq_at_pos > read1.baseq_at_pos:
            read1.overlap_flag = HAS_OVERLAP & WORSE_MATE_IN_OVERLAP
            read2.overlap_flag = HAS_OVERLAP & BETTER_MATE_IN_OVERLAP
        else: # different calls, same base qualities -> discard
            read1.overlap_flag = HAS_OVERLAP & WORSE_MATE_IN_OVERLAP
            read2.overlap_flag = HAS_OVERLAP & WORSE_MATE_IN_OVERLAP


def reverse_complement_seq(seq):
    return seq.translate(str.maketrans('CGHD', 'GCDH'))[::-1]
