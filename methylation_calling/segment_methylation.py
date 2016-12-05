import mqc

meth_dict = {'W-BC': {'C': {'C': 'methylated',
                            'T': 'unmethylated',
                            'G': 'SNP',
                            'A': 'SNP',
                            'N': 'NA'},
                      'G': {'C': 'SNP',
                            'T': 'SNP',
                            'G': 'Ref',
                            'A': 'SNP',
                            'N': 'NA'}},
             'W-BC-Rv': {'C': {'C': 'methylated',
                               'T': 'unmethylated',
                               'G': 'SNP',
                               'A': 'SNP',
                               'N': 'NA'},
                         'G': {'C': 'SNP',
                               'T': 'SNP',
                               'G': 'Ref',
                               'A': 'SNP',
                               'N': 'NA'}},
             'C-BC': {'G': {'C': 'SNP',
                            'T': 'SNP',
                            'G': 'methylated',
                            'A': 'unmethylated',
                            'N': 'NA'},
                      'C': {'C': 'Ref',
                            'T': 'SNP',
                            'G': 'SNP',
                            'A': 'SNP',
                            'N': 'NA'}},
             'C-BC-Rv': {'G': {'C': 'SNP',
                               'T': 'SNP',
                               'G': 'methylated',
                               'A': 'unmethylated',
                               'N': 'NA'},
                         'C': {'C': 'Ref',
                               'T': 'SNP',
                               'G': 'SNP',
                               'A': 'SNP',
                               'N': 'NA'}}
             }


def call_meth_at_base(bs_seq_strand, watson_ref_base, observed_watson_base):
    meth_status_str = meth_dict[bs_seq_strand][watson_ref_base][observed_watson_base]
    return meth_status_str


def call_meth_at_pos(pileup_read, watson_ref_base):
    direction_str = 'reverse' if pileup_read.alignment.is_reverse else 'forward'
    mate = 1 if pileup_read.alignment.is_read1 else 2
    bs_seq_strand = bs_strand_dict_directional_protocol[mate][direction_str]
    pos_in_read = pileup_read.query_position
    if not pos_in_read:
        meth_status_str = 'NA'
    else:
        observed_watson_base = pileup_read.alignment.query_sequence[pos_in_read]
        if observed_watson_base == 'N':
            meth_status_str = 'NA'
        else:
            meth_status_str = mqc.methylation_calling.segment_methylation.meth_dict[
                bs_seq_strand][watson_ref_base][observed_watson_base]
    return meth_status_str

def call_meth_at_base_except_when_N_was_observed2(
        bs_seq_strand, watson_ref_base, observed_watson_base):
    if observed_watson_base == 'N':
        meth_status_str = 'NA'
    else:
        meth_status_str = meth_dict[bs_seq_strand][watson_ref_base][observed_watson_base]
    return meth_status_str


def add_methylation_status_to_segment(pileup_segment: mqc.PileupSegment):
    pileup_segment.meth_status_str = (meth_dict[pileup_segment.aligned_segment.bs_seq_strand]
                                      [pileup_segment.watson_ref_base]
                                      [pileup_segment.observed_watson_base])

    """
    one boolean comparison is approximately 1.5 times faster than a string comparison
    if string comparison is done, checking for methylated first, it will be faster for CG context,
    for CH context, one would have to theck unmethylated first

        if meth_status_str == 'methylated':
            pileup_segment.has_meth_info = True
            pileup_segment.is_methylated = True
        elif meth_status_str == 'unmethylated':
            pileup_segment.has_meth_info = True
            pileup_segment.is_methylated = False
        elif meth_status_str == 'SNP':
            pileup_segment.has_snp = True
        else:
            # negative defaults are set at pileup segment construction
            pass
    """
