import mqc

# TODO: should I exchange 'None' with 'ref' to indicate a reference match? I am already using NA to indicate non-available values due to 'N's in the read sequence
meth_dict = {'W-BC': {'C': {'C': 'methylated',
                            'T': 'unmethylated',
                            'G': 'SNP',
                            'A': 'SNP'},
                      'G': {'C': 'SNP',
                            'T': 'SNP',
                            'G': 'None',
                            'A': 'SNP'}},
             'W-BC-Rv': {'C': {'C': 'methylated',
                               'T': 'unmethylated',
                               'G': 'SNP',
                               'A': 'SNP'},
                         'G': {'C': 'SNP',
                               'T': 'SNP',
                               'G': 'None',
                               'A': 'SNP'}},
             'C-BC': {'G': {'C': 'SNP',
                            'T': 'SNP',
                            'G': 'methylated',
                            'A': 'unmethylated'},
                      'C': {'C': 'None',
                            'T': 'SNP',
                            'G': 'SNP',
                            'A': 'SNP'}},
             'C-BC-Rv': {'G': {'C': 'SNP',
                               'T': 'SNP',
                               'G': 'methylated',
                               'A': 'unmethylated'},
                         'C': {'C': 'None',
                               'T': 'SNP',
                               'G': 'SNP',
                               'A': 'SNP'}}
             }


def call_meth_at_base(bs_seq_strand, watson_ref_base, observed_watson_base):
    meth_status_str = meth_dict[bs_seq_strand][watson_ref_base][observed_watson_base]
    return meth_status_str


def call_meth_at_base_except_when_N_was_observed(
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
