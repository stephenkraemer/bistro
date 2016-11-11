import mqc


def call_frag_conversion_error(pileup_segment: 'mqc.PileupSegment',
                               index_position: 'mqc.IndexPosition',
                               max_number_of_unconv_control_cyts):
    # are control positions also zero-based?
    if pileup_segment.aligned_segment.bs_seq_strand in ['W-BC', 'C-BC-Rv']:
        control_pos = index_position.watson_conv_control_cyts
        watson_ref_base = 'C'
    else:
        control_pos = index_position.crick_conv_control_cyts
        watson_ref_base = 'G'

    pos_in_read = pileup_segment.zero_based_pos_in_watson_seq

    control_pos_in_read = ((x + pos_in_read) for x in control_pos
                           if 0 <= (x + pos_in_read) < pileup_segment.aligned_segment.read_length)

    call_meth = mqc.methylation_calling.call_meth_at_base_except_when_N_was_observed
    meth_status_at_control_pos = (
    call_meth(bs_seq_strand=pileup_segment.aligned_segment.bs_seq_strand,
              watson_ref_base=watson_ref_base,
              observed_watson_base=pileup_segment.aligned_segment.watson_seq[curr_pos])
    for curr_pos in control_pos_in_read)

    n_unmethylated_chh_positions = sum(1 for mcall in meth_status_at_control_pos
                                       if mcall == 'methylated')
    has_conversion_error = (True
                            if n_unmethylated_chh_positions > max_number_of_unconv_control_cyts
                            else False)

    return has_conversion_error


def test_call_no_conv_errors(pileup_segment: 'mqc.PileupSegment',
                             index_position: 'mqc.IndexPosition',
                             max_number_of_unconv_control_cyts):
    return False
