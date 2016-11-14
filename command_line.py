import gzip

import mqc

# Trimming site arrays will have to be defined per BS-Seq strand
# is the concept of BS-Seq strands meaningful in the undirectional protocol?
trimming_site_array = [[0] * 1000, [110] * 1000]


def call():
    """
    Example:
    python3 command_line.py \
    --index_file ./test/data/chr11_16815793-16824254.cg.bed.gz \
    --bam ./test/data/b_cells_rep1_chr11_16815793-16824254.bam \
    --output_file ./test/results/methylation_calls.gz
    """
    index_file = mqc.IndexFile(args.index_file)
    # TODO: engine starts with 2nd position!
    pileups = mqc.PysamPileupEngine(index_file,
                                    bam_abspath=args.bam,
                                    mapq_threshold=0,
                                    phred_score_threshold=0,
                                    run_mode='call',
                                    overlap_handling_mode='discard_disagreeing_overlaps',
                                    trimming_site_array=trimming_site_array,
                                    # frag_conv_err_detect=mqc.conv_err.call_frag_conversion_error,
                                    frag_conv_err_detect=mqc.conv_err.test_call_no_conv_errors,
                                    max_number_of_unconv_control_cyts=2
                                    )
    with gzip.open(args.output_file, 'wt') as fout:
        c = 0
        for pileup in pileups:
            c += 1
            if c % 1000 == 0:
                print(c)
            # hand over to pileup engine
            mqc.overlap_handling.add_conservative_tags_for_calling(pileup)
            beta, n_meth, n_unmeth = pileup.get_total_meth_stats(trimming_mode='adjusted')
            print(pileup.index_position.chrom,
                  pileup.index_position.start,
                  pileup.index_position.end,
                  n_meth,
                  n_unmeth,
                  beta,
                  sep='\t',
                  end='\n',
                  file=fout)

            # try:
            #     c += 1
            # except UnboundLocalError:
            #     c = 1
            # if c > 9:
            #     break


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--index_file')
    parser.add_argument('--bam', required=True)
    parser.add_argument('--output_file', required=True)
    return parser.parse_args()

args = get_args()

def profile(foo_name_str):
    import cProfile
    import subprocess

    cProfile.runctx(foo_name_str, globals(), locals(), "Profile.prof")
    # subprocess.call(['snakeviz', "Profile.prof"])


# s = pstats.Stats("Profile.prof")
# s.strip_dirs().sort_stats("time").print_stats()


if __name__ == '__main__':
    profile('call()')
