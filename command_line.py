import mqc

# Trimming site arrays will have to be defined per BS-Seq strand
# is the concept of BS-Seq strands meaningful in the undirectional protocol?
trimming_site_array = [[0] * 1000, [110] * 1000]


def call():
    index_file = mqc.IndexFile('/home/stephen/projects/mqc/mqc/test/data/'
                               'chr11_16815793-16824254.cg.bed.gz')
    # TODO: engine starts with 2nd position!
    pileups = mqc.PysamPileupEngine(index_file,
                                    bam_abspath='/home/stephen/projects/mqc/mqc/test/data/'
                                                'b_cells_rep1_chr11_16815793-16824254.bam',
                                    mapq_threshold=0,
                                    phred_score_threshold=0,
                                    run_mode='call',
                                    overlap_handling_mode='discard_disagreeing_overlaps',
                                    trimming_site_array=trimming_site_array,
                                    # frag_conv_err_detect=mqc.conv_err.call_frag_conversion_error,
                                    frag_conv_err_detect=mqc.conv_err.test_call_no_conv_errors,
                                    max_number_of_unconv_control_cyts=2
                                    )
    for pileup in pileups:
        # hand over to pileup engine
        mqc.overlap_handling.add_conservative_tags_for_calling(pileup)
        print(pileup.index_position.chrom,
              pileup.index_position.start)
        # print(pileup)
        beta, n_meth, n_unmeth = pileup.get_total_meth_stats(trimming_mode='adjusted')
        print(beta, n_meth, n_unmeth)

        # try:
        #     c += 1
        # except UnboundLocalError:
        #     c = 1
        # if c > 9:
        #     break


def profile(foo_name):
    import cProfile
    import subprocess

    cProfile.runctx(foo_name, globals(), locals(), "Profile.prof")
    subprocess.call(['snakeviz', "Profile.prof"])


call()

# s = pstats.Stats("Profile.prof")
# s.strip_dirs().sort_stats("time").print_stats()
