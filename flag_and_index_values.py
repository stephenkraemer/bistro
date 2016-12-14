from collections import namedtuple

bsseq_strand_flags = namedtuple('BsseqStrandFlagValues',
                                ['c_bc', 'w_bc', 'w_bc_rv', 'c_bc_rv', 'mate_and_direction_bits'])(
    w_bc=96,
    c_bc=80,
    w_bc_rv=144,
    c_bc_rv=160,
    mate_and_direction_bits=240)

bsseq_strand_indices = namedtuple('BsseqIndex', 'w_bc c_bc w_bc_rv c_bc_rv')(
    c_bc=0,
    c_bc_rv=1,
    w_bc=2,
    w_bc_rv=3
)

methylation_status_flags = namedtuple('MethylationStatusFlags',
                                      ['is_na', 'is_methylated', 'is_unmethylated', 'is_snp', 'is_ref'])(
    is_na=16,
    is_methylated=8,
    is_unmethylated=4,
    is_snp=2,
    is_ref=1
)
