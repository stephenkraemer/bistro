# Mouse samples
# sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox"
sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias/sandbox"
output_rpp_dir = f"{sandbox_dir}/results_per_pid"
index_dir = f"{sandbox_dir}/genomes/GRCm38mm10_PhiX_Lambda"
reference_genome = f"{index_dir}/GRCm38mm10_PhiX_Lambda.fa"
mqc_config_file = f"{sandbox_dir}/user_config.toml"
# in_rpp_dir = "/icgc/dkfzlsdf/analysis/hs_ontogeny/results/wgbs/results_per_pid"
in_rpp_dir = "/icgc/dkfzlsdf/project/mouse_hematopoiesis/sequencing/whole_genome_bisulfite_tagmentation_sequencing/view-by-pid"
# bam_pattern_by_pid = f"{in_rpp_dir}/{{pid}}/alignment/blood_{{pid}}_merged.mdup.bam"
bam_pattern_by_pid = f"{in_rpp_dir}/{{pid}}/blood/paired/merged-alignment/blood_{{pid}}_merged.mdup.bam"

autosomes = [str(i) for i in range(1, 20)]
other_chroms = []  #['X', 'Y', 'MT', 'phix', 'L']

single_motifs = ['CG']
motifs_msv_str = 'CG'
motifs_csv_str = 'CG'

chr_prefix = ''

"""
# Human samples
sandbox_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mcall_qc/sandbox"
output_rpp_dir = f"{sandbox_dir}/results_per_pid"
# index_dir = f"{sandbox_dir}/genomes/hg38"
# reference_genome = f"{index_dir}/hg38.fa.gz"
index_dir = f"{sandbox_dir}/genomes/hg19"
reference_genome = f"{index_dir}/hs37d5.fa"
mqc_config_file = f"{sandbox_dir}/user_config.toml"
in_rpp_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias_analysis/pid_analysis/hipo35"
bam_pattern_by_pid = f"{in_rpp_dir}/{{pid}}/alignment/{{pid}}.bam"
# in_rpp_dir = "/icgc/dkfzlsdf/analysis/B080/kraemers/projects/mbias_analysis/pid_analysis/hipo35"
# bam_pattern_by_pid = f"{in_rpp_dir}/{{pid}}/alignment/{{pid}}.bam"
# blood06_H035-049K/alignment/blood06_H035-049K_merged.mdup.bam"

autosomes = [str(i) for i in range(14,20)]
other_chroms = []  #['X', 'Y', 'MT', 'phix', 'L']

single_motifs = ['CG']
motifs_msv_str = 'CG'
motifs_csv_str = 'CG'

chr_prefix = ''
"""
