""" Various processing steps not currently covered by the bistro snakefile

- merging: from strand-resolved, per (pid, chrom) to strand-merged, per pid
- create IGV tracks

Background:
I have used the demo snakefile from the bistro tests to create these meth calls. See
labbook for details. The demo snakefile is in principle ok, but in the future, I will
create a fully usable world-facing snakefile for bistro, which will then likely also
comprise these steps here. Also, in this new snakefile, I will unify the output file
paths for BED and stratified BED, which will then require adaptations in the rules/paths used
here too.

snakemake \
  --snakefile /home/kraemers/projects/mouse_hematopoiesis/src/mouse_hematopoiesis/wgbs/meth_calling2/bistro_workflow_add-ons.smk \
  --configfile /path/to/json \
  --printshellcmds \
  --latency-wait 60 \
  --jobs 1000 \
  --jobscript /home/kraemers/projects/mouse_hematopoiesis/src/mouse_hematopoiesis/wgbs/meth_calling2/jobscript.sh \
  --cluster "bsub -R rusage[mem={params.avg_mem}] -M {params.max_mem} -n {params.cores} -J {params.name} -W {params.walltime} -o /home/kraemers/temp/logs/" \
  --forcerun merge_chroms \

  --dryrun

  --forcerun create_beta_value_igv \
  --forcerun create_coverage_igv \

"""

# Implementation notes:
# - the standard bed CG calls have a different path format than the bismark and stratified bed calls,
#   this may have to be unified in the future
# - in principle, all paths should be a function of strand_resolved_calls_by_pid_chrom, currently
#   all paths are hardcoded

rpp_dir = config['rpp_dir']
chromosomes = config['chromosomes']
pids = config['pids']

strand_resolved_calls_by_pid_chrom = rpp_dir + '/{pid}/meth/meth_calls/mcalls_{pid}_CG_{chrom}.bed.gz'
# IMPORTANT: in rule 'merge_chroms' I assume that
#  strand_merged_calls_by_pid_chrom == ${strand_resolved_calls_by_pid_chrom%.bed.gz}_strands-merged.bed.gz
strand_merged_calls_by_pid_chrom = rpp_dir +  '/{pid}/meth/meth_calls/mcalls_{pid}_CG_{chrom}_strands-merged.bed.gz'
strand_merged_calls_by_pid = rpp_dir +  '/{pid}/meth/meth_calls/mcalls_{pid}_CG_chrom-merged_strands-merged.bed.gz'

coverage_igv_by_pid = strand_merged_calls_by_pid.replace('.bed.gz', '_coverage.igv')
beta_value_igv_by_pid = strand_merged_calls_by_pid.replace('.bed.gz', '_beta_value.igv')
coverage_tdf_by_pid = strand_merged_calls_by_pid.replace('.bed.gz', '_coverage.tdf')
beta_value_tdf_by_pid = strand_merged_calls_by_pid.replace('.bed.gz', '_beta_value.tdf')

rule all:
    input:
      expand(strand_merged_calls_by_pid, pid=pids),
      expand(coverage_tdf_by_pid, pid=pids),
      expand(beta_value_tdf_by_pid, pid=pids),


rule merge_chrom_strands:
    input:
        expand(strand_resolved_calls_by_pid_chrom,
               pid='{pid}',
               chrom=chromosomes)
    output:
        expand(strand_merged_calls_by_pid_chrom,
               pid='{pid}',
               chrom=chromosomes),
    params:
        walltime='00:30',
        avg_mem=4000,
        max_mem=6000,
        cores=1,
        name='{pid}_merge-strands',
    shell:
        """
        for f in {input}; do
          echo "convert $f"
          output=${{f%.bed.gz}}_strands-merged.bed.gz
          echo "write $output"
          bistro meth_calls merge \
          --strand_resolved $f \
          --merged $output \
          -v
        done
        """

rule merge_chroms:
    input:
        expand(strand_merged_calls_by_pid_chrom,
               pid='{pid}',
               chrom=chromosomes),
    params:
        walltime='00:10',
        avg_mem=1500,
        max_mem=3000,
        cores=1,
        name='{pid}_merge-chroms',
    output:
        strand_merged_calls_by_pid,
    shell:
        """
        # disable pipefail from bash strict mode
        set +o pipefail
        {{ zcat {input[0]} | head -n 1 ;
        zcat {input} | grep -v '#'; }} | gzip > {output}
        """

rule create_beta_value_igv:
    input: strand_merged_calls_by_pid
    output:
        temp(beta_value_igv_by_pid)
    params:
        avg_mem=2000,
        max_mem=4000,
        walltime='00:15',
        name='beta_igv_{pid}',
        cores=1,
    shell:
        """
        zcat {input} | awk -v pid="{wildcards.pid}" -v OFS=$'\t' '
        BEGIN {{
        print "#track viewLimits=0:1"
        print "Chromosome", "Start", "End", "Feature",
        pid "_beta_value"
        }}
        /^[^#]/ {{print $1, $2, $3, $1 "_" $2 "_" $3, $7}}' > {output}
        """

rule create_beta_value_tdf:
    input: beta_value_igv_by_pid
    output:beta_value_tdf_by_pid
    params:
        avg_mem=2000,
        max_mem=4000,
        walltime='00:15',
        name='beta_tdf_{pid}',
        cores=1,
    shell:
        """
        /home/kraemers/programs/IGVTools/igvtools toTDF {input} {output} mm10
        """

rule create_coverage_igv:
    input: strand_merged_calls_by_pid
    output:
        temp(coverage_igv_by_pid)
    params:
        avg_mem=2000,
        max_mem=4000,
        walltime='00:15',
        name='coverage_igv_{pid}',
        cores=1,
    shell:
      """
      zcat {input} | awk -v pid="{wildcards.pid}" -v OFS=$'\t' '
      BEGIN {{
      print "Chromosome", "Start", "End", "Feature",
      pid "_coverage"
      }}
      /^[^#]/ {{print $1, $2, $3, $1 "_" $2 "_" $3, $9}}' > {output}
      """

rule create_coverage_tdf:
    input: coverage_igv_by_pid
    output: coverage_tdf_by_pid
    params:
        avg_mem=2000,
        max_mem=4000,
        walltime='00:15',
        name='coverage_tdf_{pid}',
        cores=1,
    shell:
      """
      /home/kraemers/programs/IGVTools/igvtools toTDF {input} {output} mm10
      """
