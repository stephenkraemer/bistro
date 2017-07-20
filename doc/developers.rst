Guide for developers
####################

PileupRuns and the basic Iterator-Visitor pattern
*************************************************
Many algorithms parsing WGBS data can be formulated as a sequence of
:class:`~mqc.PileupRun` instances.

Every :class:`~mqc.PileupRun` uses a sort of Iterator-Visitor pattern to
process the data. First, for every index position, a :class:`~mqc.MotifPileup`
is generated. Then the :class:`~mqc.MotifPileup` is handed to the :func:`~mqc.Visitor.process`
method of the :class:`~mqc.Visitor` instances associated with the :class:`~mqc.PileupRun`.

During a :class:`~mqc.PileupRun`, for every :class:`~mqc.MotifPileup` a
sequence of visitors  (base class: :class:`~mqc.Visitor`) is activated in
order by calling their :func:`~mqc.Visitor.process` method. While a Visitor
processes a MotifPileup, it may

- modify both attributes of the MotifPileup and the reads accessible through , such as the coverage and
  methylation status at the corresponding index positon
- update its own attributes, typically statistics such as the M-bias stats

Many visitors count events observed in the MotifPileups, such as the average
methylation (resulting in a beta value distribution) or the coverage. More
complex counters count methylation events stratified by read position,
fragment length, phred score etc. For these visitors, the base class
:class:`~mqc.Counter` exists.

The data collected by such counters are ideally represented in dataframes,
with one column per identifier variable and one count column, e.g.

+-----------------+------------------+--------------------+-------+
| Fragment length | Position in read | Methylation status | Count |
+=================+==================+====================+=======+
| 100             | 1                | methylated         | 7000  |
+-----------------+------------------+--------------------+-------+
| 100             | 1                | unmethylation      | 3000  |
+-----------------+------------------+--------------------+-------+
| 100             | 2                | methylated         | 6910  |
+-----------------+------------------+--------------------+-------+
| ...             | ...              | ...                | ...   |
+-----------------+------------------+--------------------+-------+

Internally, such data are better represented by a multidimensional array with
one dimension per identifier variable, and the counts as values.
"""

More details about the design of Visitors
=========================================




Handling of configuration variables
***********************************

WGBS data processing tasks are usually run over a cohort of samples.
Typically, there a configuration values which

1. change from sample to sample within the same project
2. change between runs of the software, e.g. the number of cores available for the run
3. change from project to project, but remain the same between samples of the same project
4. can often be left at defaults across projects, but sometimes have to be adapted

This is handled as follows

1. Sample and run-specific config variables are passed through the CLI
2. Project-specific config variables are defined in a TOML config file provided by the user
3. Every mqc tool defines all its parameters from file paths to plotting parameters
   through a single default config file, nothing is hardcoded. This allows the user to work with
   sensible defaults, while at the same time allowing easy access to all parameters.

The function :func:`~mqc.config.assemble_config_vars`parses and integrates
all these configuration levels and provides them in one dict for use in the
program. The dict always contains the sections 'sample' and 'run' with
sample metadata and run parameters passed through the CLI, and additionally
all sections specified in the default config file. Config variables
specifying file paths can contain braced fields pointing to other file path
variables or sample metadata. Such fields are expanded recursively. See :func:`~mqc.config.assemble_config_vars`
for more details on how config variables are defined, assembled and (in the
case of file paths) expanded.


