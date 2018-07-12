Methylation calling and comprehensive QC with the callme tools
**************************************************************

Analysis of multidimensional M-bias stats with evaluate_mbias
=============================================================

The appropriate algorithm and parameters for detecting M-bias
cutting sites depends on your data. You can specify an algorithm
and its parameters by passing the the --plateau_detection option.
It expects a JSON dict as argument, specifying the name of the algorithm
and possible arguments to the algorithm.

See plateau_detection.rst for the available algorithms and their parameters.

Currently, only the "binomp" algorithm, which determines the most likely
true global methylation and calculates binomial P-values for deviations from
this plateau, is ready for use.

Example:
.. highlight:: bash

   callme evaluate_mbias \
     # other options
     --plateau_detection '{"algorithm": "binomp", "min_plateau_length": 30,
                          "max_slope": ' '0.0006,' '"plateau_flen": 210,
                          "plateau_bs_strands": ["w_bc", "c_bc"]}'


callme API
==========
.. click:: mqc.cli:mqc
   :prog: mqc
   :show-nested:


Advanced configuration options
==============================

The defaults for the output paths will usually be sufficient, but
can be changed. Variable expansion in braced fields ({}) similar
to the python format specification mini-language is available.
This means that you can define paths as follows, for example:

| [paths]
| abs_path = /an/absolute/path/to/a/file
| rel_path = subdir/
| nested_file = {rel_path}/{name}/file.txt

Note that you can use

- sample metadata (like the sample name), which you have passed through the CLI
- references to other paths
- absolute paths or relative paths. Relative paths are interpreted as relative to the output directory passed to the mcall command.
