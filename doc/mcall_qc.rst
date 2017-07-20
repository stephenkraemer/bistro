Methylation calling with mqc mcall
##################################

Configuring the methylation caller
**********************************
Methylation calling is typically run over a cohort of samples. Thereby

* some config variables differ from sample to sample, e.g. the output directory or the sample name and other sample metadata
* some config variables typically differ between projects, but are the same for all samples, e.g. the applied sequencing type (X10, ...) and the resulting maximal read length etc.
* some config variables only need to be changed in special cases, e.g. when the standard plotting axis limits for M-bias plots do not work because the samples have unusally strong M-bias

In general, all parameters of the program, including plotting options and file paths, are completely configurable. No relevant parameter is hardcoded.
However, most users will only need to tinker with a few parameters and  sensible defaults are provided wherever possible. We try to make the relevant
parameters easily accessible using the following configuration layers:

* All parameters which often need to be changed between samples of the same project are given as run time parameters passed through the CLI, this includes

  * sample metadata, e.g. sample name
  * the output directory for mcalls and qc files and plots for the sample

* In addition, the user may supply a config file (in the simple TOML markup language)
* A typical config file with the most useful parameters is found here. This will enough for most users.
* If you need access to more internal parameters, have a look at the default config file, which lists all parameters of the methylation calling and QC algorithms. Any of there parameters can be overwritten in the user config file.

To find all CLI parameters, checkout the help for the mcall command:

Advanced configuration options
==============================

The default config file also controls the paths for all files produced during methylation calling and qc

The defaults for the output paths will usually be sufficient, but can be changed as any other config variable.
All paths are collected in the 'paths' section of the config file. In this section (and only in this section), in addition to the standard TOML syntax, variable expansion in braced fields ({}) similar to the
python format specification mini-language is available. This means that you can define paths as follows, for example:

| [paths]
| abs_path = /an/absolute/path/to/a/file
| rel_path = subdir/
| nested_file = {rel_path}/{name}/file.txt

Note that you can use

- sample metadata (like the sample name), which you have passed through the CLI
- references to other paths
- absolute paths or relative paths. Relative paths are interpreted as relative to the output directory passed to the mcall command.
