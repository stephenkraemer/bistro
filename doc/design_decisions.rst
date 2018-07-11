Design decisions (with some discussion)
=======================================

Overlap handling
----------------

Phred score adjustment
^^^^^^^^^^^^^^^^^^^^^^

Currently, no phred score adjustment is performed by default. That is
we don't integrate the phred score information from overlapping mates.
This may be changed in the future.

Dealing with unmatched bases in overlapping reads
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Currently, any overlap event where the base calls don't match lead to
exclusion of both reads. This could be replaced by taking the better
quality base call instead. May be combined with phred score adjustment
(see above)

Discarding reads from considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A read cannot contribute to an overlap, if the pileup position is
trimmed, has a qc_fail or has no methylation status (NA). Reads are not
discarded if they have a phred score below the phred score filtering
threshold used for methylation calling. The rationale is that also
lower phred scores contain information useful for judging the reliability
of the methylation call from the read pair.
