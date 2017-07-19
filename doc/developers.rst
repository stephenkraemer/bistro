Guide for developers
####################

PileupRuns and the basic Iterator-Visitor pattern
*************************************************
Many algorithms parsing WGBS data can be formulated as a sequence of
:class:`~mqc.PileupRun` instances.

Every :class:`~mqc.PileupRun` uses a sort of Iterator-Visitor pattern to
process the data. First, for every index position, a :class:`~mqc.MotifPileup`
is generated. Then the :class:`~mqc.MotifPileup` is handed to the :method:`~mqc.Visitor.process`
method of the :class:~`mqc.Visitor` instances associated with the :class:`~mqc.PileupRun`.

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
*****************************************
.. automodule:: mqc.visitors


Important classes and functions
*******************************

.. autoclass:: mqc.PileupRun
   :members:
   :private-members:
