Plateau detection algorithms
****************************

Challenges
==========

- Plateau boundaries may differ between different groups of fragments (strata).
  In particular, the fragment length and BS-Seq strand often have a big influence.
- Different strata will differ in coverage, leading to different levels of noise
- The coverage may drop towards the end of the read, leading to higher noise
  towards the end of the read.

Provided algorithms
===================

While different algorithms are included in the code base, only one is currently
read for use. 

Binomial p-value based plateau detection (binomp)
-------------------------------------------------

Often, even if M-bias is present, it mainly affects certain BS-Seq strands,
and certain fragment lengths. This algorithm first estimates the unbiased global methylation level based on the set
of fragment lengths and BS-Seq strands defined as 'unbiased' by **plateau_flen** and
**plateau_bs_strands**. Then it detects significant deviations from this plateau and translates them into
BS-Seq strand and fragment length-specific cutting sites.

**Parameters**
min_plateau_length: int
   Minimum number of bp unaffected of M-bias required to include a stratum into
   the methylation calls.
allow_slope: [True, False]
    Whether the M-bias curve is allowed to have a slope.
max_slope: float
    Maximum allowed slope of the M-bias curve, if slope is allowed.
plateau_flen: int
    Fragment lengths >= plateau_flen are used to estimate the correct global methylation level
plateau_bs_strands: one or more of [w_bc c_bc w_bc_rv c_bc_rv]
   Only these BS-Seq strands are used to estimate the correct global methylation level



Providing your own algorithms
=============================
  




   




