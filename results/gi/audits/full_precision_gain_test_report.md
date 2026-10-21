# Full-precision (unrounded) profile -- gain sanity check

12 windows (subset of the 60-window gain-stress set), gains (0.25, 0.5, 1.3, 2.0, 4.0, 7.1), no rounding at any step.

| feature | gain | median |diff| | max |diff| | max relative |diff| |
|---|---:|---:|---:|---:|
| fp_total | 0.25 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_total | 0.50 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_total | 1.30 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_total | 2.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_total | 4.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_total | 7.10 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 0.25 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 0.50 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 1.30 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 2.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 4.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_average | 7.10 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 0.25 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 0.50 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 1.30 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 2.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 4.00 | 0.000e+00 | 0.000e+00 | 0.000e+00 |
| fp_sd | 7.10 | 0.000e+00 | 0.000e+00 | 0.000e+00 |

Verdict: PASS -- all full-precision summaries agreed with the 1x baseline within atol=1e-08, rtol=1e-06 at every tested gain.

This is a mathematical-inference check on a small subset, not a scientific claim on its own -- read alongside Section 4/5 of the provenance audit and the gain-resolution sensitivity report.
