# statsmodels cross-check of the EHG9 conditional fits

statsmodels 0.14.5; python 3.13.9; MixedLM, REML, woman intercept + recording variance component.
Tolerance: |coefficient difference| <= 0.005, |SE difference| <= 0.005.

**Overall: ALL WITHIN TOLERANCE**

```
     outcome  adjusted_log_RMS        term  statsmodels_estimate  custom_estimate  abs_diff_estimate  statsmodels_se  custom_se  abs_diff_se  sm_converged  within_tolerance
frozen_total             False UC_minus_IN               0.40918          0.40918            0.00000         0.05782    0.05782      0.00000          True              True
frozen_total             False FM_minus_IN               0.07766          0.07766            0.00000         0.05554    0.05549      0.00005          True              True
frozen_total              True UC_minus_IN               0.01967          0.01967            0.00000         0.05427    0.05408      0.00019          True              True
frozen_total              True FM_minus_IN              -0.12158         -0.12158            0.00000         0.05084    0.05081      0.00004          True              True
frozen_total              True log_RMS_std               0.52851          0.52850            0.00000         0.02027    0.01748      0.00279          True              True
    gi_total             False UC_minus_IN              -0.11306         -0.11305            0.00000         0.04609    0.04606      0.00004          True              True
    gi_total             False FM_minus_IN              -0.20034         -0.20034            0.00000         0.04463    0.04460      0.00003          True              True
    gi_total              True UC_minus_IN               0.27501          0.27501            0.00000         0.04023    0.04021      0.00002          True              True
    gi_total              True FM_minus_IN               0.00247          0.00247            0.00000         0.03815    0.03811      0.00004          True              True
    gi_total              True log_RMS_std              -0.56733         -0.56733            0.00000         0.01432    0.01426      0.00006          True              True
```
