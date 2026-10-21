# Gain sensitivity versus numerical resolution

Same 60-window gain-stress subset as Stage B and the AvgSampEn gain audit; same six non-unit gains (0.25, 0.5, 1.3, 2, 4, 7.1); Total, Average and SD SampEn only. Frozen and GI profiling evaluated at the same decimals settings used in Stage E's resolution ablation.

## Aggregate by resolution (worst case across gain x summary)

| resolution | mode | approx bins (doc) | approx bins (measured) | worst median rel.diff % | worst max rel.diff % | pass rate |
|---|---|---:|---:|---:|---:|---:|
| frozen_d3 | frozen | 80 | 84 | 630.6386 | 673.1571 | 0.000 |
| frozen_d4 | frozen | 800 | 793 | 604.1487 | 629.0289 | 0.000 |
| frozen_d5 | frozen | 7000 | 6902 | 604.2731 | 610.6597 | 0.000 |
| gi_d1 | GI | 60 | 60 | 0.0000 | 0.0000 | 1.000 |
| gi_d2 | GI | 600 | 588 | 0.0000 | 0.0000 | 1.000 |
| gi_d3 | GI | 5000 | 5258 | 0.0000 | 0.0000 | 1.000 |

## Per feature x resolution x gain

| resolution | mode | feature | gain | median rel.diff % | max rel.diff % | pass/total |
|---|---|---|---:|---:|---:|---:|
| frozen_d3 | frozen | average | 0.25 | -5.5610 | 16.9107 | 0/60 |
| frozen_d3 | frozen | average | 0.50 | -1.7075 | 5.2420 | 0/60 |
| frozen_d3 | frozen | average | 1.30 | 0.8756 | 2.7076 | 0/60 |
| frozen_d3 | frozen | average | 2.00 | 2.1196 | 9.4757 | 0/60 |
| frozen_d3 | frozen | average | 4.00 | 4.6204 | 22.0376 | 0/60 |
| frozen_d3 | frozen | average | 7.10 | 6.9842 | 31.6336 | 0/60 |
| frozen_d3 | frozen | sd | 0.25 | -6.0501 | 27.1822 | 0/60 |
| frozen_d3 | frozen | sd | 0.50 | -2.1645 | 7.7848 | 0/60 |
| frozen_d3 | frozen | sd | 1.30 | 1.4750 | 3.6651 | 0/60 |
| frozen_d3 | frozen | sd | 2.00 | 3.4401 | 9.5932 | 0/60 |
| frozen_d3 | frozen | sd | 4.00 | 7.4939 | 16.9904 | 0/60 |
| frozen_d3 | frozen | sd | 7.10 | 11.7459 | 47.7534 | 0/60 |
| frozen_d3 | frozen | total | 0.25 | -75.4433 | 77.7916 | 0/60 |
| frozen_d3 | frozen | total | 0.50 | -50.3855 | 51.0722 | 0/60 |
| frozen_d3 | frozen | total | 1.30 | 30.6864 | 31.9502 | 0/60 |
| frozen_d3 | frozen | total | 2.00 | 102.2694 | 106.8536 | 0/60 |
| frozen_d3 | frozen | total | 4.00 | 309.1764 | 320.2494 | 0/60 |
| frozen_d3 | frozen | total | 7.10 | 630.6386 | 673.1571 | 0/60 |
| frozen_d4 | frozen | average | 0.25 | -4.6197 | 18.1219 | 0/60 |
| frozen_d4 | frozen | average | 0.50 | -2.4028 | 9.1257 | 0/60 |
| frozen_d4 | frozen | average | 1.30 | 1.3317 | 7.0758 | 0/60 |
| frozen_d4 | frozen | average | 2.00 | 2.8249 | 8.4240 | 0/60 |
| frozen_d4 | frozen | average | 4.00 | 6.4896 | 15.7631 | 0/60 |
| frozen_d4 | frozen | average | 7.10 | 11.2175 | 22.9192 | 0/60 |
| frozen_d4 | frozen | sd | 0.25 | -7.5218 | 24.1345 | 0/60 |
| frozen_d4 | frozen | sd | 0.50 | -3.2441 | 18.3831 | 0/60 |
| frozen_d4 | frozen | sd | 1.30 | 1.3279 | 22.8229 | 0/60 |
| frozen_d4 | frozen | sd | 2.00 | 0.7187 | 22.1815 | 0/60 |
| frozen_d4 | frozen | sd | 4.00 | 1.7315 | 16.8880 | 0/60 |
| frozen_d4 | frozen | sd | 7.10 | 3.2287 | 16.9811 | 0/60 |
| frozen_d4 | frozen | total | 0.25 | -75.4324 | 76.3683 | 0/60 |
| frozen_d4 | frozen | total | 0.50 | -50.2428 | 51.8143 | 0/60 |
| frozen_d4 | frozen | total | 1.30 | 30.2489 | 33.4137 | 0/60 |
| frozen_d4 | frozen | total | 2.00 | 99.9024 | 106.9798 | 0/60 |
| frozen_d4 | frozen | total | 4.00 | 297.0106 | 311.3327 | 0/60 |
| frozen_d4 | frozen | total | 7.10 | 604.1487 | 629.0289 | 0/60 |
| frozen_d5 | frozen | average | 0.25 | -8.7358 | 15.2843 | 0/60 |
| frozen_d5 | frozen | average | 0.50 | -5.0569 | 8.2840 | 0/60 |
| frozen_d5 | frozen | average | 1.30 | 2.0506 | 3.9622 | 0/60 |
| frozen_d5 | frozen | average | 2.00 | 6.0041 | 12.9821 | 0/60 |
| frozen_d5 | frozen | average | 4.00 | 13.2781 | 36.7820 | 0/60 |
| frozen_d5 | frozen | average | 7.10 | 20.8841 | 65.3722 | 0/60 |
| frozen_d5 | frozen | sd | 0.25 | -4.0482 | 11.7901 | 0/60 |
| frozen_d5 | frozen | sd | 0.50 | -2.2535 | 5.9151 | 0/60 |
| frozen_d5 | frozen | sd | 1.30 | 0.7264 | 3.5285 | 0/60 |
| frozen_d5 | frozen | sd | 2.00 | 2.1896 | 6.0299 | 0/60 |
| frozen_d5 | frozen | sd | 4.00 | 4.8866 | 18.3610 | 0/60 |
| frozen_d5 | frozen | sd | 7.10 | 6.7811 | 30.1213 | 0/60 |
| frozen_d5 | frozen | total | 0.25 | -74.9803 | 75.5339 | 0/60 |
| frozen_d5 | frozen | total | 0.50 | -50.0048 | 50.5328 | 0/60 |
| frozen_d5 | frozen | total | 1.30 | 29.8921 | 30.6679 | 0/60 |
| frozen_d5 | frozen | total | 2.00 | 99.8920 | 101.0276 | 0/60 |
| frozen_d5 | frozen | total | 4.00 | 298.8050 | 301.2386 | 0/60 |
| frozen_d5 | frozen | total | 7.10 | 604.2731 | 610.6597 | 0/60 |
| gi_d1 | GI | average | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | average | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | average | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | average | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | average | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | average | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | sd | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d1 | GI | total | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | average | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | sd | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d2 | GI | total | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | average | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | sd | 7.10 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 0.25 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 0.50 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 1.30 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 2.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 4.00 | 0.0000 | 0.0000 | 60/60 |
| gi_d3 | GI | total | 7.10 | 0.0000 | 0.0000 | 60/60 |

## Questions (answer by reading the tables above, not asserted here)

1. Does finer absolute resolution (frozen_d3 -> d4 -> d5) reduce frozen gain sensitivity, and does it reach zero at any tested resolution?
2. Is GI profiling's worst-case relative deviation at or below the numerical-invariance threshold (rtol 1e-8) at every tested resolution (gi_d1, gi_d2, gi_d3), independent of bin count?
3. Is the frozen-vs-GI gap consistent with gain dependence arising from the interaction of signal scale and FIXED ABSOLUTE quantization, as opposed to quantization in general?
4. Is the frozen trend monotonic in resolution? Report what the numbers show; do not force a monotonic reading if they are not monotonic.
