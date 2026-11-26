# Benchmark results

2026-11-26 14:51; strict=False; timeout=60.0s; baseline node cap=2000000

| family | case | n | W | exp | status | check | decided_by | s | sections (s) | fallback | baseline | bl nodes | bl s | speed-up |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| random | random_n10_Wsq | 10 | 100 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n10_W1e9 | 10 | 1.00e+09 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n30_Wsq | 30 | 900 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n30_W1e9 | 30 | 1.00e+09 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n100_Wsq | 100 | 10000 | YES | YES | ok | quick | 0.009 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n100_W1e9 | 100 | 1.00e+09 | YES | YES | ok | quick | 0.007 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| random | random_n300_Wsq | 300 | 90000 | YES | YES | ok | quick | 0.045 | 0.preprocess:0.00 1.quick:0.03 | False | - | 0 | 0.000 |  |
| random | random_n300_W1e9 | 300 | 1.00e+09 | YES | YES | ok | quick | 0.081 | 0.preprocess:0.00 1.quick:0.05 | False | - | 0 | 0.000 |  |
| homometric | homometric_pair0 | 6 | 17 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_pair1 | 6 | 11 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_pair2 | 6 | 13 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_pair3 | 6 | 14 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_prod_00 | 36 | 612 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_prod_01 | 36 | 402 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_prod_10 | 36 | 402 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| homometric | homometric_prod_012 | 216 | 5862 | YES | YES | ok | quick | 0.030 | 0.preprocess:0.00 1.quick:0.02 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz3_first | 17 | 3000 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz3_alt | 17 | 3000 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz3_rand | 17 | 3000 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz5_first | 27 | 5000 | YES | YES | ok | lift | 0.108 | 0.preprocess:0.00 1.quick:0.01 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz5_alt | 27 | 5000 | YES | YES | ok | lift | 0.011 | 0.preprocess:0.00 1.quick:0.01 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz5_rand | 27 | 5000 | YES | YES | ok | lift | 0.010 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz7_first | 37 | 7000 | YES | YES | ok | lift | 0.023 | 0.preprocess:0.00 1.quick:0.02 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz7_alt | 37 | 7000 | YES | YES | ok | lift | 0.022 | 0.preprocess:0.00 1.quick:0.01 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz7_rand | 37 | 7000 | YES | YES | ok | lift | 0.027 | 0.preprocess:0.00 1.quick:0.02 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz9_first | 47 | 9000 | YES | YES | ok | lift | 0.059 | 0.preprocess:0.00 1.quick:0.04 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz9_alt | 47 | 9000 | YES | YES | ok | lift | 0.047 | 0.preprocess:0.00 1.quick:0.04 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz9_rand | 47 | 9000 | YES | YES | ok | lift | 0.038 | 0.preprocess:0.00 1.quick:0.03 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz11_first | 57 | 11000 | YES | YES | ok | lift | 0.055 | 0.preprocess:0.00 1.quick:0.04 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz11_alt | 57 | 11000 | YES | YES | ok | lift | 0.063 | 0.preprocess:0.00 1.quick:0.05 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz11_rand | 57 | 11000 | YES | YES | ok | lift | 0.066 | 0.preprocess:0.00 1.quick:0.05 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz13_first | 67 | 13000 | YES | YES | ok | lift | 0.083 | 0.preprocess:0.00 1.quick:0.07 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz13_alt | 67 | 13000 | YES | YES | ok | lift | 0.099 | 0.preprocess:0.00 1.quick:0.09 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| zhang | zhang_nz13_rand | 67 | 13000 | YES | YES | ok | lift | 0.080 | 0.preprocess:0.00 1.quick:0.07 2.positivity:0.01 3.vieta:0.00 4.factor:0.00 5.lift:0.00 | False | - | 0 | 0.000 |  |
| narrow | N20 | 20 | 42 | YES | YES | ok | factor | 0.005 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 3.vieta:0.00 4.factor:0.00 | False | - | 0 | 0.000 |  |
| narrow | N21 | 21 | 46 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| narrow | HM24B2 | 52 | 87 | YES | YES | ok | quick | 0.007 | 0.preprocess:0.00 1.quick:0.01 | False | - | 0 | 0.000 |  |
| narrow | HM12B1 | 26 | 43 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| narrow | HM18B1 | 34 | 55 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| narrow | HM24B3 | 62 | 107 | YES | YES | ok | factor | 0.079 | 0.preprocess:0.00 1.quick:0.06 2.positivity:0.00 3.vieta:0.00 4.factor:0.02 | False | - | 0 | 0.000 |  |
| narrow | HM30B2 | 60 | 99 | YES | YES | ok | quick | 0.009 | 0.preprocess:0.00 1.quick:0.01 | False | - | 0 | 0.000 |  |
| rigid | doubling_L8 | 10 | 256 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | doubling_L16 | 18 | 65536 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | doubling_L32 | 34 | 4.29e+09 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | doubling_L64 | 66 | 1.84e+19 | YES | YES | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | hyperbolic_L2 | 19 | 526002 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | hyperbolic_L4 | 29 | 2.89e+08 | YES | YES | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | hyperbolic_L8 | 49 | 8.90e+13 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| rigid | hyperbolic_L16 | 89 | 8.45e+24 | YES | YES | ok | quick | 0.005 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| symchain | symchain_HM24B2_L6 | 68 | 343 | YES | YES | ok | quick | 0.048 | 0.preprocess:0.00 1.quick:0.05 | False | - | 0 | 0.000 |  |
| symchain | symchain_HM24B2_L10 | 76 | 4183 | YES | YES | ok | quick | 0.073 | 0.preprocess:0.00 1.quick:0.07 | False | - | 0 | 0.000 |  |
| symchain | symchain_HM24B2_L40 | 136 | 4.40e+12 | YES | YES | ok | quick | 0.179 | 0.preprocess:0.00 1.quick:0.18 | False | - | 0 | 0.000 |  |
| mirprod | mirprod_N21_J8 | 57 | 82439 | YES | YES | ok | quick | 0.004 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| mirprod | mirprod_N21_J10 | 61 | 329735 | YES | YES | ok | quick | 0.005 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| skiena | skiena_jitter_n24 | 24 | 2300164 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| skiena | skiena_primes_n24 | 24 | 2.30e+07 | YES | YES | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| skiena | skiena_jitter_n60 | 60 | 5900594 | YES | YES | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| skiena | skiena_primes_n60 | 60 | 5.90e+07 | YES | YES | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| skiena | skiena_jitter_n150 | 150 | 1.49e+07 | YES | YES | ok | quick | 0.013 | 0.preprocess:0.00 1.quick:0.01 | False | - | 0 | 0.000 |  |
| skiena | skiena_primes_n150 | 150 | 1.49e+08 | YES | YES | ok | quick | 0.017 | 0.preprocess:0.00 1.quick:0.01 | False | - | 0 | 0.000 |  |
| no | no_perturb_random_n10_Wsq | 10 | 100 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_random_n10_W1e9 | 10 | 1.00e+09 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_random_n30_Wsq | 30 | 900 | NO | NO | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_random_n30_W1e9 | 30 | 1.00e+09 | NO | NO | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_N20 | 20 | 42 | NO | NO | ok | positivity | 0.003 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_N21 | 21 | 46 | NO | NO | ok | positivity | 0.003 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz3_first | 17 | 3000 | NO | NO | ok | positivity | 0.003 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz3_alt | 17 | 3000 | NO | NO | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz3_rand | 17 | 3000 | NO | NO | ok | positivity | 0.003 | 0.preprocess:0.00 1.quick:0.00 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz5_first | 27 | 5000 | NO | NO | ok | positivity | 0.009 | 0.preprocess:0.00 1.quick:0.01 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz5_alt | 27 | 5000 | NO | NO | ok | positivity | 0.014 | 0.preprocess:0.00 1.quick:0.01 2.positivity:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz5_rand | 27 | 5000 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz7_first | 37 | 7000 | NO | NO | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz9_first | 47 | 9000 | NO | NO | ok | positivity | 0.034 | 0.preprocess:0.00 1.quick:0.03 2.positivity:0.01 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz11_first | 57 | 11000 | NO | NO | ok | positivity | 0.028 | 0.preprocess:0.00 1.quick:0.02 2.positivity:0.01 | False | - | 0 | 0.000 |  |
| no | no_perturb_zhang_nz13_first | 67 | 13000 | NO | NO | ok | positivity | 0.077 | 0.preprocess:0.00 1.quick:0.07 2.positivity:0.01 | False | - | 0 | 0.000 |  |
| no | no_perturb_homometric_pair0 | 6 | 17 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_homometric_pair1 | 6 | 11 | NO | NO | ok | quick | 0.001 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_skiena_jitter_n24 | 24 | 2300164 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_skiena_primes_n24 | 24 | 2.30e+07 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_hyperbolic_L2 | 19 | 526002 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_hyperbolic_L4 | 29 | 2.89e+08 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_perturb_mirprod_N21_J8 | 57 | 82439 | NO | NO | ok | positivity | 0.134 | 0.preprocess:0.00 1.quick:0.06 2.positivity:0.08 | False | - | 0 | 0.000 |  |
| no | no_random_n7_W30_0 | 7 | 30 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n7_W30_1 | 7 | 30 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n10_W100_0 | 10 | 100 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n10_W100_1 | 10 | 100 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n11_W121_0 | 11 | 121 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n11_W121_1 | 11 | 121 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n12_W60_0 | 12 | 60 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n12_W60_1 | 12 | 60 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n31_W961_0 | 31 | 961 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n31_W961_1 | 31 | 961 | NO | NO | ok | quick | 0.000 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n101_W10000_0 | 101 | 10000 | NO | NO | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n101_W10000_1 | 101 | 10000 | NO | NO | ok | quick | 0.002 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n101_W1000000000_0 | 101 | 1.00e+09 | NO | NO | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |
| no | no_random_n101_W1000000000_1 | 101 | 1.00e+09 | NO | NO | ok | quick | 0.003 | 0.preprocess:0.00 1.quick:0.00 | False | - | 0 | 0.000 |  |

97 cases; ok 97, open 0, WRONG 0, ERROR 0
