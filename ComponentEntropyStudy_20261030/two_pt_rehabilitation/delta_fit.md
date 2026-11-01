# Post-hoc global partition-exponent calibration

This delta fit uses TI results from the same seven independent states, so it is exploratory and cannot be treated as prespecified validation. L5 is the exact mechanical duplicate of U4 and is a null check only. The published protocol used delta=1 and an unfitted delta=1.5 sensitivity.

Unweighted global delta* = **1.1989540**. Replica-index fits: +1.198609 [+1.169061, +1.228158] (nominal 95% Student-t CI, df=2). Inverse-variance weighted fit: 1.2875862; weights are fixed using the replica-mean bias variance at the unweighted optimum. Leave-one-state-out RMS held-out bias: 0.03033 kB.

## Absolute-state bias (2PT − TI), kB per particle

| State | Matching delta | Replica-root 95% CI | delta=1 | delta=1.25 | delta=1.5 | delta* |
|---|---:|---:|---:|---:|---:|---:|
| B3 | 1.30883 | [1.28832, 1.32951] | -0.10344 [-0.10708, -0.09979] | -0.01780 [-0.02340, -0.01219] | +0.05258 [+0.04499, +0.06016] | -0.03392 [-0.03910, -0.02873] |
| B5 | 1.32090 | [1.15372, 1.49187] | -0.10549 [-0.14993, -0.06104] | -0.02104 [-0.06864, +0.02656] | +0.04831 [-0.00204, +0.09867] | -0.03693 [-0.08392, +0.01006] |
| U4 | 1.15707 | [1.07037, 1.24523] | -0.04900 [-0.07291, -0.02510] | +0.02614 [-0.00077, +0.05304] | +0.08758 [+0.05815, +0.11701] | +0.01203 [-0.01431, +0.03836] |
| L5 | 1.15853 | [1.07988, 1.23842] | -0.04943 [-0.07089, -0.02796] | +0.02571 [+0.00125, +0.05018] | +0.08716 [+0.06017, +0.11414] | +0.01160 [-0.01229, +0.03549] |
| M5 | 1.15485 | [1.01092, 1.30217] | -0.04657 [-0.08525, -0.00788] | +0.02579 [-0.01630, +0.06788] | +0.08489 [+0.03990, +0.12989] | +0.01220 [-0.02924, +0.05364] |
| H5 | 1.11976 | [1.01755, 1.22373] | -0.03404 [-0.06069, -0.00740] | +0.03335 [+0.00428, +0.06242] | +0.08833 [+0.05719, +0.11948] | +0.02071 [-0.00790, +0.04931] |
| U6 | 1.18306 | [1.06971, 1.29806] | -0.05535 [-0.08565, -0.02504] | +0.01824 [-0.01355, +0.05002] | +0.07837 [+0.04517, +0.11156] | +0.00442 [-0.02707, +0.03591] |
| H7 | 1.06513 | [0.99842, 1.13268] | -0.01938 [-0.03822, -0.00054] | +0.04960 [+0.02871, +0.07049] | +0.10591 [+0.08328, +0.12854] | +0.03665 [+0.01615, +0.05715] |

| Exponent | RMS state bias | Maximum absolute state bias |
|---|---:|---:|
| 1.0 | 0.06654 | 0.10549 |
| 1.25 | 0.02931 | 0.04960 |
| 1.5 | 0.08030 | 0.10591 |
| fitted | 0.02563 | 0.03693 |

## Paired increments, kB per particle

| Comparison | TI | 2PT delta=1.25 | Bias [95% CI] | 2PT delta* | Bias [95% CI] |
|---|---:|---:|---:|---:|---:|
| M5/U4 | +0.03039 | +0.03004 | -0.00035 [-0.05602, +0.05532] | +0.03057 | +0.00018 [-0.05464, +0.05500] |
| H5/U4 | +0.10322 | +0.11043 | +0.00721 [-0.01434, +0.02877] | +0.11190 | +0.00868 [-0.01265, +0.03002] |
| H7/U6 | +0.10331 | +0.13467 | +0.03136 [+0.01067, +0.05205] | +0.13554 | +0.03223 [+0.01149, +0.05297] |
| H5/B3 | +0.30214 | +0.35329 | +0.05115 [+0.01816, +0.08414] | +0.35676 | +0.05462 [+0.02226, +0.08699] |
| H7/B5 | +0.30149 | +0.37213 | +0.07064 [+0.04374, +0.09754] | +0.37508 | +0.07358 [+0.04690, +0.10027] |
| B5/B3 | -0.00320 | -0.00645 | -0.00324 [-0.04671, +0.04023] | -0.00622 | -0.00301 [-0.04598, +0.03995] |
| L5/U4 (null) | +0.00042 | -0.00000 | -0.00042 [-0.00304, +0.00219] | -0.00000 | -0.00042 [-0.00304, +0.00219] |

## Leave-one-state-out prediction

| Held-out state | Training delta | Held-out bias |
|---|---:|---:|
| B3 | 1.17550 | -0.04155 |
| B5 | 1.17395 | -0.04496 |
| U4 | 1.20597 | +0.01400 |
| M5 | 1.20572 | +0.01404 |
| H5 | 1.20945 | +0.02335 |
| U6 | 1.20146 | +0.00511 |
| H7 | 1.21816 | +0.04159 |

Mean state-matching delta without monovalent salt (B3, U4, M5, H5): 1.18513; with monovalent salt (B5, U6, H7): 1.18970. The state roots vary from 1.065 (H7) to 1.321 (B5), while the two salt-group means are nearly equal. This is state variation without a simple salt/no-salt separation; the group means are descriptive, not an independent significance test.

Replica bands and all species entropy/fluidicity curves are in [delta_fit.json](delta_fit.json). Figures: [bias curves](delta_fit.pdf) and [absolute entropy comparison](entropy_validation_delta_fit.pdf). Intervals quantify three-replica sampling only; calibration and cross-validation use the same finite set of states.
