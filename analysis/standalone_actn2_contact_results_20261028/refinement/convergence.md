# Additional standalone validation

Supported translational toy rig and fixed-actin weak sites. Independent streams, three seeds per resolution; 6-SE comparisons, plus 1 nm tolerance for mean penetration. A pass excludes large observed resolution effects at this sample size; it does not establish equivalence or production-network stability.

| Metric | Coarse | Refined | Difference | SE of difference | Passed |
|---|---:|---:|---:|---:|---|
| any link loss | 0.83724 | 0.85938 | 0.0221 | 0.0183 | True |
| all links lost | 0.046875 | 0.045573 | 0.0013 | 0.0107 | True |
| mean live links | 2.1979 | 2.1484 | 0.0495 | 0.0544 | True |
| final penetration nm | 1033.7 | 1053.9 | 20.3 | 35.3 | True |
| barbed weak survival | 0 | 0 | 0 | 0 | True |
| inner weak survival | 0.30566 | 0.309 | 0.00334 | 0.00589 | True |

Energy-balance residuals (pN nm), dt / dt/2 / dt/4: [8.393291040022177, 4.186818817998983, 2.090460867764321].

Overall checks passed: True. Wall time: 890.7 s.
