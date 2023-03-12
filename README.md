# Tutorial for measurement alignment with R (sirt)
 
In this code (Alignment_Test.R), I will demonstrate how to perform:

1. Multigroup confirmatory factor analysis: see lines below "1. MG-CFA"
2. Measurement Alignment: see lines below "2. Measurement Alignment"
3. Monte Carlo simulations: see lines below "3. Monte Carlo Simulations"
4. Calculation of adjusted factor scores: see lines below "4. Factor Score Calculation"

to address the potential non-invariance issue in cross-national and/or cross-cultural research. To test the source code, the following R packages are required: lavaan, sirt, foreach, parallel, doParallel, psych, MASS.

This tutorial uses the open dataset used in Fischer and Karl (2019).

# References
Fischer, R., & Karl, J. A. (2019). A primer to (cross-cultural) multi-group invariance testing possibilities in R. Frontiers in psychology, 10, 1507. https://doi.org/10.3389/fpsyg.2019.01507
