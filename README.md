# IMPORTANT NOTES OF UPDATES. MUST READ BEFORE USING THE CODE!!!!!
(Special thanks to *sirt* developer, Alexander Robitzsch)
Based on comments from Alexander, I would like to pay attention to these points while using the alignment method implemented in *sirt*:<br>

1. Altering "align.scale" argument in "invariance.alignment" will not significantly influence the outcome. Please disregard my concern regarding the selection of "align.scale" in the discussion section.
2. R^2 statistics are likely to be biased when the default loss function is used in "invariance.alignment." Alexander suggested other statistics might be more appropriate to test whether alignment was done successfully. You may refer to other statistics, such as the proportions of invariance parameters across groups or the result from Monte Carlo simulation. Alternatively, if you need to use R^2, then you should use a different loss function by specifying **"align.pow=c(2,2)"** instead of "align.pow=c(.25,.25)" in "invariance.alignment."

# Citation
Han, H. (2023). Using Measurement Alignment in Research on Adolescence Involving Multiple Groups: A Brief Tutorial with R. *Journal of Research on Adolescence*. https://doi.org/10.1111/jora.12891

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
