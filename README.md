## Improved implementation of the *Neighbourhood Proportion Error*

The original Neighbourhood Proportion Error (NPE) manuscript:

Konstorum A, Jekel N, Vidal E, Laubenbacher R (2018) Comparative Analysis of Linear and Nonlinear Dimension Reduction Techniques on Mass Cytometry Data, 
bioRxiv 273862; doi: https://doi.org/10.1101/273862

Original GitHub repo:

https://github.com/akonstodata/cytof_dimred

Modified version written by David Novak (david.novak@ugent.be)

### Modifications

* This version is faster due to extended use of built-in R functions and vectorisation.

* Pre-computed *k*-nearest-neighbour indices in high-dimensional space can be passed as input.

* Points can be excluded from the final result if they come from an 'unlabelled' pseudo-population.

* In addition to the total variation distance, earth-mover's distance can be used to compute dissimilarities between the per-population distributions of counts of like cells.

* For the final result, the per-population distribution dissimilarities can be reduced to a scalar output value either by summing *or* by taking the average.

* There is an option to produce a PDF with plots of the distributions for high-dimensional and low-dimensional coordinates for each population.
