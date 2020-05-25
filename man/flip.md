**`flip`**: The argument `flip` gives an indication of if and the columns of estimated sparse component should be flipped.
Note that the estimated (sparse) loadings, i.e., the weights on original variables, are column-wise invariant to a sign flipping.
This is because flipping of a principal direction does not influence the amount of the explained variance by the component. 
If `flip=TRUE`, then the columns of loadings will be flip accordingly, such that each column is positive-skewed.
This means that for each column, the sum of cubic elements (i.e., `sum(x^3)`) are non-negative. 
