**`normalize`**: The argument `normalize` gives an indication of if and how any normalization should be done before rotation, and then undone after rotation. 
If normalize is `FALSE` (the default) no normalization is done. 
If normalize is `TRUE` then Kaiser normalization is done. 
(So squared row entries of normalized A sum to 1.0. 
This is sometimes called Horst normalization.) 
For `rotate="absmin"`, if `normalize` is a vector of length equal to the number of indicators (i.e., the number of rows of `A`), then the columns are divided by `normalize` before rotation and multiplied by `normalize` after rotation. 
Also, If `normalize` is a function then it should take `A` as an argument and return a vector which is used like the vector above.
