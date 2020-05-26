**`rotate`**: The `rotate` option specifies the rotation technique to use.
Currently, there are two build-in options---"varimax" and "absmin".
The "varimax" rotation maximizes the element-wise L4 norm of the rotated matrix. It is faster and computationally more stable. 
The "absmin" rotation minimizes the absolute sum of the rotated matrix.
It is sharper (as it directly minimizes the L1 norm) but slower and computationally less stable.
