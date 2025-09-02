# **Exercise 2**: Plotting of Radial Wave Function

## Main program
The program `plot.py` uses numerical integration to determine the normalization of the Radial Wave Function, then plots it together with the Reduced Radial Wave Function and the Radial Probability Density.

Example of how to run the program with $n=2$, $l=0$  and $Z=1$, where we use `rmax=10` to specify the maximum length to calculate for, and `prec=100` amount of points between `0` and `rmax` for the numerical integration and plotting.

```bash
python plot.py -n 2 -L 1 -Z 1 --rmax 16 --prec 100
```

Output:

![](fig.png)


## Mathematics behind code

TODO