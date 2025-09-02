"""Functions for calculating Radial Wave Function values

"""

# Third party
import numpy as np
import scipy.special as sp
from numpy.typing import NDArray

def _b0_value(n: int, L: int, Z: int) -> float:
    """Returns analytical b0 value
    
    """
    return ((2**(1 + L))
            *(1/(n*sp.factorial(2*L + 1)))
            *np.sqrt((Z*sp.factorial(n + L))/(sp.factorial(n - L - 1))))

def R_values(r_max: float, prec: int,
             n: int, L: int, Z: int) -> tuple[NDArray, NDArray]:
    """Returns `prec` evenly spaced values R of Radial Wave Function from 0 to 
    r_max and the corresponding r (distance) values.
    
    """
    # Last k-index in bk series with non-zero value. Series is zero after that
    k_max = n - L - 1

    # Create column vector with all r values used in integration
    r = np.linspace(0.0, r_max, prec)[:, np.newaxis]

    # Row vector containing all integer k-values from 0 to k_max
    k = (np.arange(0.0, k_max + 1.0,))[np.newaxis, :]

    # Recursively calculate bk values if it has more than 1 value (b0)
    bk = np.empty_like(k)
    bk[0, 0] = _b0_value(n, L, Z)

    for j in range(k_max):  # Don't use k as iterator; ruins k vector
        bk[0, j + 1] = (2*(j + L - n + 1)/((j + 1)*(j + 2*L + 2)))*bk[0, j]

    # Calculate the sums in R(r)
    R_sums = np.sum(bk*((Z*r)/(n))**(k + L), axis=1)

    # Collapse r from column vector to vector so we get a vector for R values
    r = r[:,0]  
    
    # Return R Values and r values as two vectors (1D-arrays)
    return (Z/(n))*np.exp(-(Z*r)/(n))*R_sums, r
