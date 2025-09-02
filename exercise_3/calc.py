"""Short program using Finite Central Difference in matrix form to numerically 
solve the Radial Eigen Equation

"""

import numpy as np

if __name__ == "__main__":
    
    # Quantum numbers
    L = 1

    # Discretize N r values from dr to R. 
    R = 50
    N = 100
    dr = R/N
    r = np.arange(dr, R + dr, dr)

    # Calculate discrete effective potentials for all r values
    V_eff = (L**2 + L)/(2*(r**2)) - 1/r

    # Hamiltonian matrix. First added diag is potential, three last is kinetic
    H = np.diag(V_eff)                # V_eff along main diag.
    H += np.diag(np.ones(N-1), k=-1)  # Ones along diag. 1 below main diag.
    H += np.diag(np.full(N, -2))      # Minus two's along the main diag.
    H += np.diag(np.ones(N-1), k=1)   # Ones along diag. 1 above main diag.

    # Find eigenvalues E and eigenvectors u
    E, u = np.linalg.eigh(H)
