"""Functions calculating energies of Harmonic Oscillator

"""

# Third party
import numpy as np
from numpy.typing import NDArray

def _T_vector(n_row: NDArray, L: int) -> NDArray:
    """Vector containing kinetic energies corresponding to each n value in n_row
    
    """
    # Using index 0 on n_row to make output 1D array instead of 2D row vector
    return (0.5)*(n_row[0]*np.pi/L)**2.0

def _V_matrix(x_col: NDArray, n_row: NDArray, dx: float, L: float) -> NDArray:
    """Matrix containing potential energies V_mn
    
    """
    S = (x_col - 0.5*L) * np.sin((np.pi/L)*x_col*n_row)
    return (1.0)/(L) * (S.T @ S) * dx

def numerical_energies(K: int, N: int, L: float, PREC: int) -> NDArray:
    """Returns K first numerically approximated energies of H.O. 
    
    PREC is how many points to use in numerical integration. Higher is better.
    """
    # Create column vector with all x values used in integration
    # We use midpoint Riemann summation, hence the + 0.5
    dx = L/PREC
    x_col = ((np.arange(PREC) + 0.5)*dx)[:, np.newaxis]

    # Create row-vector with all n values
    n_row = np.arange(1.0, N + 1.0)[np.newaxis, :]
    
    # First add potential to total energy matrix
    H = _V_matrix(x_col, n_row, dx, L)

    # Then add kinetic along the diagonal of the total energy matrix
    H[np.diag_indices(N)] += _T_vector(n_row, L)

    # Find eigenvalues of total energy matrix (eigenenergies)
    approx_energies = np.linalg.eigvalsh(H)

    # Return K first energies
    return approx_energies[:K]

def analytical_energies(K: int) -> NDArray:
    """Returns K first analytically calculated energies of H.O. 
    
    """
    return np.arange(0.0, K) + 0.5
