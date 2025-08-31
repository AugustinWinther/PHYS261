"""Harmonic Approximator

Approximates K first energy levels of the Harmonic Oscillator (H.O.) using the 
N >= K first eigenfunctions of the Infinite Square Well (I.S.W.) of length L.
"""

# Standard library
import argparse

# Local
import calc

def calculate_and_print_energies(K: int, N: int, L: float, PREC: int) -> None:
    """
    Calculates and prints first K energies
    """
    # Calculate analytical energies
    actual_energies = calc.analytical_energies(K)

    # Calculate numerical energies
    approx_energies = calc.numerical_energies(K, N, L, PREC)

    # Find error percentages of the approximated energies
    errors = (approx_energies/actual_energies - 1)*100

    # Print detailed information
    print()
    print(f"                 K = {K} | N = {N} | L = {L}               ")  
    print()
    print("|  State | Actual Energy | Approximate Energy |   Error %  |")
    print("|--------|---------------|--------------------|------------|")
    for i in range(K):
        print(f"| n = {i + 1:>2} | {actual_energies[i]:>13} "
              f"| {approx_energies[i]:>18.8f} | {errors[i]:>10.5f} |")
    print()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Harmonic Approximator",
                                     description=("Approximates the K first "
                                                  "energy levels of the "
                                                  "Harmonic Oscillator (H.O.) "
                                                  "using the N >= K first "
                                                  "eigenfunctions of the "
                                                  "Infinite Square Well "
                                                  "(I.S.W.) of length L."))

    parser.add_argument('-K', help="K first energies of H.O. to approximate.",
                        dest='K', 
                        type=int,
                        required=True)

    parser.add_argument('-N', help=("N first eiegenfunctions of I.S.W. to "
                                    "use for the approximation."),
                        dest='N', 
                        type=int,
                        required=True)
    
    parser.add_argument('-L', help=("Length of the I.S.W."),
                        dest='L', 
                        type=float,
                        required=True)
    
    parser.add_argument('-P', help=("Amount of points used in numerical "
                                    "integration between 0 and L."),
                        dest="PREC", 
                        type=int,
                        required=True)

    args = parser.parse_args()

    calculate_and_print_energies(args.K, args.N, args.L, args.PREC)