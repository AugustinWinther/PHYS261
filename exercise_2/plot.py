"""Radial Wavefunction Plotter

Plots the values of the Radial Wavefunction for given n, L, and Z values.
"""

# Standard library
import argparse

# Third party
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integ

# Local
import calc

def plot(rmax: float, prec: int,
         n: int, L: int, Z: int, ) -> None:
    """Plots the values of the Radial Wavefunction for given n, L, and Z values.
    
    """

    R, r = calc.R_values(rmax, prec, n, L, Z)
    u = r*R
    P = np.abs(u)**2

    # Print integral of probability density, expect. value, and most prob value
    print()
    print(f"For r ∈ [0, {rmax}] we have:")
    print()
    print(f"        ∫P = {integ.simpson(P, r):.15f}")
    print(f"       <r> = {integ.simpson(r*P, r):.15f} [1/Z]")
    print(f"MaxProb{{r}} = {r[np.argmax(P)]:.15f} [1/Z]")
    print()

    fig = plt.figure(layout="constrained")
    fig.suptitle(f"$n = {n}$ and $L = {L}$")

    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])

    ax1.set_title("Radial Wave Function")
    ax2.set_title("Reduced Radial Wave Function")
    ax3.set_title("Radial Probability Density")

    ax1.set_xlim([0, np.max(r)])
    ax2.set_xlim([0, np.max(r)])
    ax3.set_xlim([0, np.max(r)])

    ax1.set_ylim([np.min(R), np.max(R)])
    ax2.set_ylim([np.min(u), np.max(u)])
    ax3.set_ylim([np.min(P), np.max(P)])

    ax1.set_ylabel(f"$R_{{{n}{L}}}$(r)")
    ax2.set_ylabel(f"$u_{{{n}{L}}}$(r)")
    ax3.set_ylabel(f"$P_{{{n}{L}}}$(r)")

    ax3.set_xlabel("r $[a_0 / Z]$")

    ax1.plot(r, R)
    ax2.plot(r, u)
    ax3.plot(r, P)

    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Radial Wavefunction Plotter",
                                     description=("Plots the values of the"
                                                  "Radial Wavefunction for "
                                                  "given n, L, and Z values."))

    parser.add_argument("-n", help="Quantum number n.",
                        dest="n", 
                        type=int,
                        required=True)

    parser.add_argument("-L", help=("Quantum number L (same as small L)."),
                        dest="L", 
                        type=int,
                        required=True)
    
    parser.add_argument("-Z", help=("Amount of protons."),
                        dest="Z", 
                        type=int,
                        required=True)
    
    parser.add_argument("--rmax", help=("Maximum r value"),
                        dest="rmax", 
                        type=float,
                        default=6)
    
    parser.add_argument("--prec", help=("Amount of points used in numerical "
                                        "integration between 0 and rmax."),
                        dest="prec", 
                        type=int,
                        default=200)

    args = parser.parse_args()

    plot(args.rmax, args.prec, args.n, args.L, args.Z)