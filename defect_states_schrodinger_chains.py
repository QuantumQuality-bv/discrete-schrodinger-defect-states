#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Numerical experiments for a finite one-dimensional discrete Schrodinger / tight-binding chain
with a single on-site defect.

The script generates the two project figures:
1. the out-of-band defect eigenvalue as a function of defect strength, and
2. the semi-log localization profile of the associated defect eigenvector.

All figures are written to the local figures/ directory.
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as npl
import matplotlib.ticker as mticker  # for inset tick control

OUTPUT_DIR = Path(__file__).resolve().parent / "figures"

# Helper: finite-N cosine band edges

def band_edges_finite_N(a: float, b: float, N: int):
    """
    Finite-N cosine band edges for the Dirichlet chain:
      E_max^(N) = a + 2|b| cos(π/(N+1)),   E_min^(N) = a - 2|b| cos(π/(N+1)).
    """
    t = np.cos(np.pi / (N + 1.0))
    e_max = a + 2.0 * abs(b) * t
    e_min = a - 2.0 * abs(b) * t
    return e_min, e_max

# Figure 1 — Defect level vs α (bulk and near-edge sites)

def figure1():
    a, b = 0.0, 1.0
    N = 200
    ks = [100, 10]   # bulk site and near-edge site
    alpha_vals = np.arange(0.1, 8.05, 0.05)

    # Finite-N band edges
    eminN, emaxN = band_edges_finite_N(a, b, N)

    # Analytic upper bound (finite-N)
    U_bound = emaxN + alpha_vals

    # Base Toeplitz chain (Dirichlet)
    H0 = np.diag(np.full(N, a))
    idx = np.arange(N - 1)
    H0[idx, idx + 1] = b
    H0[idx + 1, idx] = b

    # For each α and k, build H(α,k) and take the largest eigenvalue
    lambdas_by_k = {k: [] for k in ks}

    for alpha in alpha_vals:
        for k in ks:
            H = H0.copy()
            H[k - 1, k - 1] += alpha            # rank-one defect at site k
            w = npl.eigvalsh(H)                 # symmetric eigensolver
            lambdas_by_k[k].append(w[-1])       # outlier = largest eigenvalue

    lam_bulk = np.asarray(lambdas_by_k[100])
    lam_edge = np.asarray(lambdas_by_k[10])

    # Plot
    fig, ax = plt.subplots(figsize=(6.5, 4.0))

    # Solid blue (thicker): bulk k=100
    ax.plot(
        alpha_vals, lam_bulk,
        color="C0", lw=2.2, ls="-", zorder=3,
        label=r"$\lambda(\alpha,100)$ (bulk)"
    )
    # Orange dashed (slightly thinner): near-edge k=10
    ax.plot(
        alpha_vals, lam_edge,
        color="C1", lw=1.6, ls="--", dashes=(6, 2), zorder=3.2,
        label=r"$\lambda(\alpha,10)$ (near edge)"
    )


    # Dashed curve: provable finite-N upper bound
    ax.plot(alpha_vals, U_bound, 'k--', lw=1.0, zorder=1.5,
            label=r"$U_{\mathrm{bound}}(\alpha)=E_{\max}^{(N)}+\alpha$")

    # Dotted finite-N band edges (upper labeled for legend; lower dimmed)
    ax.axhline(eminN, color='0.7', ls=':', lw=0.8, zorder=0.4)  # lower edge (dim)
    ax.axhline(emaxN, color='k', ls=':', lw=1.0, zorder=0.5,
               label=r"finite-$N$ band edges")

    # Axis limits: show only relevant range
    ax.set_xlim(alpha_vals.min(), alpha_vals.max())
    ax.set_ylim(eminN - 0.05, U_bound.max() + 0.2)

    ax.set_xlabel(r"$\alpha$")
    ax.set_ylabel(r"$\lambda(\alpha,k)$")
    ax.set_title(r"Defect eigenvalue $\lambda(\alpha,k)$ vs. $\alpha$ ($a=0,b=1,N=200$)")
    ax.legend(loc='upper left', fontsize=8)

    # -------- Inset to reveal tiny bulk vs near-edge difference near band edge
    # Focus on small α where separation is smallest
    mask = (alpha_vals >= 0.2) & (alpha_vals <= 1.5)
    ax_in = ax.inset_axes([0.57, 0.15, 0.38, 0.38])  # [x0,y0,w,h] in axes coords
    ax_in.plot(alpha_vals[mask], lam_bulk[mask], color="C0", lw=1.9, ls="-", zorder=3)
    ax_in.plot(alpha_vals[mask], lam_edge[mask], color="C1", lw=1.4, ls="--", zorder=3.2)
    # Band edge in inset
    ax_in.axhline(emaxN, color='k', ls=':', lw=0.8, zorder=0.5)
    # Tight y-limits around E_max^(N)
    y_lo = emaxN - 5e-4
    y_hi = max(lam_bulk[mask].max(), lam_edge[mask].max()) + 5e-3
    ax_in.set_xlim(alpha_vals[mask].min(), alpha_vals[mask].max())
    ax_in.set_ylim(y_lo, y_hi)
    ax_in.set_xticks([0.2, 0.8, 1.4])
    # Show only the band-edge label on the inset y-axis with small font
    ax_in.yaxis.set_major_locator(mticker.FixedLocator([emaxN]))
    ax_in.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.4f'))
    ax_in.tick_params(axis='both', labelsize=7)
    # Label the inset band edge explicitly
    x0, x1 = alpha_vals[mask].min(), alpha_vals[mask].max()
    ax_in.text(x0 + 0.03*(x1 - x0), emaxN + 0.0001*(y_hi - y_lo),
               r"$E_{\max}^{(N)}$", fontsize=7, ha='left', va='bottom')
    ax_in.set_title("zoom near band edge", fontsize=8)

    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = OUTPUT_DIR / "fig1_defect_vs_alpha.pdf"
    png_path = OUTPUT_DIR / "fig1_defect_vs_alpha.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path

# Figure 2 — Semi-log localization and envelope for a bulk defect

def figure2():
    """
    Semi-log plot of the bulk defect eigenvector corresponding to the
    defect eigenvalue λ(α,k) together with the localization envelope of
    Theorem 3, for a=0, b=1, N=200, k=100, α=2; here λ(α,k) is the unique
    out-of-band eigenvalue of H(α,k) lying above the surrogate band
    [a-2|b|, a+2|b|].
    """
    a, b = 0.0, 1.0
    N = 200
    k = 100
    alpha = 2.0

    # Base Toeplitz chain (Dirichlet)
    H0 = np.diag(np.full(N, a))
    idx = np.arange(N - 1)
    H0[idx, idx + 1] = b
    H0[idx + 1, idx] = b

    # Add defect at site k
    H = H0.copy()
    H[k - 1, k - 1] += alpha

    # Eigenvalues/eigenvectors
    w, V = npl.eigh(H)
    eminN, emaxN = band_edges_finite_N(a, b, N)

    # Select outlier strictly above the finite-N band edge
    mask = w > emaxN + 1e-9
    if np.any(mask):
        j_out = np.where(mask)[0][0]
    else:
        # Fallback: largest eigenvalue
        j_out = np.argmax(w)

    lam = w[j_out]
    v = V[:, j_out].real
    v /= npl.norm(v)

    # Spatial index relative to defect site
    j = np.arange(1, N + 1)
    x = np.abs(j - k)
    y = np.abs(v)

    # Quadratic for r_-
    #    b r^2 + (a - λ) r + b = 0
    disc = (a - lam) ** 2 - 4.0 * (b ** 2)
    if disc < 0.0:
        disc = 0.0  # safeguard against tiny negative due to rounding
    root = np.sqrt(disc)
    r1 = (-(a - lam) + root) / (2.0 * b)
    r2 = (-(a - lam) - root) / (2.0 * b)
    r_minus = r1 if abs(r1) < 1.0 else r2
    r_abs = abs(r_minus)

    d = min(k - 1, N - k)

    # Envelope (only if |r_-| < 1)
    envelope_available = r_abs < 1.0 - 1e-12
    if envelope_available:
        C = 1.0 / np.sqrt(
            (1.0 - r_abs ** 2) * (1.0 - r_abs ** (2 * d + 2))
        )
        y_env = C * (r_abs ** x)

    # Plot (semi-log)
    fig, ax = plt.subplots(figsize=(6.2, 4.0))

    # Restrict x-range for visibility
    mask_x = x <= 60
    ax.semilogy(x[mask_x], y[mask_x], 'o', ms=3,
                label=r"$|v_j^{(\alpha,k)}|$ (bulk defect mode)")

    if envelope_available:
        ax.semilogy(x[mask_x], y_env[mask_x], 'k-', lw=1.4,
                    label="envelope bound")
    else:
        ax.text(
            0.05, 0.9,
            r"$|r_-|\ge 1$: no exponential envelope",
            transform=ax.transAxes,
            fontsize=8,
            ha="left",
            va="center",
        )

    ax.set_xlim(0, 60)
    ax.set_xlabel(r"$|j-k|$")
    ax.set_ylabel(r"$|v_j|$")
    ax.set_title(
        r"Semi-log localization at defect eigenvalue $\lambda(\alpha,k)$" "\n"
        r"($a=0,b=1,N=200,k=100,\alpha=2$)",
        fontsize=11,
    )
    ax.legend(fontsize=8)

    fig.tight_layout()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = OUTPUT_DIR / "fig2_semilog_env.pdf"
    png_path = OUTPUT_DIR / "fig2_semilog_env.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight")
    plt.close(fig)
    return pdf_path, png_path


def main():
    figure1_paths = figure1()
    figure2_paths = figure2()

    print("Generated figures:")
    for path in (*figure1_paths, *figure2_paths):
        print(f"  {path}")


if __name__ == "__main__":
    main()
