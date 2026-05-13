# Defect States in Discrete Schrodinger Chains

This repository contains the Python code for a Math 420 final project on defect states in finite one-dimensional discrete Schrodinger / tight-binding chains. The project studies how a single on-site defect produces an out-of-band eigenvalue and a geometrically localized eigenvector.

## What the script does

- Builds a finite tridiagonal Toeplitz chain with a single rank-one diagonal defect.
- Computes the defect eigenvalue for bulk and near-edge defect locations.
- Plots the defect eigenvalue as a function of defect strength alpha.
- Computes the localized defect eigenvector for a representative bulk defect.
- Plots the semi-log localization profile with the geometric envelope bound.
- Saves all generated figures to the `figures/` directory.

## Repository structure

```text
discrete-schrodinger-defect-states/
├── README.md
├── defect_states_schrodinger_chains.py
├── requirements.txt
├── figures/
│   ├── fig1_defect_vs_alpha.pdf
│   ├── fig1_defect_vs_alpha.png
│   ├── fig2_semilog_env.pdf
│   └── fig2_semilog_env.png
└── report/
    └── final_project_report.pdf
```

## Setup

Create a virtual environment, then install the dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell, activate the environment with:

```powershell
.venv\Scripts\Activate.ps1
```

## Run

```bash
python defect_states_schrodinger_chains.py
```

The script writes these files to `figures/`:

```text
figures/fig1_defect_vs_alpha.pdf
figures/fig1_defect_vs_alpha.png
figures/fig2_semilog_env.pdf
figures/fig2_semilog_env.png
```

## Dependencies

The project uses NumPy for matrix construction and symmetric eigensolvers, and Matplotlib for plotting.
