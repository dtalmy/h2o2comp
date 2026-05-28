#!/usr/bin/env bash
# Reproduce every figure in the paper.
# Run from the src/ directory:
#     bash run_all.sh

set -euo pipefail

cd "$(dirname "$0")"

echo "[1/4] intro_plot.py  (figure 1)"
python intro_plot.py

echo "[2/4] fit_kdam.py    (figures S1, S2, S3, kdam tables)"
python fit_kdam.py

echo "[3/4] simulations.py (figures 2, 3, 4, 5, 6)"
python simulations.py

echo "[4/4] sensitivity.py (figure S4)"
python sensitivity.py

echo "done. outputs in ../figures/ and ../data/"
