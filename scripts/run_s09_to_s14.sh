#!/usr/bin/env bash
# Gene calling -> MMseqs2 families -> matrix -> baselines -> MLP -> transfer.
# Each step logs to outputs/<step>.log; stops at the first failure.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=/home/rl/mambaforge/envs/karl_datasci/bin/python
TAG=${TAG:-id50_c80}
step() { echo "[$(date +%H:%M:%S)] $1"; shift; "$@"; }
step s09 $PY scripts/s09_call_genes.py              > outputs/s09_call_genes.log 2>&1
step s10 $PY scripts/s10_mmseqs_cluster.py          > outputs/s10_mmseqs_${TAG}.log 2>&1
step s11 $PY scripts/s11_family_matrix.py --tag $TAG > outputs/s11_family_matrix_${TAG}.log 2>&1
step s12 $PY scripts/s12_baseline_models.py --tag $TAG > outputs/s12_${TAG}_genus.log 2>&1
step s13 $PY scripts/s13_multitask_mlp.py --tag $TAG --final > outputs/s13_${TAG}_genus.log 2>&1
step s14 $PY scripts/s14_transfer_lab_genomes.py --tag $TAG > outputs/s14_${TAG}.log 2>&1
echo "[$(date +%H:%M:%S)] all done"
