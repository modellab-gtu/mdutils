#!/bin/bash
# ABFE / hydration free energy analysis of an abfe_setup folder (automd/fe/abfe/analyze.py);
# options: abfe_analyze.sh -h
ROOT="$(cd "$(dirname "$(readlink -f "$0")")/.." && pwd)"
export PYTHONPATH="$ROOT" PYTHONNOUSERSITE=1 MPLBACKEND=Agg
# the ligparam env's scipy needs its own libstdc++ (CXXABI_1.3.15)
export LD_LIBRARY_PATH=/cta/capps/envs/ligparam/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
exec "${MDUTILS_FE_PYTHON:-/cta/capps/envs/ligparam/bin/python}" -m automd.fe.abfe.analyze "$@"
