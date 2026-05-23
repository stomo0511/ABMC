#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-g++}"
TARGET="${TARGET:-leiden_cpm}"

CXXFLAGS=(
  -O3
  -march=native
  -std=c++17
  -I.
)

LDFLAGS=()

uname_s="$(uname)"
if [[ "${uname_s}" == "Darwin" ]]; then
  CXXFLAGS+=(
    -I/opt/homebrew/include
    -Xpreprocessor
    -fopenmp
    -I/opt/homebrew/opt/libomp/include
  )
  LDFLAGS+=(
    -L/opt/homebrew/lib
    -L/opt/homebrew/opt/libomp/lib
    -lomp
  )
else
  CXXFLAGS+=(
    -fopenmp
  )
  LDFLAGS+=(
    -fopenmp
  )
fi

"${CXX}" "${CXXFLAGS[@]}" -o "${TARGET}" gve-leiden-cpm.cpp "${LDFLAGS[@]}"
