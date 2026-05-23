#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-g++}"
CXXFLAGS=(
  -O2
  -std=c++17
  -I.
)

mkdir -p tests/bin

"${CXX}" "${CXXFLAGS[@]}" tests/test_cpm_properties.cpp -o tests/bin/test_cpm_properties
"${CXX}" "${CXXFLAGS[@]}" tests/test_cpm_leiden.cpp -o tests/bin/test_cpm_leiden

tests/bin/test_cpm_properties
tests/bin/test_cpm_leiden
