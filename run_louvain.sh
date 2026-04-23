#!/bin/bash

export OMP_NUM_THREADS=32

LOUV="./louvain"
TARGET_DIR="../0_florida"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        $LOUV $file
    fi
done
