#!/bin/bash

export OMP_NUM_THREADS=32

CMD="./leiden_cpm"
RES="0.001"
COL="1"


TARGET_DIR="../1_florida"

for file in $TARGET_DIR/*.mtx; do
    for col in $COL; do
        mat=$(basename "$file" .mtx)
        if [ -f "$file" ]; then
            $CMD $RES $file $col
        fi
    done
done
