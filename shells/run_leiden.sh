#!/bin/bash

export OMP_NUM_THREADS=32

CMD="./leiden"
RES="1.0"
COL="1 2 3"


TARGET_DIR="../0_florida"

for file in $TARGET_DIR/*.mtx; do
    for col in $COL; do
        mat=$(basename "$file" .mtx)
        if [ -f "$file" ]; then
            $CMD $RES $file $col
        fi
    done
done
