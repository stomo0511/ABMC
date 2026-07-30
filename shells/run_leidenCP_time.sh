#!/bin/bash

export OMP_NUM_THREADS=32

# 実行回数
N=5

CMD="./leiden_cpm"

TARGET_DIR="../1_florida"
COLS="1"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for c in $COLS; do
            for i in $(seq 1 $N); do
                echo "$mat"
                $CMD 0.001 $file $p $c
            done
        done
    fi
done
