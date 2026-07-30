#!/bin/bash

export OMP_NUM_THREADS=32

# 実行回数
N=5

CMD="./abmc"

TARGET_DIR="../1_florida"
BLKS="32 64 128 256 512 1024"
POLS="1"
COLS="1"

for file in $TARGET_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for bs in $BLKS; do
            for p in $POLS; do
                for c in $COLS; do
                    for i in $(seq 1 $N); do
                        echo "$mat"
                        $CMD $file $bs $p $c
                    done
                done
            done
        done
    fi
done
