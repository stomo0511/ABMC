#!/bin/bash

export OMP_NUM_THREADS=32

CMD="./abmc"
BLOCKS="32 64 128 256 512 1024"
POL="1 2 3"
COL="1 2 3"


TARGET_DIR="../1_florida"

for file in $TARGET_DIR/*.mtx; do
    for blk in $BLOCKS; do
        for pol in $POL; do
            for col in $COL; do
                mat=$(basename "$file" .mtx)
                if [ -f "$file" ]; then
                    $CMD $file $blk $pol $col
                fi
            done
        done
    done
done
