#!/bin/bash

CMD="./eblock"

MATRIX_DIR="../0_florida"
BLOCK_DIR="../0_blocks/Leiden"

GAMMA=(05 075 1 125 15 2)
COLS=(1 2 3)

for file in $MATRIX_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for g in "${GAMMA[@]}"; do
            for c in "${COLS[@]}"; do
                BLK_FILE="$(basename "$file" .mtx)_leiden_cpm_r0p00${g}_c${c}.blk"
                COL_FILE="$(basename "$file" .mtx)_leiden_cpm_r0p00${g}_c${c}.bcol"

                echo; echo "-------------------------------------------------"
                echo $mat "gamma = $g, coloring = $c"

                $CMD $file $BLOCK_DIR/$BLK_FILE $BLOCK_DIR/$COL_FILE
            done
        done
    fi
done
