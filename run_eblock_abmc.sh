#!/bin/bash

CMD="./eblock"

MATRIX_DIR="../0_florida"
BLOCK_DIR="../0_blocks/ABMC"

BS=(32 64 128 256 512 1024)
POLYC=(1 2 3)
COLS=(1 2 3)

for file in $MATRIX_DIR/*.mtx; do
    mat=$(basename "$file" .mtx)
    if [ -f "$file" ]; then
        for bs in "${BS[@]}"; do
			for p in "${POLYC[@]}"; do
                for c in "${COLS[@]}"; do
                    BLK_FILE="$(basename "$file" .mtx)_abmc_B${bs}_p${p}_c${c}.blk"
                    COL_FILE="$(basename "$file" .mtx)_abmc_B${bs}_p${p}_c${c}.bcol"

                    echo; echo "-------------------------------------------------"
                    echo $mat "bs = $bs, policy = $p, coloring = $c"

                    $CMD $file $BLOCK_DIR/$BLK_FILE $BLOCK_DIR/$COL_FILE
                done
            done
        done
    fi
done
