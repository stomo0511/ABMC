#!/bin/bash

export OMP_NUM_THREADS=32

CMD="./leiden_cpm"
RES="0.0005 0.00075 0.001 0.00125 0.0015 0.002"
COL="2"

TARGET_DIR="../0_florida"

for file in $TARGET_DIR/*.mtx; do
    for col in $COL; do
    	for res in $RES; do
            mat=$(basename "$file" .mtx)
            if [ -f "$file" ]; then
                echo -n $mat ", "
        		$CMD $res $file $col
            fi
	    done
    done
done
