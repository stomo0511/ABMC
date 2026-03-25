#!/usr/bin/sh

CMD="./abmc"
DAT="./0_florida/525825_parabolic_fem.mtx"
BLKS="32 64 128 256 512 1024 2048"
POL="1"

for bk in $BLKS
do
    echo "=== block size: $bk ==="
    $CMD $DAT $bk $POL
done