#! /bin/bash

s="$(<tr2.txt)"
rm -f  1a1xA.pssm input.txt tr2.txt tr3.txt
cd ${s}
rm -f  1a1xA.pssm 1a1xA.seq 1a1xA.hsa2 1a1xA.hsb2  1a1xA.spd3;
