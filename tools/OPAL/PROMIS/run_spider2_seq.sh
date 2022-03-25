#! /bin/bash

cd ..
gawk '{if(NR!=1)if(NR!=3)if(NR!=4)if(NR!=5)if(NR!=6)if(NR!=7)if(NR!=8){print}}' properties.txt > PROMIS/tr2.txt #get path for spider2
gawk '{if(NR!=1)if(NR!=2)if(NR!=3)if(NR!=5)if(NR!=6)if(NR!=7)if(NR!=8){print}}' properties.txt > PROMIS/tr3.txt #get path for LIBsvm for later use
gawk {print} input.fasta > PROMIS/input.txt 
cd PROMIS
s2="$(<tr2.txt)"
gawk {print} input.txt  > ${s2}/1a1xA.seq 

cd ${s2}
rm -f  1a1xA.hsa2 1a1xA.hsb2  1a1xA.spd3;
bash run_local.sh 1a1xA.seq

