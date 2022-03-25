#! /bin/bash

gawk '{if(NR!=1)if(NR!=2)if(NR!=3)if(NR!=4)if(NR!=5)if(NR!=6)if(NR!=7){print}}' properties.txt > tr4.txt #get path for MoRFchibi file 
s2="$(<tr4.txt)" 
gawk '{if(NR!=1)if(NR!=2)if(NR!=3)if(NR!=4)if(NR!=5)if(NR!=7)if(NR!=8){print}}' properties.txt > tr5.txt #get path for PROMIS file for later use
gawk {print} input.fasta > ${s2}/input.fasta 
cd ${s2}
./mc
gawk '{if(NR!=1)if(NR!=2)if(NR!=3)if(NR!=4)if(NR!=5)if(NR!=6)if(NR!=7)if(NR!=8)if(NR!=9)if(NR!=10)if(NR!=11)if(NR!=12)if(NR!=13)if(NR!=14){print}}' output.txt > output1.txt

