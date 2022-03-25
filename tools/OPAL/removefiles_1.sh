#! /bin/bash

s="$(<tr4.txt)"
rm -f  tr4.txt tr5.txt;
rm -f input.fasta
rm -f list.txt
cd ${s}
rm -f  input.fasta output.txt output1.txt ;
