#!/bin/bash
#source /opt/intel/oneapi/setvars.sh intel64
for a in {1..21} 
do
mv Tem$a/BUFFER/VTOTAL.DAT DATA/V$a.dat
done