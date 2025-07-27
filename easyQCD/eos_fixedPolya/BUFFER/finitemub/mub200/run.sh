#!/bin/bash
source /opt/intel/oneapi/setvars.sh intel64
rm *.mod
rm *.o
rm *.out
rm eqcd
cd BUFFER
rm *.DAT
rm *.dat
cd ..
make
rm *.o
rm *.mod
./eqcd
