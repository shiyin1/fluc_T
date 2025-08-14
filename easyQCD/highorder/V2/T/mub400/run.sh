#!/bin/bash
#source /opt/intel/oneapi/setvars.sh intel64
rm exam/eqcd
rm exam/*.out
rm exam/*.o
rm exam/*.mod
rm exam/BUFFER/*
make
rm *.o
rm *.mod
./eqcd