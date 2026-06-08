#!/bin/bash
cd exam
rm *.o
rm *.out
rm qm2p1
rm BUFFER/*
make
rm *.o
./qm2p1
#python3 mf_T.py