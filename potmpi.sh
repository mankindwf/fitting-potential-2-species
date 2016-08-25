#!/bin/bash
mpiCC -c -g -std=c++11 potential.cc -I../libinfer-master/ -fmax-errors=1
mpiCC potential.o -o potential -L../libinfer-master/ -linfer
