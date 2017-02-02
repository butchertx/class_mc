#!/bin/bash
gcc -c -I../../lib ../../lib/random.c
g++ -c -std=c++11 -I../../lib ../../lib/MemTimeTester.cpp
g++ -c -std=c++11 -I../../lib ../../lib/Matrix.cpp
g++ -c -std=c++11 -I../../lib ../../lib/class_mc_io.cpp
g++ -c -std=c++11 -I../../lib ../../lib/IsingLattice2D.cpp
g++ -c -std=c++11 -I../../lib ../../lib/LongRangeWolff2D.cpp
g++ -c -std=c++11 -I../../lib spin_boson_class_mc.cpp -fopenmp
g++ -std=c++11 -o spin_boson.exe spin_boson_class_mc.o LongRangeWolff2D.o IsingLattice2D.o class_mc_io.o Matrix.o MemTimeTester.o random.o -fopenmp
