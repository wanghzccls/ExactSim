EXSim: SFMT.c main.cpp 
	g++ -march=core2 -pthread -std=c++11 -DSFMT_MEXP=607 -I SFMT-src-1.4.1/ -O3 -o EXSim  SFMT.c  main.cpp
