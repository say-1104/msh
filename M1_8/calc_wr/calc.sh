#!/bin/sh
g++ -o  calcwr calcwr.cpp -std=c++11
if [ ! -e ./calcwr ]; then
	echo "compile error!"
	exit
fi

./calcwr <neffsi_apcm.dat