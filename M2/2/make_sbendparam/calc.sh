#!/bin/sh
g++ -o  make_sbend make_sbend.cpp -std=c++11
if [ ! -e ./make_sbend ]; then
	echo "compile error!"
	exit
fi

w=0.400
wpcm=0.300
gst=0.200
gfi=1.200
l=10.0
./make_sbend ${w} ${wpcm} ${gst} ${gfi} ${l}

