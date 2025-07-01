#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD50 calc_MD50.cpp -std=c++11
if [ ! -e ./calc_MD50 ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.300

L1=15
L2=39
Ldel1=14.2

echo -n > forMD50.dat

l1=`echo "scale=3;${L1} + ${Ldel1}" | bc`
l2=`echo "scale=3;${l1} + ${L2}" | bc`


./make_cfg ${dz} ${g} 1 ${l2} ${g} 3 ${L1} 2 ${l1} 4 ${l2} 2 

for wl in `seq 1.530 0.010 1.570`
do
    ./cmt -pcm2 -wl ${wl}
    echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
done
