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
echo -n > forMD.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> forMD.dat

for Ldel in `seq 0.0 0.1 42.0`
do
    echo -n "${Ldel}" >> forMD.dat

    l1=`echo "scale=3;84.0 - ${Ldel}" | bc`

    ./make_cfg ${dz} ${g} 1 84.0 ${g} 2 ${l1} 4 84.0 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ./cmt -pcm2 -wl ${wl}
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
    done
    echo >> forMD.dat
done