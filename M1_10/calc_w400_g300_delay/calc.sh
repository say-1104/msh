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


for Ldel in `seq 0.0 0.1 42.0`
do
    echo -n > forMD50.dat

    Lhalf=`echo "scale=3;84.0 - ${Ldel}" | bc`

    ./make_cfg ${dz} ${g} 1 84.0 ${g} 3 ${Ldel} 3 ${Lhalf} 4 84.0 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ./cmt -pcm2 -wl ${wl}
        echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
    done
    ./calc_MD50 <forMD50.dat
    echo "${Ldel}\t`tail -n1 MDlambda50 | head -n1 | cut -f1 `" >> forMD.dat
done