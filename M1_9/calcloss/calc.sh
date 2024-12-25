#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD100 calc_MD100.cpp -std=c++11
if [ ! -e ./calc_MD100 ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
echo -n > Loss.dat

for g in `seq 0.400 0.010 0.400`
do  
    echo -n "${g}" >> Loss.dat   
    echo -n > forMD.dat
    for L in `seq 50.0 0.1 50.0`
    do
        echo -n "${L}" >> forMD.dat
        echo "L: ${L}\tdw: ${g}"

        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 2

        for wl in `seq 1.550 0.010 1.550`
        do
            ./cmt -pcm1 -wl ${wl}
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        ./calc_MD0 <forMD.dat
    done

done