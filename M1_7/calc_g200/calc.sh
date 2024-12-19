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
echo -n > PSR.dat
for wl in `seq 1.530 0.010 1.570`
do
    echo  -n "\t${wl}" >> PSR.dat
done
echo >> PSR.dat

for L in `seq 0.0 0.1 34.6`
do
    for g in 0.436
    do
        echo "L: ${L}\tdw: ${g}"

        ./make_cfg ${dz} ${g} 1 34.6 ${g} 2 ${L} 1 34.6 2
        echo -n "${L}" >> PSR.dat

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> PSR.dat
        done
        echo >> PSR.dat
    done

done