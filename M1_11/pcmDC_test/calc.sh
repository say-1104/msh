#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.200
echo -n > result.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> result.dat

for L in `seq 0.0 0.1 30.0`
do
    echo -n "${L}" >> result.dat

    ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 7

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_taper/cmt -pcm9 -wl ${wl}
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result.dat
    done
    echo >> result.dat
done