#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1

for g in `seq 0.200 0.1 1.200`
do
    echo -n > result_${g}.dat
    echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> result_${g}.dat

    for L in `seq 0.0 0.1 30.0`
    do
        echo -n "${L}" >> result_${g}.dat

        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 2

        for wl in `seq 1.530 0.010 1.570`
        do
            ~/msh/CMT_taper/cmt -pcm4 -wl ${wl}
            echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result_${g}.dat
        done
        echo >> result_${g}.dat
    done
done