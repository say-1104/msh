#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
filename=result_si
dz=0.1
g=0.200
w=0.400
echo -n > ${filename}.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> ${filename}.dat

for L in `seq 0.0 0.1 80.0`
do
    echo -n "${L}" >> ${filename}.dat

    ./make_cfg ${dz} ${g} ${w} 1 dc ${L} 1 
    for wl in `seq 1.550 0.010 1.550`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl}
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> ${filename}.dat
    done
    echo >> ${filename}.dat
done
