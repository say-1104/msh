#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi

#DIR=$(cd $(dirname $0);pwd)
w=0.400
p=0.300

dz=0.01
g=0.200
Lsbend=10.0
filename="result_ca"
echo -n > ${filename}.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> ${filename}.dat
for L in `seq 0.0 0.1 25.0`
do  
    echo -n "${L}" >> ${filename}.dat

    ./make_cfg ${dz} ${g} ${w} 1 dc ${L} 8

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> ${filename}.dat
    done
    echo >> ${filename}.dat
done