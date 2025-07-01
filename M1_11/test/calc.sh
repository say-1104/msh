#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g1=0.250
g2=0.300
w=0.400
Lsbend=10.0
echo -n > result.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> result.dat

for L in `seq 39.1 0.1 39.1`
do
    echo -n "${L}" >> result.dat

    ./make_cfg ${dz} ${g2} ${w} 3 sbend ${Lsbend} ${g1} 1 dc ${L} 1 sbend ${Lsbend} ${g2} 1

    for wl in `seq 1.570 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm4 -wl ${wl} -zw
        #echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result.dat
        echo 
    done
    echo >> result.dat
done