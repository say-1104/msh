#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi

dz=0.01
g1=0.180
g2=1.180
w=0.410
Lsbend=10.0
echo -n > result.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> result.dat
L1=14.4
L2=8.0
for L in `seq 0.0 0.1 30.0`
do  
    echo -n "${L}" >> result.dat

    Ldel=`echo "scale=3;30.0 - ${L}" | bc`
    ./make_cfg ${dz} ${g2} ${w} 10 sbend ${Lsbend} ${g1} 3 dc ${L1} 8 dc ${L2} 7 sbend ${Lsbend} ${g2} 3 dc ${L} 9 dc ${Ldel} 7 sbend ${Lsbend} ${g1} 3 dc ${L2} 7 dc ${L1} 9 sbend ${Lsbend} ${g2} 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result.dat
    done
    echo >> result.dat
done
