#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi

dz=0.05
g1=0.300

w=0.400

echo -n > result2.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> result2.dat

for L in `seq 11.1 0.1 50.0`
do  
    echo -n "${L}" >> result2.dat
    L1=`echo "scale=3;${L} - 11.1" | bc`
    L2=`echo "scale=3;${L} + 11.1" | bc`
    ./make_cfg ${dz} ${g1} ${w} 2 dc ${L1} 8 dc ${L2} 9 

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result2.dat
    done
    echo >> result2.dat
done
