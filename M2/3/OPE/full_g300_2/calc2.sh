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
L1=19
L2=11.6
L3=28.2
L=28.8
for L4 in `seq 11.2 0.05 50`
do  

    L42=`echo "scale=3;${L4} + ${L4}" | bc`
    L5=`echo "scale=3;${L4} - 11.2" | bc`
    L6=`echo "scale=3;${L4} + 11.2" | bc`
    echo -n "${L42}" >> result2.dat

    ./make_cfg ${dz} ${g1} ${w} 4 dc ${L} 7 dc ${L5} 9 dc ${L6} 8 dc ${L} 7
    #./make_cfg ${dz} ${g1} ${w} 4 dc ${L1} 8 dc ${L2} 7 dc ${L2} 7 dc ${L1} 9 

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> result2.dat
    done
    echo >> result2.dat
done
