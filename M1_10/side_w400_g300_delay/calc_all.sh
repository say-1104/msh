#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi

dz=0.1
g=0.300
echo -n > output_side_w400_g300.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> output_side_w400_g300.dat
L1=22.5
L2=38
Ldel1=12.5

for Ldel in `seq 0.0 0.1 17.5`
do
    echo -n "${Ldel}" >> output_side_w400_g300.dat

    Ldel2=`echo "scale=3;17.5 - ${Ldel}" | bc`
    l1=`echo "scale=3;${L1} + ${Ldel1}" | bc`
    l2=`echo "scale=3;${l1} + ${L2}" | bc`
    l3=`echo "scale=3;${l2} + ${Ldel}" | bc`
    l4=`echo "scale=3;${l3} + ${Ldel2}" | bc`
    l5=`echo "scale=3;${l4} + 17.5" | bc`
    l6=`echo "scale=3;${l5} + ${L2}" | bc`
    l7=`echo "scale=3;${l6} + ${Ldel1}" | bc`
    l8=`echo "scale=3;${l7} + ${L1}" | bc`


    ./make_cfg ${dz} ${g} 1 ${l8} ${g} 9 ${L1} 5 ${l1} 3 ${l2} 5 ${l3} 4 ${l4} 3 ${l5} 1 ${l6} 5 ${l7} 1 ${l8} 5

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT/cmt -pcm5 -wl ${wl}
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> output_side_w400_g300.dat
    done
    echo >> output_side_w400_g300.dat
done