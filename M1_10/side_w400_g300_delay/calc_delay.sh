#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.300
echo -n > side_delay_w400_g300.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> side_delay_w400_g300.dat

for Ldel in `seq 0.1 0.1 50.0`
do
    echo -n "${Ldel}" >> side_delay_w400_g300.dat

    l1=`echo "scale=3;${Ldel} + ${Ldel}" | bc`

    ./make_cfg ${dz} ${g} 1 ${l1} ${g} 2 ${Ldel} 1 ${l1} 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT/cmt -pcm5 -wl ${wl}
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> side_delay_w400_g300.dat
    done
    echo >> side_delay_w400_g300.dat
done