#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD50 calc_MD50.cpp -std=c++11
if [ ! -e ./calc_MD50 ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.300
echo -n > forMD.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> forMD.dat
L1=15
L2=39
Ldel1=14.2

for Ldel in `seq 0.0 0.1 21.0`
do
    echo -n > forMD50.dat
    echo -n "${Ldel}" >> forMD.dat

    Ldel2=`echo "scale=3;42.0 - ${Ldel}" | bc`
    l1=`echo "scale=3;${Ldel} + ${Ldel}" | bc`
    l2=`echo "scale=3;${l1} + ${L2}" | bc`
    l3=`echo "scale=3;${l2} + ${Ldel}" | bc`
    l4=`echo "scale=3;${l3} + ${Ldel2}" | bc`
    l5=`echo "scale=3;${l4} + ${L2}" | bc`
    l6=`echo "scale=3;${l5} + ${Ldel1}" | bc`
    l7=`echo "scale=3;${l6} + ${L1}" | bc`


    ./make_cfg ${dz} ${g} 1 ${Ldel2} ${g} 2 ${L1} 3 ${l1} 4 

    for wl in `seq 1.530 0.010 1.570`
    do
        ./cmt -pcm2 -wl ${wl}
        echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
    done
    echo >> forMD.dat
done