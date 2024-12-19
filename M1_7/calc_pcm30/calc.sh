#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD100 calc_MD100.cpp -std=c++11
if [ ! -e ./calc_MD100 ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
echo -n > PSR2.dat
for wl in `seq 1.530 0.010 1.570`
do
    echo  -n "\t${wl}" >> PSR2.dat
done
echo >> PSR2.dat
wcent=0.500
for L in `seq 0.0 0.1 147.1`
do
    for dw in 0.009
    do
        echo "L: ${L}\tdw: ${dw}"

        wst=`echo "scale=3;${wcent} - ${dw}" | bc`
        wfi=`echo "scale=3;${wcent} + ${dw}" | bc`

        ./make_cfg ${dz} ${wst} 1 147.1 ${wfi} 2 ${L} 1 147.1 2
        echo -n "${L}" >> PSR2.dat

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> PSR2.dat
        done
        echo >> PSR2.dat
    done

done