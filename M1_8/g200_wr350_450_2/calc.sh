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
st=451
en=475
echo -n > heatmap100_g200_wr${st}_${en}.dat
echo -n > heatmap0_g200_wr${st}_${en}.dat
for L in `seq 10 0.1 110`
do
    for g in `seq 0.451 0.001 0.475`
    do
        echo "L: ${L}\tdw: ${g}"
        echo -n > forMD.dat

        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 2

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        
        ./calc_MD100 <forMD.dat
        echo "${L}\t${g}\t`head -n1 MDlambda`" >> heatmap100_g200_wr${st}_${en}.dat

        echo -n > forMD.dat
        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 1

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        
        ./calc_MD0 <forMD.dat
        echo "${L}\t${g}\t`head -n1 MDlambda`" >> heatmap0_g200_wr${st}_${en}.dat
    done

done