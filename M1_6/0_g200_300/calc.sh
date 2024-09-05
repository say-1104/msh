#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_MD calc_MD.cpp -std=c++11
if [ ! -e ./calc_MD ]; then
	echo "compile error!"
	exit
fi
dz=0.1
echo -n > heatmap.dat
for L in `seq 20 0.1 180`
do
    for g in `seq 0.200 0.001 0.300`
    do
        echo "L: ${L}\tdw: ${g}"
        echo -n > forMD.dat

        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 2

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        
        ./calc_MD <forMD.dat
        echo "${L}\t${g}\t`head -n1 MDlambda`" >> heatmap.dat
    done

done