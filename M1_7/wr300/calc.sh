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
echo -n > heatmap_wr300.dat
for g in `seq 0.200 0.010 0.200`
do
    echo -n "\t${g}" >> heatmap_wr300.dat
done
echo >> heatmap_wr300.dat

for L in `seq 24.8 0.1 24.8`
do
    echo -n "${L}" >> heatmap_wr300.dat
    for g in `seq 0.200 0.010 0.200`
    do
        echo "L: ${L}\tdw: ${g}"
        echo -n > forMD.dat

        ./make_cfg ${dz} ${g} 1 ${L} ${g} 1 ${L} 1

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            cp output ./output_${wl}
        done
    done
    echo >> heatmap_wr300.dat
done