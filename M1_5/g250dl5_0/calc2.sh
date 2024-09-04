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
wcent=0.450
dL=5
echo -n > heatmap.dat
for L in `seq 20 0.1 160`
do
    for dw in `seq 0.000 0.001 0.050`
    do
        echo "L: ${L}\tdw: ${dw}"
        echo -n > forMD.dat

        wst=`echo "scale=3;${wcent} - ${dw}" | bc`
        wfi=`echo "scale=3;${wcent} + ${dw}" | bc`
        L2=`echo "scale=3;${L} - ${dL}" | bc`

        ./make_cfg ${dz} ${wst} 1 ${L} ${wfi} 3 ${dL} 1 ${L2} 2 ${L} 1

        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        
        ./calc_MD <forMD.dat
        echo "${L}\t${dw}\t`head -n1 MDlambda`" >> heatmap.dat
    done

done