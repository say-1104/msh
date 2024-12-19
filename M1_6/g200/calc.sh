#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1
wcent=0.450
echo -n > graph.dat
for L in `seq 38.1 0.1 38.1`
do
    #echo -n "${L}" >> graph.dat
    for dw in `seq 0.000 0.001 0.000`
    do
        echo -n > forMD.dat

        wst=`echo "scale=3;${wcent} - ${dw}" | bc`
        wfi=`echo "scale=3;${wcent} + ${dw}" | bc`

        #./make_cfg ${dz} ${wst} 1 ${L} ${wfi} 3 22.0 2 70.0 1 ${L} 2
        ./make_cfg ${dz} ${wst} 1 ${L} ${wfi} 3 4.0 1 8.0 2 ${L} 1
        
        for wl in `seq 1.530 0.010 1.570`
        do
            ./cmt -pcm1 -wl ${wl}
            cp output ${wl}_output
            echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
        done
        echo >> graph.dat
    done
done