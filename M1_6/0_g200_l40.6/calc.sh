#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
dz=0.1
wcent=0.450
echo -n > heatmap.dat
for L in `seq 40.6 0.1 40.6`
do
    #echo -n "${L}" >> graph.dat
    for dw in `seq 0.000 0.001 0.000`
    do

        wst=`echo "scale=3;${wcent} - ${dw}" | bc`
        wfi=`echo "scale=3;${wcent} + ${dw}" | bc`
        for L1 in `seq 0.0 0.1 40.6`
        do
            
            for L2 in `seq 0.0 0.1 ${L1}`
            do
                echo "${L1}\t${L2}\t1.0" >> heatmap.dat
            done
            L1end=`echo "scale=1;${L1} + 0.1" | bc`
            for L2 in `seq ${L1end} 0.1 40.6`
            do  
                echo -n > forMD.dat
                #./make_cfg ${dz} ${wst} 1 ${L} ${wfi} 3 22.0 2 70.0 1 ${L} 2
                ./make_cfg ${dz} ${wst} 1 ${L} ${wfi} 3 ${L1} 2 ${L2} 1 ${L} 2
                
                for wl in `seq 1.530 0.010 1.570`
                do
                    ./cmt -pcm1 -wl ${wl}
                    echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD.dat
                done
            
                ./calc_MD <forMD.dat
                echo "${L1}\t${L2}\t`head -n1 MDlambda`" >> heatmap.dat
            done
        done
    done
done