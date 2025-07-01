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

dz=0.01
g1=0.180
g2=1.180
w=0.410
Lsbend=10.0
Lall=22.4
echo -n > heatmap.dat

for L1 in `seq 0 0.1 ${Lall}`
do  

    echo -n > forMD50.dat
    L2=`echo "scale=3;${Lall} - ${L1}" | bc`
    ./make_cfg ${dz} ${g2} ${w} 4 sbend ${Lsbend} ${g1} 3 dc ${L1} 8 dc ${L2} 7 sbend ${Lsbend} ${g2} 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
    done
    ./calc_MD50 <forMD50.dat
    echo -n "${L1}\t${L2}\t`tail -n1 MDlambda50 | head -n1 | cut -f1 `" >> heatmap.dat
    echo "\t${Lall}" >> heatmap.dat 
done

