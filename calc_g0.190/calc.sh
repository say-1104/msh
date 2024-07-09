#!/bin/sh
dz=0.1
z0=0.0
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
#make all
for wl in `seq 1.540 0.010 1.570`
do
    echo -n > PSR_${wl}.dat
    count=1
    for g in 0.190
    do
        z1=37.6
        ./make_cfg ${dz} ${z0} ${g} ${z1} ${g}
        Leff_fi=`tail -n1 structure.cfg | head -n1 | cut -f1 `
        for Leff in `seq 0.0 0.2 ${Leff_fi}`
        do
            ./cmt -pcm1 -wl ${wl} -leff ${Leff} 

            echo -e -n "${Leff}\t`tail -n1 output | head -n1 | cut -f3 `" >> PSR_${wl}.dat
            echo -e "\t`tail -n1 output | head -n1 | cut -f4 `" >> PSR_${wl}.dat
        done
    done
done