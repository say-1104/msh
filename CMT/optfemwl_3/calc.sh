#!/bin/sh
dz=0.1
z0=0.0
w0=0.450
z4=37.6
w4=0.450
flag=1  #1: wlの最適化
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  calc_FOM calc_FOM.cpp -std=c++11
if [ ! -e ./calc_FOM ]; then
	echo "compile error!"
	exit
fi
#make all

echo -n > allFOM.dat
for z1 in `seq 0.2 0.2 37.4`
do
    for w1 in `seq 0.400 0.005 0.500`
    do  
        echo -e "z1: ${z1}\tw1: ${w1}"
        echo -n > PSR.dat
        echo -e "${flag}" >> PSR.dat
        ./make_cfg ${dz} ${z0} ${w0} ${z1} ${w1} ${z4} ${w4}
        for Leff in 0.0 37.6
        do
            for wl in `seq 1.530 0.010 1.570`
            do
                ./cmt -pcm1 -wl ${wl} -leff ${Leff}
                echo -e -n "${Leff}\t${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> PSR.dat
                echo -e "\t`tail -n1 output | head -n1 | cut -f4 `" >> PSR.dat
            done
        done
        ./calc_FOM <PSR.dat
        echo -e "${z1}\t${w1}\t`head -n1 FOM`" >> allFOM.dat
    done
done