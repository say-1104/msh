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
g++ -o  calc_MD0 calc_MD0.cpp -std=c++11
if [ ! -e ./calc_MD0 ]; then
	echo "compile error!"
	exit
fi
dz=0.1
g=0.300
echo -n > heatmap.dat

for L1 in `seq 10.0 0.5 40.0`
do  
    for L2 in `seq 10.0 0.5 40.0`
    do  
        echo -n "${L1}\t${L2}\t" >> heatmap.dat   
        echo -n > forMD.dat
        for Ldel in `seq 10.0 0.2 40.0`
        do
            echo -n > forMD50.dat

            Lall=`echo "scale=3;${L1} + ${Ldel} + ${L2}" | bc`
            Lhalf=`echo "scale=3;${L1} + ${Ldel}" | bc`

            ./make_cfg ${dz} ${g} 1 ${Lall} ${g} 3 ${L1} 2 ${Lhalf} 4 ${Lall} 2

            for wl in `seq 1.530 0.010 1.570`
            do
                ./cmt -pcm2 -wl ${wl}
                echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
            done
            ./calc_MD50 <forMD50.dat
            echo "${Ldel}\t`tail -n1 MDlambda50 | head -n1 | cut -f1 `" >> forMD.dat
        done
        ./calc_MD0 <forMD.dat
        echo "`tail -n1 MDlambda | head -n1 | cut -f2 `\t`tail -n1 MDlambda | head -n1 | cut -f1 `" >> heatmap.dat
    done
done