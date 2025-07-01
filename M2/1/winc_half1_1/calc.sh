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

dz=0.1
g1=0.200
g2=1.200
w=0.400
Lsbend=10.0
echo -n > heatmap.dat

for L1 in `seq 0.0 0.2 5.0`
do  
    for L2 in `seq 0.0 0.2 30.0`
    do  
        for Ldel in `seq 0.0 0.2 30.0`
        do  
            echo -n > forMD50.dat
            
            ./make_cfg ${dz} ${g1} ${w} 4 dc ${L1} 1 dc ${Ldel} 4 dc ${L2} 1 sbend ${Lsbend} ${g2} 1

            for wl in `seq 1.530 0.010 1.570`
            do
                ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl}
                echo "${wl}\t`tail -n1 output | head -n1 | cut -f3 `" >> forMD50.dat
            done
            rm MDlambda50
            ./calc_MD50 <forMD50.dat
            if [ -f "MDlambda50" ]; then
                echo -n "${L1}\t${L2}\t${Ldel}\t`tail -n1 MDlambda50 | head -n1 | cut -f1 `" >> heatmap.dat
                Lall=`echo "scale=3;${L1} + ${del} + ${L2}" | bc`
                echo "\t${Lall}" >> heatmap.dat
            fi
        done    
    done
done
