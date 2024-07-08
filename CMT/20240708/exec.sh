#!/bin/sh

for wl in `seq 1.530 0.010 1.570`
do
    echo -n > T_${wl}.dat
    for Leff in `seq 0.0 0.2 26.6`
    do
        ./cmt -pcm1 -wl ${wl} -leff ${leff} #2>output_err

        echo -e -n "${Leff}\t`tail -n1 output | head -n1 | cut -f3 `" >> ./T_${wl}.dat
        echo -e "\t`tail -n1 output | head -n1 | cut -f4 `" >> ./T_${wl}.dat
    done
done