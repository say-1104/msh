#!/bin/sh
export OMP_NUM_THREADS=12
export MKL_NUM_THREADS=12
wl=1.550
icc -o make_lps make_lps.cpp
for g in `seq 0.100 0.010 0.300`
do
    echo -n > Lc.dat
    ./make_lps 0.420 0.450 0.322 $g 0.020
    


for wl in `seq 1.530 0.010 1.570`
do
    echo -n > input_${wl}.pre
    echo -e "1" >> ./input_${wl}.pre

    cd beo_c/wave_${wl}
    echo -e -n "${w}\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    echo -e -n "\t`tail -n2 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

    cd b1/wave_${wl}
    echo -e -n "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

    cd b2_c/wave_${wl}
    echo -e "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

    echo -e "1" >> ./input_${wl}.pre

    cd beo_a/wave_${wl}
    echo -e -n "${w}\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    echo -e -n "\t`tail -n2 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

    cd b1/wave_${wl}
    echo -e -n "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

    cd b2_a/wave_${wl}
    echo -e "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../input_${wl}.pre
    cd ../..

done