#!/bin/sh
for wl in `seq $1 $2 $3`
do
    echo -n > input_${wl}.pre
    echo -e "51" >> ./input_${wl}.pre

    for w in `seq 0.400 0.002 0.440`
    #for g in 0.150 0.160
    do
        cd beo_c/${wl}/${w}
        echo -e -n "${w}\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        echo -e -n "\t`tail -n2 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..

        cd b1/${wl}/${w}
        echo -e -n "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..

        cd b2_c/${wl}/${w}
        echo -e "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..
    done

    echo -e "51" >> ./input_${wl}.pre
    for w in `seq 0.400 0.002 0.440`
    #for g in 0.150 0.160
    do
        cd beo_a/${wl}/${w}
        echo -e -n "${w}\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        echo -e -n "\t`tail -n2 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..

        cd b1/${wl}/${w}
        echo -e -n "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..

        cd b2_a/${wl}/${w}
        echo -e "\t`tail -n1 _results.dat | head -n1 | cut -f4 `" >> ../../../input_${wl}.pre
        cd ../../..
    done

done