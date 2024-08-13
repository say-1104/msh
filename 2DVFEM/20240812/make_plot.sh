#!/bin/sh
nsio2=1.4443882489

echo -n > curve.dat
echo -e -n "wide[um]\tneff\tpol-X" >> curve.dat
cd result
for wide in `seq 0.200 0.005 0.500`
do  
    rownum=$((`cat ${wide}_results.dat | wc -l` - 1))
    for n in `seq ${rownum}`

    echo -e -n "${wide}" >> ../curve.dat
    neff=`tail -n${n} ${wide}_results.dat | head -n1 | cut -f2`
    if [ `echo "${neff} > ${nsio2}" | bc` == 1 ]; then
        echo -e -n "\t${neff}\t`tail -n${n} ${wide}_results.dat | head -n1 | cut -f15`" >> ../curve.dat
    fi
done