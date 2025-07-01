#!/bin/sh
g++ -o  make_cfg make_cfg.cpp -std=c++11
if [ ! -e ./make_cfg ]; then
	echo "compile error!"
	exit
fi
g++ -o  forExcel forExcel.cpp -std=c++11
if [ ! -e ./forExcel ]; then
	echo "compile error!"
	exit
fi
DIR=$(cd $(dirname $0);pwd)
w=`echo ${DIR} | awk -F "/" '{ print $8 }'`
p=`echo ${DIR} | awk -F "/" '{ print $9 }'`
d=`echo ${DIR} | awk -F "/" '{ print $10 }'`

dz=0.01
g1=`echo "scale=3;0.200 - ${w} + 0.400" | bc`
g2=`echo "scale=3;1.200 - ${w} + 0.400" | bc`
Lsbend=10.0
filename="result_w${w}p${p}d${d}"
echo -n > ${filename}.dat
echo "\t1.530\t1.540\t1.550\t1.560\t1.570" >> ${filename}.dat
L1=14.4
L2=8.0
for L in `seq 0.0 0.1 25.0`
do  
    echo -n "${L}" >> ${filename}.dat

    Ldel=`echo "scale=3;25.0 - ${L}" | bc`
    ./make_cfg ${dz} ${g2} ${w} 10 sbend ${Lsbend} ${g1} 3 dc ${L1} 8 dc ${L2} 7 sbend ${Lsbend} ${g2} 3 dc ${L} 9 dc ${Ldel} 7 sbend ${Lsbend} ${g1} 3 dc ${L2} 7 dc ${L1} 9 sbend ${Lsbend} ${g2} 3

    for wl in `seq 1.530 0.010 1.570`
    do
        ~/msh/CMT_sbend/cmt -pcm9 -wl ${wl} -zw
        echo -n "\t`tail -n1 output | head -n1 | cut -f3 `" >> ${filename}.dat
    done
    echo >> ${filename}.dat
done
filename="w${w}p${p}d${d}"
./forExcel ${filename}