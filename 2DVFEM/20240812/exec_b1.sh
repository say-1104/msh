#!/bin/sh
export OMP_NUM_THREADS=24
export MKL_NUM_THREADS=24
# Usage: edge [file] [options]
#Guided mode solver
#  -d R    : PML depth
#  -p RRRR : PML position (Xmin, Ymin, Xmax, Ymax)
#  -w R    : Operating wavelength
#  -s      : Use Sellmeir's formula
#  -a R    : Assumed effective index
#  -r R    : Bending radius (micro meter)
#  -t R    : Attenuation parameter (tanD)
#  -l      : Use linear element
#  -b      : Output field data for BPM
#  -g      : Output field data for GiD
#  -o R    : Select output data format
#  -m      : Input refractive index from .mat file
#  -detail : Show details
#  -h      : Print this help and exist


#32行までは物理的な係数等なので今はとりあえず気にしなくていい
a0=0.696750
a1=0.408218
a2=0.890815
l0=0.069066
l1=0.115662
l2=9.900559
A=0.939816
B=0.00810461
epsilon=11.6858
lamda=1.1071 
nair=1.0
A0=3.7296
A1=0.5424
A2=0.1651
C0=5.1563
C1=1.2388
C2=0.3387

rm _results.dat
for wavelength in `seq $1 $2 $3`
#for wavelength in 1.550
do
    #wavelengthは計算する波長[um]
    rm -R wave_$wavelength
    mkdir wave_$wavelength
    rm -r b1/${wavelength}/${4}
    mkdir b1/${wavelength}/${4}
    sed -i -e "2d" wire.v2.msh
    sed -i -e "2i $wavelength" wire.v2.msh
    #50行まででシリコンシリカの屈折率計算
    pi=`echo "scale=10; 4*a(1)" | bc -l`
    k=`echo "scale=6; 2.0*$pi/$wavelength" | bc -l`
    nsi=`echo "scale=10; sqrt($epsilon+$A/($wavelength*$wavelength)+($B*$lamda*$lamda)/($wavelength*$wavelength-$lamda*$lamda))" | bc -l`
    eps0=`echo "scale=10; $a0*$wavelength*$wavelength/($wavelength*$wavelength-$l0*$l0)" | bc -l`
    eps1=`echo "scale=10; $a1*$wavelength*$wavelength/($wavelength*$wavelength-$l1*$l1)" | bc -l`
    eps2=`echo "scale=10; $a2*$wavelength*$wavelength/($wavelength*$wavelength-$l2*$l2)" | bc -l`
    nsillica=`echo "scale=10; sqrt($eps0+$eps1+$eps2+1)" | bc -l`
    echo $nsi
    echo $nsillica
    napcm=`echo "scale=10; $A0-$wavelength*$A1+$wavelength*$wavelength*$A2" | bc -l`
    ncpcm=`echo "scale=10; $C0-$wavelength*$C1+$wavelength*$wavelength*$C2" | bc -l`
    echo $napcm
    echo $ncpcm

    sed -i -e "6d" wire.v2.msh
    sed -i -e "6i $nsillica  0.00000000" wire.v2.msh
    sed -i -e "7d" wire.v2.msh
    sed -i -e "7i $nair  0.00000000" wire.v2.msh
    sed -i -e "8d" wire.v2.msh
    sed -i -e "8i $nsi  0.00000000" wire.v2.msh
    sed -i -e "9d" wire.v2.msh
    sed -i -e "9i $nair  0.00000000" wire.v2.msh
    sed -i -e "10d" wire.v2.msh
    sed -i -e "10i $nair  0.00000000" wire.v2.msh

    /home/okazaki/Solver/2D_VFEM/edge `ls *v2.msh | sed -e "s/\.msh//g"` -o 0 -w $wavelength -d 1.0 -F 70
    mv `ls *.slv | sed -e "s/\.msh//g"` wave_$wavelength/
    cp _results.dat b1/${wavelength}/${4}/
    mv _results.dat wave_$wavelength/
    mv loss wave_$wavelength/
    cd wave_$wavelength/
    #67行目はとりあえず理解しなく良い
    #ls | grep slv | grep -v "edge" | grep -v "fravia" | grep -v "flavia" | sed -e "s/.slv//g" | sort | runfiles slv2gif  $
    cd ..
done


