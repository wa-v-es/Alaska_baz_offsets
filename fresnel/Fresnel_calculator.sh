#! /bin/bash

dtotal=98 #Degrees
dpoint=49 #Degrees
ddepth=2850
d1_km=`echo $dpoint | awk '{print $1*111.1}'`
d2_km=`echo $dtotal $dpoint | awk '{print ($1-$2)*111.1}'`
freq=1 #Hz
vel=13.6 #km.s-1 #P-wave at CMB
vel=7.26 #km.s-1 #S-wave at CMB
#vel=4.41 #km.s-1 #S-wave at 200 km depth (prem)

lambda=`echo $freq $vel | awk '{print $2/$1}'`
n=1
k
fres_rad=`echo $n $lambda $d1_km $d2_km | awk '{print sqrt($1*$2* (($3*$4)/($3+$4)) )}'`
echo dtotal $dtotal dpoint $dpoint
echo d1_km $d1_km d2_km $d2_km
echo freq $freq vel $vel lambda $lambda
echo ----
echo fres_rad $fres_rad
echo $n $lambda $d1_km $d2_km | awk '{print sqrt($1*$2* (($3*$4)/($3+$4)) )}'
