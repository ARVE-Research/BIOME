#!/usr/bin/env bash

# this script plots some daily meteorological variables using the Generic Mapping Tools (Classic Mode)

gmt gmtset PS_MEDIA a3

infile=${1}
output=${infile%%.*}.ps

bounds=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg\n",$1,$2,$3,$4,$5}' $infile | gmt gmtinfo -I1/5/5/20 -C) )

t0=${bounds[0]}
t1=${bounds[1]}

tmin=${bounds[4]}
tmax=${bounds[3]}

pmax=${bounds[7]}

echo $t0/$t1/$tmin/$tmax

# ---
# temperature

gmt psbasemap -R$t0/$t1/$tmin/$tmax -JX20/10 -Bpx1o -Bpya10f2 -BWSn -P -K > $output

# night
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$4}' $infile | gmt psxy -R -J -Wthin,blue -O -P -K >> $output

# day
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$3}' $infile | gmt psxy -R -J -Wthin,red -O -P -K >> $output

# ---
# precipitation

gmt psbasemap  -R$t0/$t1/0/$pmax -JX20/10 -Bpya -BE -O -P -K >> $output

awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$5}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,dodgerblue -O -P >> $output

# ---

gmt psconvert -A -Tf -Z $output
