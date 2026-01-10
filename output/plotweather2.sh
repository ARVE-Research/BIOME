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

gmt psbasemap -R$t0/$t1/$tmin/$tmax -JX20/10 -Bpx1o -Bpya10f2 -BWSn -Y15 -P -K > $output

# night
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$4}' $infile | gmt psxy -R -J -Wthin,blue -O -P -K >> $output

# day
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$3}' $infile | gmt psxy -R -J -Wthin,red -O -P -K >> $output

# ---
# precipitation

gmt psbasemap  -R$t0/$t1/0/$pmax -JX20/10 -Bpya -BE -O -P -K >> $output

awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$5}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,dodgerblue -O -P -K >> $output

# ---
# LOWER GRAPH
# snow water equivalent (mm) 

gmt psbasemap  -R$t0/$t1/0/400 -JX20/10 -Bpya -BESn -Y-11 -O -P -K >> $output

awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$8}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,slategray1 -O -P -K >> $output

# snow melt (mm)

gmt psbasemap  -R$t0/$t1/0/400 -JX20/10 -Bpya -BESn -O -P -K >> $output

awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$7}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,lightred -O -P -K >> $output

# daily snow accumulation (mm)
gmt psbasemap  -R$t0/$t1/0/400 -JX20/10 -Bpya -BESn -O -P -K >> $output

awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$6}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,dodgerblue -O -P -K >> $output
 
# fsnow fractional snow cover (0–1)
gmt psbasemap -R$t0/$t1/0/1 -JX20/10 -Bpya0.1f0.1 -BW -O -P -K >> "$output"

awk '{printf "2021-%02i-%02iT6:00:00 %lg\n",$1,$2,$9}' $infile | gmt psxy -R -J -Wthick,slateblue2 -O -P >> "$output"  # last layer, no -K

# ---

 gmt psconvert -A -Tf -Z $output

echo "Plotting complete. PDF saved"