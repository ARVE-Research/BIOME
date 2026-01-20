#!/usr/bin/env bash
# this script plots some daily meteorological variables using the Generic Mapping Tools (Classic Mode)
gmt gmtset PS_MEDIA a3
gmt gmtset FORMAT_DATE_MAP o
gmt gmtset FORMAT_TIME_PRIMARY_MAP Character
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black

infile=${1}
output=${place_clean}_weather.ps

# Get lat/lon from netCDF file
lon=$(gmt convert singlepixel.nc?lon | tail -1)
lat=$(gmt convert singlepixel.nc?lat | tail -1)

# Reverse geocode to get place name
place=$(curl -s "https://nominatim.openstreetmap.org/reverse?format=json&lat=${lat}&lon=${lon}&zoom=10" | \
        grep -o '"display_name":"[^"]*' | \
        cut -d'"' -f4 | \
        cut -d',' -f1-2)

if [ -z "$place" ]; then
    place="Unknown location"
fi

echo "Plotting data for: ${place} (Lon: ${lon}, Lat: ${lat})"

bounds=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg\n",$1,$2,$3,$4,$5}' $infile | gmt gmtinfo -I1/5/5/20 -C) )
t0=${bounds[0]}
t1=${bounds[1]}
tmin=${bounds[4]}
tmax=${bounds[3]}
pmax=${bounds[7]}
echo $t0/$t1/$tmin/$tmax

# ---
# GRAPH 1
# temperature
gmt psbasemap -R$t0/$t1/$tmin/$tmax -JX19/10 -Bpxa1O -Bpya10f2+l"Temperature (C)" -BWSn+t"Daily Min/Max Temperatures and Precipitation - ${place} (${lon}, ${lat})" -X4 -Y30 -P -K > $output
# night
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$4}' $infile | gmt psxy -R -J -Wthin,blue -O -P -K >> $output
# day
awk '{printf "2021-%02i-%02iT %lg\n",$1,$2,$3}' $infile | gmt psxy -R -J -Wthin,red -O -P -K >> $output 
# ---
# precipitation
gmt psbasemap  -R$t0/$t1/0/$pmax -JX19/10 -Bpya+l"Precipitation (mm)" -BE -O -P -K >> $output
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$5}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,dodgerblue -O -P -K >> $output

# ---
# GRAPH 2
# snow water equivalent (mm) 
gmt psbasemap  -R$t0/$t1/0/450 -JX19/10 -Bpxa1O -Bpya+l"SWE (mm)" -BWSn+t"Snow Water Equivalent and Snow Cover Fraction" -Y-12 -O -P -K >> $output
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$8}' $infile | gmt psxy -R -J -Sb0.02 -Wthin,slategray1 -O -P -K >> $output
# fsnow fractional snow cover (0–1)
gmt psbasemap -R$t0/$t1/0/1.1 -JX19/10 -Bpya0.1f0.1+l"Snow Cover Fraction" -BE -O -P -K >> "$output"
awk '{printf "2021-%02i-%02iT6:00:00 %lg\n",$1,$2,$9}' $infile | gmt psxy -R -J -Wthick,violetred -O -P -K >> "$output"

# GRAPH 3
#snow melt + accumulation (same panel, transparent overlap)
gmt psbasemap  -R$t0/$t1/0/155 -JX19/10 -Bpxa1O -Bpya+l"Snow (mm)" -BWSen+t"Daily Snow Accumulation and Melt" -Y-12 -O -P -K >> $output
# melt bars (transparent coral)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$7}' $infile | gmt psxy -R -J -Sb0.02 -Glightcoral@50 -Wthin,lightcoral -O -P -K >> $output
# accumulation bars (transparent blue) - same time so they overlap
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$6}' $infile | gmt psxy -R -J -Sb0.02 -Gdodgerblue@50 -Wthin,dodgerblue -O -P -K >> $output
 
# Bsw albedo (0–1)
gmt psbasemap -R$t0/$t1/0/1 -JX19/10 -Bpya0.1f0.1+l"Albedo" -BE -O -P -K >> "$output"
awk '{printf "2021-%02i-%02iT6:00:00 %lg\n",$1,$2,$11}' $infile | gmt psxy -R -J -Wthin,orange -O -P >> "$output"  # last layer, no -K
# ---
 gmt psconvert -A+m0.5c -Tf -Z $output
 
echo "Plotting complete. PDF saved as ${place}weather.pdf"