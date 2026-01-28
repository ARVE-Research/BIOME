#!/usr/bin/env bash
# this script plots radiation and PET variables using the Generic Mapping Tools (Classic Mode)

gmt gmtset PS_MEDIA a3
gmt gmtset FORMAT_DATE_MAP o
gmt gmtset FORMAT_TIME_PRIMARY_MAP Character
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black

infile=${1}
output=${place_clean}_radiation.ps

# Get lat/lon from .nc file
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

# get bounds - need columns 18 (swrad), 19 (lw_rad), 20 (lw_rad2), 17 (dpet)
bounds=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg %lg\n",$1,$2,$18,$19,$20,$17}' $infile | gmt gmtinfo -I1/10/1/1 -C) )

t0=${bounds[0]}
t1=${bounds[1]}
swrad_max=${bounds[3]}
lw_min=$(echo "${bounds[4]} ${bounds[6]}" | awk '{print ($1 < $2) ? $1 : $2}')
lw_max=$(echo "${bounds[5]} ${bounds[7]}" | awk '{print ($1 > $2) ? $1 : $2}')
pet_max=${bounds[9]}

# round up swrad_max to nearest 50
swrad_max=$(echo "scale=0; ((${swrad_max%.*}/50)+1)*50" | bc)

# round lw_min down and lw_max up to nearest 10
lw_min=$(echo "scale=0; ((${lw_min%.*}/10)-1)*10" | bc)
lw_max=$(echo "scale=0; ((${lw_max%.*}/10)+1)*10" | bc)

# round up pet_max to nearest 5
pet_max=$(echo "scale=0; ((${pet_max%.*}/5)+1)*5" | bc)

echo "Time range: $t0 to $t1"
echo "Shortwave max: $swrad_max W m-2"
echo "Longwave range: $lw_min to $lw_max W m-2"
echo "PET max: $pet_max mm"

# ---
# GRAPH 1: Net Longwave Radiation (W m-2) - two methods
# ---

gmt psbasemap -R$t0/$t1/$lw_min/$lw_max -JX19/7 -Bpxa1O -Bpya20f10+l"Net LW (W m@+-2@+)" -BWSen+t"Net Longwave Radiation - ${place} (${lon}, ${lat})" -X4 -Y26 -P -K > $output

# lw_rad line - Sandoval/SPLASH 2.0 method (column 19)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$19}' $infile | gmt psxy -R -J -Wthick,dodgerblue -O -P -K >> $output

# lw_rad2 line - Josey/Kimball dewpoint method (column 20)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$20}' $infile | gmt psxy -R -J -Wthick,firebrick -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w5c+o0.2c -F+p0.5p+gwhite -O -P -K << EOF >> $output
S 0.2c - 0.5c - thick,dodgerblue 0.6c LW Sandoval (SPLASH 2.0)
S 0.2c - 0.5c - thick,firebrick 0.6c LW Josey (dewpoint)
EOF

# ---
# GRAPH 2: Shortwave Radiation (W m-2)
# ---

gmt psbasemap -R$t0/$t1/0/$swrad_max -JX19/7 -Bpxa1O -Bpya50f10+l"SW (W m@+-2@+)" -BWSen+t"Downwelling Shortwave Radiation" -Y-10 -P -O -K >> $output

# swrad line (column 18)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$18}' $infile | gmt psxy -R -J -Wthick,orange -O -P -K >> $output

# ---
# GRAPH 3: Daily Potential Evapotranspiration (mm)
# ---

gmt psbasemap -R$t0/$t1/0/$pet_max -JX19/7 -Bpxa1O -Bpya5f1+l"PET (mm d@+-1@+)" -BWSen+t"Daily Potential Evapotranspiration" -Y-10 -P -O -K >> $output

# AET bars (column 14)
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$14}' $infile | gmt psxy -R -J -Sb0.02 -Ggold2@30 -Wthin,lightblue -O -P -K >> $output

# dpet line (column 17)
# awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$17}' $infile | gmt psxy -R -J -Sb0.02 -Gvioletred1@30 -Wthin,violetred1 -O -P >> $output
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$17}' $infile | gmt psxy -R -J -Wthick,violetred1 -O -P >> $output

# ---
# Convert to PDF
# ---

gmt psconvert -A+m0.5c -Tf -Z $output

echo "Plotting complete. PDF saved as ${place_clean}_radiation.pdf"
