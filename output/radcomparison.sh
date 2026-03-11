#!/usr/bin/env bash

gmt gmtset PS_MEDIA a3
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black
gmt gmtset TIME_SYSTEM UNIX
gmt gmtset TIME_EPOCH 2000-01-01T00:00:00
gmt gmtset TIME_UNIT d

infile_model=${1}
infile_obs=ricefarm.txt
output=radiation_comparison.ps

# ============================================================
# 1. Observed (top graph)
# ============================================================

# get bounds
lw_min=$(awk '{print $6}' $infile_obs | sort -n | head -1)
lw_max=$(awk '{print $5}' $infile_obs | sort -n | tail -1)

lw_min=$(echo "scale=0; ((${lw_min%.*}/50)-1)*50" | bc)
lw_max=$(echo "scale=0; ((${lw_max%.*}/50)+1)*50" | bc)

t0=$(awk 'NR==1{print $1}' $infile_obs)
t1=$(awk 'END{print $1}' $infile_obs)

gmt psbasemap -R1/365/$lw_min/$lw_max \
  -JX19t/10 \
  -Bpxa1O -Bpya50f10+l"Radiation (W m@+-2@+)" \
  -BWSen+t"Observed - Glenn County Rice Farm, FLUXNET Station US-RGo, 2021-2025" \
  -X4 -Y16 -P -K > $output
  
# SW down (col 2)
awk '{print $1,$2}' $infile_obs | gmt psxy -R -J -Wthick,orange -O -K >> $output

# LW down (col 4)
awk '{print $1,$4}' $infile_obs | gmt psxy -R -J -Wthick,darkblue -O -K >> $output

# LW up (col 5)
awk '{print $1,$5}' $infile_obs | gmt psxy -R -J -Wthick,deepskyblue -O -K >> $output

# Net LW (col 6)
awk '{print $1,$6}' $infile_obs | gmt psxy -R -J -Wthick,purple -O -K >> $output

# Net radiation (col 7)
awk '{print $1,$7}' $infile_obs | gmt psxy -R -J -Wthick,green -O -K >> $output

# Legend
gmt pslegend -R -J -DjTL+w3.2c+o0.2c -F+p0.5p+gwhite -O -K << EOF >> $output
S 0.2c - 0.5c - thick,orange 0.6c SW downwelling
S 0.2c - 0.5c - thick,darkblue 0.6c LW downwelling
S 0.2c - 0.5c - thick,deepskyblue 0.6c LW upwelling
S 0.2c - 0.5c - thick,purple 0.6c Net LW
S 0.2c - 0.5c - thick,green 0.6c Net rad
EOF

# ============================================================
# 2. Modeled (bottom graph)
# ============================================================

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

# get bounds
bounds1=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg %lg %lg %lg\n",$1,$2,-$22,-$24,$18,$18-$22,$31,$32}' $infile_model | gmt gmtinfo -I1/10 -C) )
t0=${bounds1[0]}
t1=${bounds1[1]}
lw_min=$(echo "${bounds1[2]} ${bounds1[4]}" | awk '{print ($1 < $2) ? $1 : $2}')
lw_max=$(echo "${bounds1[3]} ${bounds1[5]} ${bounds1[7]} ${bounds1[9]} ${bounds1[11]} ${bounds1[13]}" | awk '{m=$1; if($2>m)m=$2; if($3>m)m=$3; if($4>m)m=$4; if($5>m)m=$5; if($6>m)m=$6; print m}')
lw_min=$(echo "scale=0; ((${lw_min%.*}/50)-1)*50" | bc)
lw_max=$(echo "scale=0; ((${lw_max%.*}/50)+1)*50" | bc)

echo "Time range: $t0 to $t1"
echo "Longwave/Shortwave range: $lw_min to $lw_max W m-2"

# psbasemap - append to same file, shift down
gmt psbasemap -R$t0/$t1/$lw_min/$lw_max -JX19/10 \
  -Bpxa1O -Bpya50f10+l"Radiation (W m@+-2@+)" \
  -BWSen+t"Single Pixel Model Output (${lon}, ${lat})" \
  -Y-12 -O -P -K >> $output

# swrad (column 18) - orange
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$18}' $infile_model | gmt psxy -R -J -Wthick,orange -O -P -K >> $output

# lwnight (column 24) - purple, negative
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,-$24}' $infile_model | gmt psxy -R -J -Wthick,purple -O -P -K >> $output

# net rad (sw - lw daytime) - green
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$18-$22}' $infile_model | gmt psxy -R -J -Wthick,green -O -P -K >> $output

# lw_down (column 31) - darkblue
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$31}' $infile_model | gmt psxy -R -J -Wthick,darkblue -O -P -K >> $output

# lw_up (column 32) - deepskyblue
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$32}' $infile_model | gmt psxy -R -J -Wthick,deepskyblue -O -P >> $output

# ============================================================
# Convert to PDF
# ============================================================

gmt psconvert -A+m0.5c -Tf -Z $output

echo "Plot complete: ${output%.ps}.pdf"