#!/usr/bin/env bash
# this script plots radiation and PET variables using the Generic Mapping Tools (Classic Mode)
# 
# Column reference from tsoutputmod.f90:
#  1 - month
#  2 - day
#  3 - daytime temperature (C)
#  4 - nighttime temperature (C)
#  5 - precipitation (mm)
#  6 - snowfall (mm)
#  7 - snowmelt (mm)
#  8 - snow water equivalent (mm)
#  9 - fractional snow cover (0-1)
# 10 - snow albedo
# 11 - shortwave radiation albedo
# 12 - AET/PET ratio (0-1)
# 13 - dewpoint temperature (C)
# 14 - actual evapotranspiration (mm)
# 15 - soil water content (mm)
# 16 - relative saturation w/whc (0-1)
# 17 - daily PET (mm)
# 18 - total surface shortwave radiation (W m-2)
# 19 - net longwave, Sandoval method (W m-2)
# 20 - sunshine fraction (0-1)
# 21 - daytime temperature (C)
# 22 - daytime longwave radiation (W m-2)
# 23 - nighttime temperature (C)
# 24 - nighttime longwave radiation (W m-2)
# 25 - hour angle of sunset (hours)
# 26 - hour angle of net radiation crossover (hours)
# 27 - positive (daytime) net radiation (J m-2 d-1)
# 28 - negative (nighttime) net radiation (J m-2 d-1)

gmt gmtset PS_MEDIA a3
gmt gmtset FORMAT_DATE_MAP o
gmt gmtset FORMAT_TIME_PRIMARY_MAP Character
gmt gmtset FONT_ANNOT_PRIMARY 9p,Helvetica,black
gmt gmtset FONT_LABEL 11p,Helvetica-Bold,black
gmt gmtset FONT_TITLE 13p,Helvetica-Bold,black

infile=${1}
output=${infile%.txt}_radiation.ps

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

place_clean=$(echo "$place" | tr ' ,' '_' | tr -d "'")

echo "Plotting data for: ${place} (Lon: ${lon}, Lat: ${lat})"

# get bounds for Graph 1: lwday (col 22), lwnight (col 24), swrad (col 18)
bounds1=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg %lg\n",$1,$2,$22,$24,$18}' $infile | gmt gmtinfo -I1/10 -C) )

t0=${bounds1[0]}
t1=${bounds1[1]}
lw_min=$(echo "${bounds1[2]} ${bounds1[4]}" | awk '{print ($1 < $2) ? $1 : $2}')
lw_max=$(echo "${bounds1[3]} ${bounds1[5]} ${bounds1[7]}" | awk '{m=$1; if($2>m)m=$2; if($3>m)m=$3; print m}')

# round lw_min down and lw_max up to nearest 50
lw_min=$(echo "scale=0; ((${lw_min%.*}/50)-1)*50" | bc)
lw_max=$(echo "scale=0; ((${lw_max%.*}/50)+1)*50" | bc)

# get bounds for Graph 2: hour_sw (col 25), hour_net (col 26)
bounds2=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg\n",$1,$2,$25,$26}' $infile | gmt gmtinfo -I1/1 -C) )

hour_min=0
hour_max=14

# get bounds for Graph 3: HNpos (col 27), HNneg (col 28) - divide by 1e6 to get MJ m-2 d-1
bounds3=( $(awk '{printf "2021-%02i-%02iT12:00:00 %lg %lg\n",$1,$2,$27/1e6,$28/1e6}' $infile | gmt gmtinfo -I1/1 -C) )

hn_min=$(echo "${bounds3[2]} ${bounds3[4]}" | awk '{print ($1 < $2) ? $1 : $2}')
hn_max=$(echo "${bounds3[3]} ${bounds3[5]}" | awk '{print ($1 > $2) ? $1 : $2}')

# round hn_min down and hn_max up to nearest 5
hn_min=$(echo "scale=0; ((${hn_min%.*}/5)-1)*5" | bc)
hn_max=$(echo "scale=0; ((${hn_max%.*}/5)+1)*5" | bc)

echo "Time range: $t0 to $t1"
echo "Longwave/Shortwave range: $lw_min to $lw_max W m-2"
echo "Hour range: $hour_min to $hour_max hours"
echo "HN range: $hn_min to $hn_max MJ m-2 d-1"

# ---
# GRAPH 1: Longwave (day/night) and Shortwave Radiation (W m-2)
# ---

gmt psbasemap -R$t0/$t1/$lw_min/$lw_max -JX19/7 -Bpxa1O -Bpya50f10+l"Radiation (W m@+-2@+)" -BWSen+t"Radiation Components - ${place} (${lon}, ${lat})" -X4 -Y26 -P -K > $output

# swrad line (column 18) - orange
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$18}' $infile | gmt psxy -R -J -Wthick,orange -O -P -K >> $output

# lwday line (column 22) - firebrick
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$22}' $infile | gmt psxy -R -J -Wthick,firebrick -O -P -K >> $output

# lwnight line (column 24) - navy
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$24}' $infile | gmt psxy -R -J -Wthick,navy -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w5c+o0.2c -F+p0.5p+gwhite -O -P -K << EOF >> $output
S 0.2c - 0.5c - thick,orange 0.6c SW downwelling
S 0.2c - 0.5c - thick,firebrick 0.6c LW daytime
S 0.2c - 0.5c - thick,navy 0.6c LW nighttime
EOF

# ---
# GRAPH 2: Hour angles (sunset and net radiation crossover)
# ---

gmt psbasemap -R$t0/$t1/$hour_min/$hour_max -JX19/7 -Bpxa1O -Bpya2f1+l"Hours" -BWSen+t"Hour Angles" -Y-10 -P -O -K >> $output

# hour_sw line (column 25) - gold
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$25}' $infile | gmt psxy -R -J -Wthick,gold2 -O -P -K >> $output

# hour_net line (column 26) - purple
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$26}' $infile | gmt psxy -R -J -Wthick,purple -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w5c+o0.2c -F+p0.5p+gwhite -O -P -K << EOF >> $output
S 0.2c - 0.5c - thick,gold2 0.6c Hour of sunset
S 0.2c - 0.5c - thick,purple 0.6c Hour of net rad crossover
EOF

# ---
# GRAPH 3: Net Radiation (MJ m-2 d-1)
# ---

gmt psbasemap -R$t0/$t1/$hn_min/$hn_max -JX19/7 -Bpxa1O -Bpya5f1+l"Net Radiation (MJ m@+-2@+ d@+-1@+)" -BWSen+t"Accumulated Net Radiation" -Y-10 -P -O -K >> $output

# HNpos line (column 27) - divide by 1e6 for MJ
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$27/1e6}' $infile | gmt psxy -R -J -Wthick,springgreen3 -O -P -K >> $output

# HNneg line (column 28) - divide by 1e6 for MJ
awk '{printf "2021-%02i-%02iT06:00:00 %lg\n",$1,$2,$28/1e6}' $infile | gmt psxy -R -J -Wthick,violetred3 -O -P -K >> $output

# legend
gmt pslegend -R -J -DjTR+w4c+o0.2c -F+p0.5p+gwhite -O -P << EOF >> $output
S 0.2c - 0.5c - thick,springgreen3 0.6c HNpos (daytime)
S 0.2c - 0.5c - thick,violetred3 0.6c HNneg (nighttime)
EOF

# ---
# Convert to PDF
# ---

gmt psconvert -A+m0.5c -Tf -Z $output

echo "Plotting complete. PDF saved as ${output%.ps}.pdf"